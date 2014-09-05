/**
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#define _CRT_SECURE_NO_WARNINGS

/*#include "config.h"*/
#include "exception.h"
#include "jni_helper.h"

#include <stdio.h>
#include <string.h>

#ifdef WIN32
#include "uthash.h"
#endif

static void hdfsThreadDestructor(void *v);

/** Pthreads thread-local storage for each library thread. */
typedef struct HdfsTls {
    JNIEnv *env;
} HDFSTLS;

static JavaVM * hdfs_JVM = NULL;
static short hdfs_InitLib = 0;

#ifndef WIN32

static void* hdfs_dl_handle = NULL;

/*static pthread_rwlock_t hdfs_HashLock;*/

static pthread_mutex_t hdfs_HashMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t hdfs_JvmMutex = PTHREAD_MUTEX_INITIALIZER;

static int hdfs_hashTableInited = 0;

/* SetUp Control variables */
static pthread_once_t hdfs_threadInit_Once = PTHREAD_ONCE_INIT;
static pthread_once_t hdfs_hashTable_Once = PTHREAD_ONCE_INIT;

/** nonzero if we succeeded in initializing gTlsKey */
static int hdfs_gTlsKeyInitialized = 0;

/** Key that allows us to retrieve thread-local storage */
static pthread_key_t hdfs_gTlsKey = 0;

#define LOCK_HASH_TABLE() pthread_mutex_lock(&hdfs_HashMutex)
#define UNLOCK_HASH_TABLE() pthread_mutex_unlock(&hdfs_HashMutex)

#else

typedef struct {
    const char * key;
    void * cls;
    UT_hash_handle hh;
} HASHITEM, * PHASHITEM;

PHASHITEM hdfs_HashTls = NULL;

static DWORD hdfs_dwTlsIndex1 = 0;
static HINSTANCE hdfs_hinstLib = 0;

static HANDLE hdfs_JvmMutex = 0;
static PSRWLOCK hdfs_HashLock = NULL;

#endif

#ifdef WIN32
#define LOCK_JVM_MUTEX() \
dwWaitResult = WaitForSingleObject(hdfs_JvmMutex,INFINITE) 
#else    
#define LOCK_JVM_MUTEX() \
pthread_mutex_lock(&hdfs_JvmMutex)
#endif

#ifdef WIN32
#define UNLOCK_JVM_MUTEX() \
    ReleaseMutex(hdfs_JvmMutex) 
#else    
#define UNLOCK_JVM_MUTEX() \
    pthread_mutex_unlock(&hdfs_JvmMutex)
#endif

/** The Native return types that methods could return */
#define JVOID         'V'
#define JOBJECT       'L'
#define JARRAYOBJECT  '['
#define JBOOLEAN      'Z'
#define JBYTE         'B'
#define JCHAR         'C'
#define JSHORT        'S'
#define JINT          'I'
#define JLONG         'J'
#define JFLOAT        'F'
#define JDOUBLE       'D'

/**
 * MAX_HASH_TABLE_ELEM: The maximum no. of entries in the hashtable.
 * It's set to 4096 to account for (classNames + No. of threads)
 */
#define MAX_HASH_TABLE_ELEM 4096

#ifndef WIN32

/* Invoked By "pthread_once" to create the thread key */
static void Make_Thread_Key()
{
    int ret = pthread_key_create(&hdfs_gTlsKey, hdfsThreadDestructor);
    if (ret) {
        fprintf(stderr, "getJNIEnv: pthread_key_create failed with "
            "error %d\n", ret);
        return;
    }
    hdfs_gTlsKeyInitialized = 1;    
}

/* Invoked By "pthread_once" to create the Hash Table */
static void hashTableInit ()
{

 /* if ( pthread_rwlock_init (&hdfs_HashLock, NULL) != 0 ) {
        fprintf (stderr, "can't create rwlock for hash table");
        return;
    } */

    if ( hcreate(MAX_HASH_TABLE_ELEM) == 0 ) {
        fprintf ( stderr, "error creating hashtable, <%d>: %s\n",
                  errno, strerror(errno) );
        return;
    }
    hdfs_hashTableInited = 1; 
}

static int insertEntryIntoTable ( const char *key, void *data )
{

    ENTRY e, *ep = NULL;
    if (key == NULL || data == NULL) {
        return 0;
    }

    pthread_once ( &hdfs_hashTable_Once, hashTableInit );    
    if ( !hdfs_hashTableInited ) {
      return -1;
    }

    e.data = data;
    e.key = (char*) key;
    
 /* if ( pthread_rwlock_wrlock(&hdfs_HashLock) != 0 ) {
        fprintf (stderr, "can't get hash table wlock");
        return -1;
    } */

    LOCK_HASH_TABLE(); 
    ep = hsearch(e, ENTER);
    UNLOCK_HASH_TABLE(); 
    
    /*pthread_rwlock_unlock(&hdfs_HashLock); */
    
    if (ep == NULL) {
        fprintf(stderr, "warn adding key (%s) to hash table, <%d>: %s\n",
                key, errno, strerror(errno));
    }

    return 0;

}

static void* searchEntryFromTable ( const char *key )
{

    ENTRY e, *ep = NULL;
    if (key == NULL) {
        return NULL;
    }

    pthread_once ( &hdfs_hashTable_Once, hashTableInit );
    if ( !hdfs_hashTableInited ) {
      return NULL;
    }

    e.key = (char*)key;
    
 /* if ( pthread_rwlock_rdlock(&hdfs_HashLock) != 0 ) {
        fprintf (stderr, "can't get hash table rdlock");
        return NULL;
    }*/

    LOCK_HASH_TABLE(); 
    ep = hsearch(e, FIND);
    UNLOCK_HASH_TABLE(); 
    
    /*pthread_rwlock_unlock(&hdfs_HashLock);*/
    
    if (ep != NULL) {
        return ep->data;
    }

    return NULL;

}
#else

static int insertEntryIntoTable ( const char *key, void *cls )
{

    PHASHITEM item = NULL;
    DWORD dwWaitResult;
    
    if (key == NULL || cls == NULL) {
        return 0;
    }
               
    item = (PHASHITEM) calloc ( 1, sizeof(HASHITEM) );

    item->key = key;
    item->cls = cls;
    
    AcquireSRWLockExclusive ( hdfs_HashLock );
             
    HASH_ADD_KEYPTR ( hh, hdfs_HashTls, item->key, strlen(item->key), item );

    ReleaseSRWLockExclusive ( hdfs_HashLock );    
        
    return 0;

}

static void* searchEntryFromTable ( const char *key )
{

    PHASHITEM item = NULL;
    DWORD dwWaitResult;
    
    if (key == NULL) {
        return NULL;
    }
            
    AcquireSRWLockShared ( hdfs_HashLock );
    
    if (!hdfs_HashTls) {
    
        // char *mykey = "dummy";    
                        
        ReleaseSRWLockShared ( hdfs_HashLock );
        
        /* add dummy entry for testing */
        // insertEntryIntoTable ( mykey, -1 );        
        // HASH_FIND_STR( hdfs_HashTls, mykey, item);        
        // printf ("Test Hash mykey: %s\n", item->key);
        
        return NULL;
        
    }
        
    HASH_FIND_STR( hdfs_HashTls, key, item);
    
    ReleaseSRWLockShared ( hdfs_HashLock );    
    
    if (item != NULL) {        
        return item->cls;
    }
        
    return NULL;

}

#endif


/**
 * The function that is called whenever a thread with libhdfs thread local data
 * is destroyed.
 *
 * @param v         The thread-local data
 */
static void hdfsThreadDestructor(void *v)
{
    HDFSTLS *tls = v;
    JavaVM *vm;
    JNIEnv *env = tls->env;
    jint ret;

    ret = (*env)->GetJavaVM(env, &vm);
    if (ret) {
        fprintf(stderr, "hdfsThreadDestructor: GetJavaVM failed with "
                "error %d\n", ret);
        (*env)->ExceptionDescribe(env);
    } else {
        (*vm)->DetachCurrentThread(vm);
    }

    free(tls);

}

#ifdef WIN32
BOOL WINAPI DllMain (
    HINSTANCE hinstDLL,  // handle to DLL module
    DWORD fdwReason,     // reason for calling function
    LPVOID lpReserved )  // reserved
{

    HDFSTLS *tls = NULL;
    
    // Perform actions based on the reason for calling.
    switch( fdwReason )
    {
        case DLL_PROCESS_ATTACH:
        
            if ((hdfs_dwTlsIndex1 = TlsAlloc()) == TLS_OUT_OF_INDEXES) 
                return FALSE; 

            hdfs_HashLock = (PSRWLOCK) calloc ( 1, sizeof ( SRWLOCK ) );
            if (!hdfs_HashLock) {
                fprintf (stderr, " Could not allocate Hash Table Lock\n" );
                return FALSE;
            }
            
            InitializeSRWLock ( hdfs_HashLock );
                
            hdfs_JvmMutex = CreateMutex ( 
                NULL,              // default security attributes
                FALSE,             // initially not owned
                NULL );            // unnamed mutex
                
            if (hdfs_JvmMutex == NULL) 
            {
                fprintf (stderr, "Create JVM Mutex error: %d\n", GetLastError() );
                return FALSE;
            }

            fprintf (stderr, "dll attached\n" );
            fprintf (stderr, "dll: tls1=%d\n", hdfs_dwTlsIndex1 ); 
            
            // Initialize once for each new process.
            // Return FALSE to fail DLL load.
            break;

        case DLL_THREAD_ATTACH:
            
            // Do thread-specific initialization.
            fprintf (stderr, "dll: thread attach\n" );
            
            break;

        case DLL_THREAD_DETACH:
            
            // Do thread-specific cleanup.
            fprintf (stderr, "dll: detach thread\n" );
            
            tls = TlsGetValue(hdfs_dwTlsIndex1); 
            if (tls) {
                fprintf (stderr, "dll thread: invoke thread destructor\n" );
                hdfsThreadDestructor(tls);
                TlsSetValue(hdfs_dwTlsIndex1, NULL);
            }
                        
            break;

        case DLL_PROCESS_DETACH:
            
            // Perform any necessary cleanup.
            fprintf (stderr, "dll: detach process\n" );
            
            if (hdfs_HashTls) {
                PHASHITEM item, tmp;
                fprintf (stderr, "dll: clean up hash table\n" );                
                HASH_ITER ( hh, hdfs_HashTls, item, tmp ) {
                    fprintf (stderr, "dll: delete key %s\n", item->key );
                    HASH_DEL ( hdfs_HashTls, item );
                    free ( item );
                }
            }
            
            tls = TlsGetValue(hdfs_dwTlsIndex1);             
            if (tls) {
                fprintf (stderr, "dll: invoke thread destructor\n" );
                hdfsThreadDestructor(tls);
            }
 
            // Release the TLS index.
 
            TlsFree(hdfs_dwTlsIndex1);  
            
            if (hdfs_JvmMutex) CloseHandle(hdfs_JvmMutex);            
            if (hdfs_HashLock) free(hdfs_HashLock);
            
            fprintf (stderr, "dll detached\n" );
            
            break;
    }

    return TRUE;  // Successful DLL_PROCESS_ATTACH.

}

// jelson changes: start
INIT_ONCE g_InitOnce = INIT_ONCE_STATIC_INIT;

BOOL CALLBACK staticLibInit(PINIT_ONCE InitOnce, PVOID Parameter, PVOID *lpContext)
{
	DllMain(NULL, DLL_PROCESS_ATTACH, NULL);
}

void maybePerformStaticLibInit()
{
	InitOnceExecuteOnce(&g_InitOnce, staticLibInit, NULL, NULL);
}
// jelson changes: end


#endif

int hdfsLibInit ( void * parms )
{

    JNIEnv* env = getJNIEnv();
    
    if (!env) return 1;
    
    hdfs_InitLib = 1;
    
    return 0;

}

// test compiler 

void destroyLocalReference(JNIEnv *env, jobject jObject)
{
    if (jObject)
        (*env)->DeleteLocalRef(env, jObject);
}

static jthrowable validateMethodType(JNIEnv *env, MethType methType)
{
    if (methType != STATIC && methType != INSTANCE) {
        return newRuntimeError(env, "validateMethodType(methType=%d): "
            "illegal method type.\n", methType);
    }
    return NULL;
}

jthrowable newJavaStr(JNIEnv *env, const char *str, jstring *out)
{
    jstring jstr;

    if (!str) {
        /* Can't pass NULL to NewStringUTF: the result would be
         * implementation-defined. */
        *out = NULL;
        return NULL;
    }

    jstr = (*env)->NewStringUTF(env, str);
    if (!jstr) {
        /* If NewStringUTF returns NULL, an exception has been thrown,
         * which we need to handle.  Probaly an OOM. */
        return getPendingExceptionAndClear(env);
    }

    *out = jstr;
    return NULL;

}

jthrowable newCStr(JNIEnv *env, jstring jstr, char **out)
{
    const char *tmp;

    if (!jstr) {
        *out = NULL;
        return NULL;
    }

    tmp = (*env)->GetStringUTFChars(env, jstr, NULL);
    if (!tmp) {
        return getPendingExceptionAndClear(env);
    }

    *out = _strdup(tmp);
    (*env)->ReleaseStringUTFChars(env, jstr, tmp);
    return NULL;

}


jthrowable invokeMethod ( JNIEnv *env, jvalue *retval, MethType methType,
                 jobject instObj, const char *className,
                 const char *methName, const char *methSignature, ... )
{

    va_list args;
    jclass cls;
    jmethodID mid = NULL;
    jthrowable jthr;
    const char *str;
    char returnType;

    jthr = validateMethodType(env, methType);
    if (jthr)
        return jthr;

    jthr = globalClassReference(className, env, &cls);
    if (jthr) 
        return jthr;

    jthr = methodIdFromClass(className, methName, methSignature,
                            methType, env, &mid);
    if (jthr) {
        return jthr;
    }
        
    str = methSignature;
    while (*str != ')') str++;
    str++;
    returnType = *str;

//    printf ("Begin Method Invokation:%s ## %s\n", className, methName );
    
    va_start(args, methSignature);
    if (returnType == JOBJECT || returnType == JARRAYOBJECT) {
        jobject jobj = NULL;
        if (methType == STATIC) {
            jobj = (*env)->CallStaticObjectMethodV(env, cls, mid, args);
        }
        else if (methType == INSTANCE) {
            jobj = (*env)->CallObjectMethodV(env, instObj, mid, args);
        }
        retval->l = jobj;
    }
    else if (returnType == JVOID) {
        if (methType == STATIC) {
            (*env)->CallStaticVoidMethodV(env, cls, mid, args);
        }
        else if (methType == INSTANCE) {
            (*env)->CallVoidMethodV(env, instObj, mid, args);
        }
    }
    else if (returnType == JBOOLEAN) {
        jboolean jbool = 0;
        if (methType == STATIC) {
            jbool = (*env)->CallStaticBooleanMethodV(env, cls, mid, args);
        }
        else if (methType == INSTANCE) {
            jbool = (*env)->CallBooleanMethodV(env, instObj, mid, args);
        }
        retval->z = jbool;
    }
    else if (returnType == JSHORT) {
        jshort js = 0;
        if (methType == STATIC) {
            js = (*env)->CallStaticShortMethodV(env, cls, mid, args);
        }
        else if (methType == INSTANCE) {
            js = (*env)->CallShortMethodV(env, instObj, mid, args);
        }
        retval->s = js;
    }
    else if (returnType == JLONG) {
        jlong jl = -1;
        if (methType == STATIC) {
            jl = (*env)->CallStaticLongMethodV(env, cls, mid, args);
        }
        else if (methType == INSTANCE) {
            jl = (*env)->CallLongMethodV(env, instObj, mid, args);
        }
        retval->j = jl;
    }
    else if (returnType == JINT) {
        jint ji = -1;
        if (methType == STATIC) {
            ji = (*env)->CallStaticIntMethodV(env, cls, mid, args);
        }
        else if (methType == INSTANCE) {
            ji = (*env)->CallIntMethodV(env, instObj, mid, args);
        }
        retval->i = ji;
    }
    va_end(args);

//    printf ("End Method Invokation\n");
    
    jthr = (*env)->ExceptionOccurred(env);
    if (jthr) { 
        (*env)->ExceptionClear(env);
        return jthr;
    }

//    printf ("Method success\n");
    
    return NULL;

}

jthrowable constructNewObjectOfClass(JNIEnv *env, jobject *out, const char *className,
                                  const char *ctorSignature, ...)
{

    va_list args;
    jclass cls;
    jmethodID mid = NULL;
    jobject jobj;
    jthrowable jthr;

    jthr = globalClassReference(className, env, &cls);
    if (jthr)
        return jthr;

    jthr = methodIdFromClass(className, "<init>", ctorSignature,
                            INSTANCE, env, &mid);
    if (jthr)
        return jthr;

    va_start(args, ctorSignature);
    jobj = (*env)->NewObjectV(env, cls, mid, args);
    va_end(args);

    if (!jobj)
        return getPendingExceptionAndClear(env);

    *out = jobj;
    return NULL;

}


jthrowable methodIdFromClass(const char *className, const char *methName,
                            const char *methSignature, MethType methType,
                            JNIEnv *env, jmethodID *out)
{

    jclass cls;
    jthrowable jthr;
    jmethodID mid = NULL;

    jthr = validateMethodType(env, methType);
    if (jthr)
        return jthr;
    
    jthr = globalClassReference(className, env, &cls);
    if (jthr)
        return jthr;

    if (methType == STATIC) {
        mid = (*env)->GetStaticMethodID(env, cls, methName, methSignature);
    }
    else if (methType == INSTANCE) {
        mid = (*env)->GetMethodID(env, cls, methName, methSignature);
    }

    if (mid == NULL) {
        fprintf(stderr, "could not find method %s from class %s with "
            "signature %s\n", methName, className, methSignature);
        return getPendingExceptionAndClear(env);
    }

    *out = mid;
    return NULL;

}

jthrowable globalClassReference(const char *className, JNIEnv *env, jclass *out)
{
    jclass clsLocalRef;

    jclass cls = searchEntryFromTable(className);
    if (cls) {
        *out = cls;
        return NULL;
    }

    clsLocalRef = (*env)->FindClass(env,className);
    if (clsLocalRef == NULL) {
        return getPendingExceptionAndClear(env);
    }
    
    cls = (*env)->NewGlobalRef(env, clsLocalRef);
    if (cls == NULL) {
        (*env)->DeleteLocalRef(env, clsLocalRef);
        return getPendingExceptionAndClear(env);
    }

    (*env)->DeleteLocalRef(env, clsLocalRef);
    
    insertEntryIntoTable(className, cls);

    *out = cls;
    return NULL;

}


jthrowable classNameOfObject(jobject jobj, JNIEnv *env, char **name)
{

    jthrowable jthr;
    jclass cls, clsClass = NULL;
    jmethodID mid = NULL;
    jstring str = NULL;
    const char *cstr = NULL;
    char *newstr;

    cls = (*env)->GetObjectClass(env, jobj);
    if (cls == NULL) {
        jthr = getPendingExceptionAndClear(env);
        goto done;
    }

    clsClass = (*env)->FindClass(env, "java/lang/Class");
    if (clsClass == NULL) {
        jthr = getPendingExceptionAndClear(env);
        goto done;
    }

    mid = (*env)->GetMethodID(env, clsClass, "getName", "()Ljava/lang/String;");
    if (mid == NULL) {
        jthr = getPendingExceptionAndClear(env);
        goto done;
    }

    str = (*env)->CallObjectMethod(env, cls, mid);
    if (str == NULL) {
        jthr = getPendingExceptionAndClear(env);
        goto done;
    }

    cstr = (*env)->GetStringUTFChars(env, str, NULL);
    if (!cstr) {
        jthr = getPendingExceptionAndClear(env);
        goto done;
    }

    newstr = _strdup(cstr);
    if (newstr == NULL) {
        jthr = newRuntimeError(env, "classNameOfObject: out of memory");
        goto done;
    }

    *name = newstr;
    jthr = NULL;

done:
    destroyLocalReference(env, cls);
    destroyLocalReference(env, clsClass);
    if (str) {
        if (cstr)
            (*env)->ReleaseStringUTFChars(env, str, cstr);
        (*env)->DeleteLocalRef(env, str);
    }
    return jthr;

}


/**
 * Get the global JNI environemnt.
 *
 * We only have to create the JVM once.  After that, we can use it in
 * every thread.  You must be holding the jvmMutex when you call this
 * function.
 *
 * @return          The JNIEnv on success; error code otherwise
 */
static JNIEnv* getGlobalJNIEnv(void)
{

    const jsize vmBufLength = 1;
    JavaVM* vmBuf[1]; 
    JNIEnv *env = NULL;
    jint rv = 0;
    jint noVMs = 0;
    jthrowable jthr;    
    /*JavaVM *vm = NULL; */
    char *error = NULL;
    
    typedef jint (*FGetVMS) (JavaVM**, const jsize, jint* );
    FGetVMS  fpGetVM = NULL;   

    typedef jint (*FCreateVM) (JavaVM**, void**, JavaVMInitArgs* );
    FCreateVM  fpCreateVM = NULL;

//    printf ( "Get Global JNI\n" );
    
    #ifndef WIN32

    char *JVMPath = getenv("LIBHDFS_JVM_PATH");
    char jvmPath [2000] = "";

    if (JVMPath) {
        strcpy (jvmPath, JVMPath);
    } else {
        strcpy (jvmPath, "libjvm.so");
    }

    if (!hdfs_dl_handle) {    
        hdfs_dl_handle = (void*) dlopen( jvmPath, 0 );
        if (!hdfs_dl_handle) {
            printf( "!!! %s\n", dlerror() );
            return NULL;
        }
    }  
    
    if (hdfs_dl_handle) {
    
        fpGetVM = (FGetVMS) dlsym( hdfs_dl_handle, "JNI_GetCreatedJavaVMs" );
        error = (char*) dlerror();
        if (error != NULL) {
            fprintf(stderr, "!!! %s\n", error );
            return NULL;
        }   
        
        fpCreateVM = (FCreateVM) dlsym( hdfs_dl_handle, "JNI_CreateJavaVM" );
        error = (char*) dlerror();
        if (error != NULL) {
            fprintf(stderr, "!!! %s\n", error );
            return NULL;
        }

    }
        
    #else

    if (hdfs_hinstLib == NULL) {
		wchar_t *env_libhdfs_jvm_path;
		size_t env_libhdfs_jvm_path_size;

		wchar_t jvmPath[2000];
    
        fprintf (stderr, "libhdfs: loading jvm\n" );
        
		_wdupenv_s(&env_libhdfs_jvm_path, &env_libhdfs_jvm_path_size, L"LIBHDFS_JVM_PATH");

        if (env_libhdfs_jvm_path != NULL) {
			fprintf(stderr, "using value from LIBHDFS_JVM_PATH: %S\n", env_libhdfs_jvm_path);
			wcscpy_s (jvmPath, _countof(jvmPath), env_libhdfs_jvm_path);
        } else {
			wchar_t *env_java_home;
			size_t env_java_home_size;

			_wdupenv_s(&env_java_home, &env_java_home_size, L"JAVA_HOME");

			if (env_java_home != NULL) {
				_snwprintf(jvmPath, _countof(jvmPath), L"%s\\jre\\bin\\server\\jvm.dll", env_java_home);
				fprintf(stderr, "Found JAVA_HOME of %S; trying jvm path of '%S'\n", env_java_home, jvmPath);
			} else {
				wcscpy_s (jvmPath, _countof(jvmPath), L"c:\\program files\\java\\jre\\bin\\server\\jvm.dll");
				fprintf(stderr, "LIBHDFS_JVM_PATH and JAVA_HOME both not set; blindingly trying %S\n", jvmPath);
			}
        }
        
        hdfs_hinstLib = LoadLibrary ( jvmPath );
        if (!hdfs_hinstLib) {
                    
            LPVOID lpMsgBuf;    
            DWORD dw = GetLastError(); 

            FormatMessage (
                FORMAT_MESSAGE_ALLOCATE_BUFFER | 
                FORMAT_MESSAGE_FROM_SYSTEM |
                FORMAT_MESSAGE_IGNORE_INSERTS,
                NULL,
                dw,
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                (LPTSTR) &lpMsgBuf,
                0, NULL );
        
            fprintf (stderr, "jvm load failed\n" );
            fprintf (stderr, "Error Code:%d %s\n", dw, lpMsgBuf );

            LocalFree(lpMsgBuf);

            return NULL;
            
        }       
        
    } /* endif load the dll */
        
    if (hdfs_hinstLib)
    {
    
        // fprintf (stderr, "dll: get proc addresses\n" );
        
        fpGetVM = (FGetVMS) GetProcAddress ( hdfs_hinstLib,
                               "JNI_GetCreatedJavaVMs" );

        if (!fpGetVM) {
            fprintf (stderr, "dll: could not get proc 'JNI_GetCreatedJavaVMs'\n" );
            return NULL;
        }   
                               
        fpCreateVM = (FCreateVM) GetProcAddress ( hdfs_hinstLib,
                                   "JNI_CreateJavaVM" );
        if (!fpCreateVM) {
            fprintf (stderr, "dll: could not get proc 'JNI_CreateJavaVM'\n" );
            return NULL;
        }
                                   
    }

    #endif
    
    /*rv = JNI_GetCreatedJavaVMs(vmBuf, vmBufLength, &noVMs);*/
    rv = (fpGetVM) (vmBuf, vmBufLength, &noVMs);
    if (rv != 0) {
        fprintf(stderr, "JNI_GetCreatedJavaVMs failed with error: %d\n", rv);
        return NULL;
    }   
    
    if (noVMs == 0) {  

        /*Get the environment variables for initializing the JVM */
        char *hadoopClassPath = getenv("LIBHDFS_CLASSPATH");
        char *hadoopClassPathVMArg = "-Djava.class.path="; 
        /*char *hadoopClassPathVMArg = "-cp ";*/
        int  optHadoopClassPathLen;
        char *optHadoopClassPath = NULL;

        int  noArgs = 1, cnt = 0;
        char *hadoopJvmArgs = NULL;
        char jvmArgDelims[] = " ";
        char *str, *token, *savePtr;
        JavaVMOption *options = NULL;
        JavaVMInitArgs vm_args;

        if (hadoopClassPath == NULL) {
            fprintf(stderr, "libhdfs: Environment variable LIBHDFS_CLASSPATH not set!\n");
            return NULL;
        }

        optHadoopClassPathLen = strlen(hadoopClassPath) +
          strlen(hadoopClassPathVMArg) + 1;

        optHadoopClassPath = malloc(sizeof(char)*optHadoopClassPathLen);
        /* snprintf ( optHadoopClassPath, optHadoopClassPathLen,
                     "%s%s", hadoopClassPathVMArg, hadoopClassPath );*/

        #if 0        
        savePtr = strdup(hadoopClassPath);
        str     = strdup(hadoopClassPath);
        savePtr[0] = 0;        
        
        cnt = strlen (str);
        token = strtok (str, "\\");
        
        while (token) {                
            strcat (savePtr,token);
            if (strlen(savePtr) == cnt) break;
            strcat (savePtr,"/");
            token = strtok (NULL, "\\");        
        }
       
        sprintf ( optHadoopClassPath,
                   "%s%s", hadoopClassPathVMArg, savePtr );
        
        free (savePtr);
        free (str);        
        #endif
        
        sprintf ( optHadoopClassPath,
                   "%s%s", hadoopClassPathVMArg, hadoopClassPath );        
                           
        /* Determine the # of LIBHDFS_OPTS args */
        hadoopJvmArgs = getenv("LIBHDFS_OPTS");

        if (hadoopJvmArgs != NULL)  {
          hadoopJvmArgs = _strdup(hadoopJvmArgs);
          for (noArgs = 1, str = hadoopJvmArgs; ; noArgs++, str = NULL) {
            /* token = strtok_r(str, jvmArgDelims, &savePtr);*/
            token = strtok (str, jvmArgDelims);
            if (NULL == token) {
              break;
            }
          }
          free(hadoopJvmArgs);
        }

        /* Now that we know the # args, populate the options array */
        options = calloc(noArgs, sizeof(JavaVMOption));
        options[0].optionString = optHadoopClassPath;

        hadoopJvmArgs = getenv("LIBHDFS_OPTS");

        if (hadoopJvmArgs != NULL)  {
          hadoopJvmArgs = _strdup(hadoopJvmArgs);
          for (noArgs = 1, str = hadoopJvmArgs; ; noArgs++, str = NULL) {
            /* token = strtok_r(str, jvmArgDelims, &savePtr);*/
            token = strtok (str, jvmArgDelims);
            if (NULL == token) {
              break;
            }
            options[noArgs].optionString = token;
          }
        }

        /*Create the VM */
        vm_args.version = JNI_VERSION_1_2;
        vm_args.options = options;
        vm_args.nOptions = noArgs;
        vm_args.ignoreUnrecognized = 1;

        /*rv = JNI_CreateJavaVM(&vm, (void**) &env, &vm_args);*/
        rv = (fpCreateVM) (&hdfs_JVM, (void**) &env, &vm_args);

        if (hadoopJvmArgs != NULL)  {
          free(hadoopJvmArgs);
        }

        free(optHadoopClassPath);
        free(options);

        if (rv != 0) {
            fprintf(stderr, "Call to JNI_CreateJavaVM failed "
                    "with error: %d\n", rv);
            return NULL;
        }

        fprintf (stderr, "dll: jvm created\n" );
        
        /*This is not backwards comaptible */
        /*
        jthr = invokeMethod ( env, NULL, STATIC, NULL,
                         "org/apache/hadoop/fs/FileSystem",
                         "loadFileSystems", "()V" );
        if (jthr) {            
            printExceptionAndFree ( env, jthr, PRINT_EXC_ALL,
                                    "loadFileSystems" );
            return NULL;
        }
        
        printf ( "dll: return from GetEnv\n" ); */
        
        return env;

    }

    // fprintf (stderr, "dll: attach current thread \n" );
    
    /*Attach this thread to the VM */       
  /*vm = vmBuf[0];
    rv = (*vm)->AttachCurrentThread(vm, (void**) &env, 0); */

    if (!hdfs_JVM) hdfs_JVM = vmBuf[0];    

    rv = (*hdfs_JVM)->AttachCurrentThread(hdfs_JVM, (void**) &env, 0);
    if (rv != 0) {
        fprintf(stderr, "Call to AttachCurrentThread "
                "failed with error: %d\n", rv);
        return NULL;
    }

    // fprintf (stderr, "dll: return from GetEnv attach \n" );
    return env;
    
}

/**
 * getJNIEnv: A helper function to get the JNIEnv* for the given thread.
 * If no JVM exists, then one will be created. JVM command line arguments
 * are obtained from the LIBHDFS_OPTS environment variable.
 *
 * Implementation note: we rely on POSIX thread-local storage (tls).
 * This allows us to associate a destructor function with each thread, that
 * will detach the thread from the Java VM when the thread terminates.  If we
 * failt to do this, it will cause a memory leak.
 *
 * However, POSIX TLS is not the most efficient way to do things.  It requires a
 * key to be initialized before it can be used.  Since we don't know if this key
 * is initialized at the start of this function, we have to lock a mutex first
 * and check.  Luckily, most operating systems support the more efficient
 * __thread construct, which is initialized by the linker.
 *
 * @param: None.
 * @return The JNIEnv* corresponding to the thread.
 */

JNIEnv* getJNIEnv(void)
{
    JNIEnv *env = NULL;
    HDFSTLS *tls = NULL;
    int ret = 0;
    jint rv = 0;
    
#ifdef WIN32
    DWORD dwWaitResult; 

	maybePerformStaticLibInit();

	tls = TlsGetValue(hdfs_dwTlsIndex1); 
    if (tls) return tls->env;
#endif

#ifdef HAVE_BETTER_TLS
    static __thread HDFSTLS *quickTls = NULL;
    if (quickTls) return quickTls->env;
#endif

#ifndef WIN32

    pthread_once(&hdfs_threadInit_Once, Make_Thread_Key);
    
    if (!hdfs_gTlsKeyInitialized)
        return NULL;
    
    tls = pthread_getspecific(hdfs_gTlsKey);
    if (tls) {
        return tls->env;
    }
    
#endif

    if (!hdfs_InitLib) { 
        LOCK_JVM_MUTEX();
        env = getGlobalJNIEnv();
        UNLOCK_JVM_MUTEX();
    } else {
        rv = (*hdfs_JVM)->AttachCurrentThread(hdfs_JVM, (void**) &env, 0);
        if (rv != 0) {
            fprintf(stderr, "Call to AttachCurrentThread "
                    "failed with error: %d\n", rv);
            return NULL;
        }    
    }
    
    if (!env) {
        fprintf(stderr, "getJNIEnv: getGlobalJNIEnv failed\n");
        return NULL;
    }
        
    tls = calloc ( 1, sizeof(HDFSTLS) );
    if (!tls) {
        fprintf(stderr, "getJNIEnv: OOM allocating %zd bytes\n",
                sizeof(HDFSTLS) );
        return NULL;
    }

    tls->env = env;

#ifdef WIN32
    // fprintf (stderr, "dll: save environment\n" );
    if (!TlsSetValue(hdfs_dwTlsIndex1, tls))
         return NULL;    
    return env;
#endif
            
#ifdef HAVE_BETTER_TLS
    quickTls = tls;
    return env;
#endif

#ifndef WIN32
    ret = pthread_setspecific(hdfs_gTlsKey, tls);
    if (ret) {
        fprintf(stderr, "getJNIEnv: pthread_setspecific failed with "
            "error code %d\n", ret);
        hdfsThreadDestructor(tls);
        return NULL;
    }
#endif

    return env;

}
