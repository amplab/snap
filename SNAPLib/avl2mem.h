// The AVLElement class would normally be declared in the avl.cpp file, except that because it's
// a template, it needs to be in the header file.  It can only be accessed (including creation and
// destruction) by the AVLTree friend class.

#include <limits.h>

#pragma warning(disable: 4311)
#pragma warning(disable: 4312)

template<class elementClass, class keyClass> class AVLElement {
public:
				 AVLElement(void);
				~AVLElement(void);

public:
	inline void setParent(AVLElement<elementClass,keyClass>*parent) {parentX = parent;}
	inline void setBalance(AVLBalance balance) {balanceX = balance;}

	inline AVLElement<elementClass,keyClass>*getParent() {return parentX;}
	inline AVLBalance						 getBalance() {return balanceX;}

private:
    AVLBalance	 							 balanceX;
    AVLElement<elementClass,keyClass>		*parentX;

public:
	inline AVLElement<elementClass,keyClass> *getLeft(void) { return left; }
	inline void setLeft(AVLElement<elementClass,keyClass> *val) { left = val; }
	inline AVLElement<elementClass,keyClass> *getRight(void) { return right; }
	inline void setRight(AVLElement<elementClass,keyClass> *val) { right = val; }
	inline elementClass *getElement(void) { return element; }
	inline void setElement(elementClass *val) { element = val; }
	inline keyClass getKey(void) { return getElement()->getKey(); }

private:
    AVLElement<elementClass,keyClass>		*left;
    AVLElement<elementClass,keyClass>		*right;
    elementClass							*element;
};


    template<class elementClass, class keyClass> 
AVLElement<elementClass,keyClass>::AVLElement(void)
{
    setBalance(AVLNew);
    left = right = NULL;
	setParent(NULL);
}

    template<class elementClass, class keyClass> 
AVLElement<elementClass,keyClass>::~AVLElement(void)
{
    _ASSERT(getBalance() == AVLNew);
    _ASSERT(left == NULL && right == NULL && getParent() == NULL);
}

//
//	The standard in-memory allocator:
//
//	Uses a pool if supplied,
//	allocates one if requested,
//	falls back to 'new' if a NULL pool is supplied.
//

template<class nodeClass,class poolClass>
class OptPoolMemAllocator
{
public:
	OptPoolMemAllocator(
		unsigned preallocateSize);
	OptPoolMemAllocator(
		poolClass *avlElementPool);
	~OptPoolMemAllocator();

	nodeClass *Allocate();
	void Free(nodeClass *avlElement);
	void DumpPoolStats(nodeClass *avlElement);
private:
    poolClass		*avlElementPool;
	bool			 poolExternal;
};

template<class nodeClass,class poolClass>
OptPoolMemAllocator<nodeClass,poolClass>::OptPoolMemAllocator(
	unsigned preallocateSize)
{
	poolExternal = FALSE;
	if (ULONG_MAX == preallocateSize) {
		avlElementPool = NULL;
	} else {
	    avlElementPool = new Pool(sizeof(nodeClass));
		CERTIFY_ALLOCATION(avlElementPool);
	    if (preallocateSize && (NULL != avlElementPool)) {
			avlElementPool->preAllocate(preallocateSize);
		}
	}
}

template<class nodeClass,class poolClass>
OptPoolMemAllocator<nodeClass,poolClass>::OptPoolMemAllocator(
	poolClass				*avlElementPool)
{
	poolExternal = TRUE;
	this->avlElementPool = avlElementPool;

    _ASSERT(avlElementPool->getObjectSize() == sizeof(nodeClass));
}

template<class nodeClass,class poolClass>
OptPoolMemAllocator<nodeClass,poolClass>::~OptPoolMemAllocator()
{
	if (NULL != avlElementPool && !poolExternal) {
		delete avlElementPool;
	}
}


template<class nodeClass,class poolClass>
nodeClass *OptPoolMemAllocator<nodeClass,poolClass>::Allocate()
{
	nodeClass *avlElement;
	if (NULL == avlElementPool) {
		avlElement = new nodeClass;
		CERTIFY_ALLOCATION(avlElement);
	} else {
	    avlElement = (nodeClass *)avlElementPool->allocate();
	}
	return avlElement;
}

template<class nodeClass,class poolClass>
void OptPoolMemAllocator<nodeClass,poolClass>::Free(nodeClass *avlElement)
{
	if (avlElementPool) {
		avlElementPool->free((void *)avlElement);
	} else {
		delete avlElement;
	}
}

template<class nodeClass,class poolClass>
void OptPoolMemAllocator<nodeClass,poolClass>::DumpPoolStats(nodeClass *avlElement)
{
	if (NULL == avlElementPool) {
		fprintf(stderr,"Unable to allocate avlElementPool; this AVL tree is essentially useless\n");
	} else {
	    fprintf(stderr,"AVLTree AVLElement pool: %d allocations, %d frees, %d news, objectSize %d\n",
			avlElementPool->numAllocations(),
			avlElementPool->numFrees(),
			avlElementPool->numNews(),
			avlElementPool->getObjectSize());
	}
}

//
//
// Gather definitions together to describe the usual in-memory implementation
//
//

template<class elementClass, class keyClass>
class MemImpl
{
public:
	typedef AVLElement<elementClass,keyClass> *nodeClass;
	typedef Pool poolClass;
	typedef OptPoolMemAllocator<AVLElement<elementClass,keyClass>,poolClass> allocatorClass;
};

