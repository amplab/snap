// ComputeMD5.cpp : Compute the MD5 hash of a file.
//

#include "stdafx.h"

void usage()
{
    fprintf(stderr, "usage: ComputeMD5 inputFilename\n");
}

HANDLE hFile;
int main(int argc, char* argv[])
{
    if (2 != argc) usage();

    HCRYPTPROV hProv = NULL;
    BOOL CACWorked;

    CACWorked = CryptAcquireContext(&hProv, "ComputeMD5Container", NULL, PROV_RSA_FULL, CRYPT_MACHINE_KEYSET);

    if (!CACWorked && NTE_BAD_KEYSET == GetLastError()) {
        CACWorked = CryptAcquireContext(&hProv, "ComputeMD5Container", NULL, PROV_RSA_FULL, CRYPT_NEWKEYSET | CRYPT_MACHINE_KEYSET);
    }

    if (!CACWorked) {
        fprintf(stderr, "Unable to CryptAcquireKeyset, %d\n", GetLastError());
        goto done;
    }

    HANDLE hFile = CreateFile(argv[1], GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        fprintf(stderr, "Unable to open file '%s', %d\n", argv[1], GetLastError());
        goto done;
    }

    BY_HANDLE_FILE_INFORMATION fileInfo;

    if (!GetFileInformationByHandle(hFile, &fileInfo)) {
        fprintf(stderr, "Unable to GetFileInformationByHandle, %d\n", GetLastError());
        goto done;
    }

    ULARGE_INTEGER liFileSize;
    liFileSize.HighPart = fileInfo.nFileSizeHigh;
    liFileSize.LowPart = fileInfo.nFileSizeLow;

    HANDLE hMapping = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (NULL == hMapping) {
        fprintf(stderr, "Unable to CreateFileMapping, %d\n", GetLastError());
        goto done;
    }


    HCRYPTHASH hHash;
    if (!CryptCreateHash(hProv, CALG_MD5, NULL, 0, &hHash)) {
        fprintf(stderr, "CryptCreateHash failed, %d\n", GetLastError());
        goto done;
    }


    //
    // CryptHashData takes a 32 bit input size, so we need to feed the file to it in chunks.
    //
    ULARGE_INTEGER fileOffset;
    fileOffset.QuadPart = 0;
    while (fileOffset.QuadPart < liFileSize.QuadPart) {
        DWORD amountToHash;
        if (liFileSize.QuadPart - fileOffset.QuadPart < 0x70000000) {
            amountToHash = (DWORD)(liFileSize.QuadPart - fileOffset.QuadPart);
        } else {
            amountToHash = 0x70000000;
        }

        const BYTE *fileData = (const BYTE *)MapViewOfFile(hMapping, FILE_MAP_READ, fileOffset.HighPart, fileOffset.LowPart, amountToHash);
        if (NULL == fileData) {
            fprintf(stderr, "Unable to MapViewOfFile, %d\n", GetLastError());
            goto done;
        }

        if (!CryptHashData(hHash, fileData, amountToHash, 0)) {
            fprintf(stderr, "CryptHashData failed, %d\n", GetLastError());
            goto done;
        }

        if (!UnmapViewOfFile(fileData)) {
            fprintf(stderr, "Unmap view of file failed, %d\n", GetLastError());
        }

        fileOffset.QuadPart += amountToHash;
    }


    const DWORD hashValueSizeInBytes = 16;
    BYTE hashValue[hashValueSizeInBytes];
    DWORD hashValueSize = hashValueSizeInBytes;

    if (!CryptGetHashParam(hHash, HP_HASHVAL, hashValue, &hashValueSize, 0)) {
        fprintf(stderr, "CryptGetHash failed, %d\n", GetLastError());
        goto done;
    }

    for (DWORD i = 0; i < hashValueSizeInBytes; i++) {
        printf("%02x", hashValue[i]);
    }


done:
    if (NULL != hProv) {
        CryptReleaseContext(hProv, 0);
    }
}

