/*++

Module Name:

    WindowsFileMapper.h

Abstract:

    Header for support code for file mapping on Windows

Authors:

    Bill Bolosky, November, 2012

Environment:

    User mode service.

Revision History:

    
--*/

#pragma once

#ifdef  _MSC_VER
class WindowsFileMapper {
public:
    WindowsFileMapper();

    bool init(const char *fileName);
    const _int64 getFileSize();

    char *createMapping(size_t offset);
    void deleteMapping();
#endif  // _MSC_VER