#ifndef WINCOMPAT_H
#define WINCOMPAT_H
// This header fixes some incompatibilities when compiling for Windows

#if defined(WIN32) || defined(WIN64)

#include <memory>
#include <direct.h>

// Support for use of non-standard type uint
// h5utils.h needs this (and H5Writer.h by extension)
typedef unsigned int uint;

// Support for differences in mkdir function between Linux and Windows
#define mkdir(A, B) _mkdir(A)

// Copied from linux libc sys/stat.h
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)

// Support for getting file stats
typedef struct _stat64 filestat;
#define filestat(A, B) _stati64(A, B)

#else

// Support for getting file stats
typedef struct stat filestat;
#define filestat(A, B) stat(A, B)

#endif

#endif // WINCOMPAT_H
