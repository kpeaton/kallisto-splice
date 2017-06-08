#ifndef WINCOMPAT_H
#define WINCOMPAT_H
// This header should be included when compiling for Windows

#include <memory>
#include <direct.h>

// Support for use of non-standard type uint
// h5utils.h needs this (and H5Writer.h by extension)
typedef unsigned int uint;

// Support for differences in mkdir function between Linux and Windows
#define mkdir(A, B) _mkdir(A)

#if defined(WIN32) || defined(WIN64)
// Copied from linux libc sys/stat.h:
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#endif

#endif // WINCOMPAT_H
