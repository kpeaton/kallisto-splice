To compile on a Windows 64-bit machine using Visual Studio 2013:

 - Install HDF5 (64-bit) compiled in VS2013: https://www.hdfgroup.org/HDF5/release/obtain5.html

   - In "VC++ Directories": - Add HDF5 include directory to "Include Directories"
                            - Add HDF5 lib directory to "Library Directories"
   - In "C/C++" -> "General": - Add HDF5 include directory to "Additional Include Directories"
   - In "Linker" -> "General": - Add HDF5 lib directory to "Additional Library Directories"
   - In "Linker" -> "Input": - Add szip.lib, zlib.lib, and hdf5.lib to "Additional Dependencies"

 - Install the Microsft SDK v7.1 for Windows: https://www.microsoft.com/en-us/download/details.aspx?id=8279

   - In "Linker" -> "Input": - Add WS2_32.Lib and Psapi.Lib to "Additional Dependencies"

 - The pthread folder included in src\win64 is from: ftp://sourceware.org/pub/pthreads-win32/prebuilt-dll-2-9-1-release/

   - In "C/C++" -> "General": - Add src\win64\pthread\include to "Additional Include Directories"
   - In "Linker" -> "General": - Add src\win64\pthread\lib\x64 to "Additional Library Directories"
   - In "Linker" -> "Input": - Add pthreadVC2.lib to "Additional Dependencies"

 - The headers getopt.h, strings.h, and unistd.h included in src\win64 are drop-in replacements for Unix headers
   (along with dirent.c/.h, dlfcn.c/.h, and times.c/.h, which have to be compiled as well)

   - In "VC++ Directories": - Add src\win64 to "Include Directories"

 - In "C\C++" -> "Preprocessor": - Add _WINSOCKAPI_ to "Preprocessor Definitions"

 - In "Linker" -> "Command Line": -Add /FORCE:MULTIPLE to "Additional Options"