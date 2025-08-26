////////////////////////////////////////////////////////////////////// 
// libsrc/InputFile.cpp 
// (c) 2000-2008 Goncalo Abecasis
// 
// This file is distributed as part of the QTDT source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile QTDT.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Friday August 29, 2008
// 
 
#include "InputFile.h"

#ifdef __ZLIB_AVAILABLE__

IFILE::IFILE(const char * filename, const char * mode)
   {
   // Some implementations of zlib will not open files that are
   // larger than 2Gb. To ensure support for large (uncompressed)
   // files, we fall-back on the regular fopen when the initial
   // gzopen call fails and the filename does not end in .gz

   gzMode = true;
   gzHandle = gzopen(filename, mode);

   if (gzHandle == NULL)
      {
      int lastchar = 0;

      while (filename[lastchar] != 0) lastchar++;

      if (lastchar >= 3 && filename[lastchar - 3] == '.' &&
                           filename[lastchar - 2] == 'g' &&
                           filename[lastchar - 1] == 'z')
         return;

      gzMode = false;
      handle = fopen(filename, mode);
      }
   };

#endif
 
