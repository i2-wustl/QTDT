////////////////////////////////////////////////////////////////////// 
// libsrc/Error.h 
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
 
#ifndef _ERROR_H_
#define _ERROR_H_

// #ifdef __cplusplus
// extern "C" {
// #endif

void error(const char * msg, ...);
void warning(const char * msg, ...);
void numerror(const char * msg, ...);

// #ifdef __cplusplus
//   };
// #endif


#endif
 
