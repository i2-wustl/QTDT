////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtUi.h 
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
 
#ifndef _QTDTUI_H_
#define _QTDTUI_H_

#include "Constant.h"
#include "Parameters.h"
#include "Genetics.h"
#include "Pedigree.h"
#include "Qtdt.h"
#include "QtdtParameters.h"

void ShowBanner(bool showReferences);
void Documentation();
void ReadInputFiles(QtdtParameters & p, Pedigree & ped);
void RunTests(QtdtParameters & p, Pedigree & ped);
void OutputSummary(QtdtParameters & p, QtdtResult * best);
void PrintPValue(double p);

#endif
 
