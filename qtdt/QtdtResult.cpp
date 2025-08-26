////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtResult.cpp 
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
 
#include "QtdtResult.h"

void QtdtResult::Clear()
   {
   success = false;
   variance = statisticFull = statisticNull = testValue = 0.0;
   p = 1.0;
   dfNull = dfFull = probandCount = informativeProbands = 0;
   t = m = al = 0;
   ped = NULL;
   }

void QtdtResult::Copy(QtdtResult & r)
   {
   t = r.t;
   m = r.m;
   al = r.al;
   success = r.success;
   variance = r.variance;
   statisticNull = r.statisticNull;
   statisticFull = r.statisticFull;
   testValue = r.testValue;
   p = r.p;
   dfNull = r.dfNull;
   dfFull = r.dfFull;
   probandCount = r.probandCount;
   informativeProbands = r.informativeProbands;
   }

QtdtResult & QtdtResult::operator = (QtdtResult & rhs)
   {
   Copy(rhs);
   return (*this);
   }

int QtdtResult::CompareSignificance(QtdtResult * l, QtdtResult * r)
   {
   if (l->p > r->p)
      return 1;
   if (l->p < r->p)
      return -1;
   return r->probandCount - l->probandCount;
   }

int QtdtResult::CompareAssay(QtdtResult * l, QtdtResult * r)
   {
   if (l->t != r->t)
      return l->t - r->t;
   if (l->m != r->m)
      return l->m - r->m;
   return l->al - r->al;
   }
 
