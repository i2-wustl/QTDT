////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtResult.h 
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
 
#ifndef __QTDTRESULT_H__
#define __QTDTRESULT_H__

#include "QtdtOptions.h"
#include "Pedigree.h"

class QtdtResult : protected QtdtOptions
   {
   protected:
      void Copy (QtdtResult & r);
      void Clear();

   public:
      bool     success;
      double   variance;
      double   statisticNull;
      double   statisticFull;
      double   testValue;
      double   p;
      int      dfNull;
      int      dfFull;
      int      probandCount, informativeProbands;
      int      t, m, al;
      Pedigree * ped;

      QtdtResult()
         { Clear(); }
      QtdtResult(QtdtResult & original)
         { Copy(original); }

      QtdtResult & operator = (QtdtResult & rhs);

      bool isInteresting()
         { return p < .1; }

      static int CompareSignificance(QtdtResult *, QtdtResult *);
      static int CompareAssay(QtdtResult *, QtdtResult *);
   };

#endif 
