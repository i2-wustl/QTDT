////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtBroker.h 
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
 
#ifndef __QTDTBROKER_H__

// The QtdtBroker class merges binary and quantitative traits seamlessly
//

#define  UNDERLYING_QUANTITATIVE    0
#define  UNDERLYING_BINARY          1

#include "Pedigree.h"
#include "QtdtOptions.h"

class QtdtBroker : public QtdtOptions
   {
   public:
      static int    CountPhenotypes();

      static String PhenotypeLabel(int t);

      static bool   isFullyPhenotyped(Person & p);
      static bool   isPhenotyped(Person & p, int t);
      static double GetPhenotype(Person & p, int t);

      static bool   isBinaryPhenotype(int t);
   };

#endif 
