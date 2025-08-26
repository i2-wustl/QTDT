////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtParameters.h 
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
 
#ifndef __QTDT_PARAMETERS__
#define __QTDT_PARAMETERS__

#include "Parameters.h"
#include "QtdtOptions.h"

class QtdtParameters : public ParameterList, public QtdtOptions
   {
   public:
      static OptionList linkOptions[], associationOptions[],
                        imprintingOptions[], varianceOptions[],
                        flavourOptions[], covariateOptions[],
                        minimizerOptions[], transmissionOptions[];
      static LongParameterList longOptions[];

      String dataFile;
      String pedFile;
      String ibdFile;
      String rawFile;
      String summaryFile;

      int bestAssocs;

      static bool showReferences;
      static bool allProbabilities;
      static bool skipRawFile;

      QtdtParameters();

      virtual void Read(int argc, char ** argv, int start = 1);
   };

#endif
 
