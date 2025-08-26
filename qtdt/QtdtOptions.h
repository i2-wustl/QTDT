////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtOptions.h 
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
 
#ifndef __QTDT_OPTIONS__
#define __QTDT_OPTIONS__

#include "IBD.h"
#include "Random.h"

#include "stdio.h"

// Constants for selecting regression models
#define QTDT_NONE              0
#define QTDT_ALLISON           1
#define QTDT_FULKER            2
#define QTDT_ORTHOGONAL        4
#define QTDT_GENERAL           8
#define QTDT_WITHIN           16
#define QTDT_STRATIFICATION   32
#define QTDT_RABINOWITZ       64
#define QTDT_MONKS           128
#define QTDT_TOTAL           256
#define AFF_TDT_ZTEST        512
#define AFF_TDT_LINKAGE     1024

// Convenient groupings and bit-masks for similar tests
#define RABINOWITZ_FAMILY     (QTDT_RABINOWITZ | QTDT_MONKS)
#define ORTHOGONAL_FAMILY     (QTDT_GENERAL | QTDT_ORTHOGONAL | QTDT_WITHIN | \
                               QTDT_STRATIFICATION)
#define FULKER_FAMILY         (ORTHOGONAL_FAMILY | QTDT_FULKER)
#define ZTEST_FAMILY          (AFF_TDT_ZTEST | AFF_TDT_LINKAGE)

// Properties of various tests
#define QTDT_ALLOW_DOMINANCE    (QTDT_ALLISON | QTDT_TOTAL | FULKER_FAMILY)
#define QTDT_ALLOW_IMPRINTING   (QTDT_ALLISON | QTDT_TOTAL | ORTHOGONAL_FAMILY)
#define QTDT_ALLOW_MULTIALLELIC (QTDT_TOTAL | FULKER_FAMILY)
#define QTDT_ALLOW_FOUNDERS     (FULKER_FAMILY | QTDT_MONKS | QTDT_TOTAL | ZTEST_FAMILY)
#define QTDT_ALLOW_PEDIGREE     (ORTHOGONAL_FAMILY | ZTEST_FAMILY | QTDT_MONKS)

// Constants for selecting linkage models
#define LINK_NONE              0
#define LINK_SIMPLE            1

// Constants for selecting variance components
#define VAR_OFF                0
#define VAR_NON_SHARED         1
#define VAR_COMMON_ENV         2
#define VAR_NUCLEAR            4
#define VAR_TWIN               8
#define VAR_POLYGENIC         16
#define VAR_ADDITIVE          32
#define VAR_DOMINANCE         64

// Constants for selecting test flavours
#define FLAVOUR_ASSOC          0
#define FLAVOUR_LINKAGE        1
#define FLAVOUR_VARIANCE       2

// Constants for selecting covariates
#define COVAR_NONE             0
#define COVAR_USER             1
#define COVAR_SEX              2
#define COVAR_PARENTS          4

// Allelic transmission scoring
#define TRANSMISSION_PEDIGREE  0
#define TRANSMISSION_NUCLEAR   1

class QtdtOptions
   {
   public:
      static int imprinting;
      static int associationModel;
      static int linkageModel;
      static int flavour;
      static int covariates;
      static int minProbands;
      static int varComponentNull;
      static int varComponentFull;
      static int mcPermutations;
      static int nullSamples;
      static int numericMinimizer;
      static int transmissionScoring;
      static int randomSeed;
      static bool firstAlleleOnly;
      static bool multiAllelic;
      static bool dominance;
      static bool composite;
      static bool diagnosticQ;
      static bool excludeFounderPhenotypes;
      static bool useDeviates;
      static char * whichStatistic;
      static char * whichTest;
      static FILE * rawDataFile;
      static Random rand;

      static bool traitWise, markerWise, alleleWise;

      static IBDTable ibdRepository;

      QtdtOptions();
      ~QtdtOptions();

      static void ConsistencyCheck();
   };

#endif
 
