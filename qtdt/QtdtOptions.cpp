////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtOptions.cpp 
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
 
#include "QtdtOptions.h"
#include "Error.h"

int            QtdtOptions::imprinting = 0;
int            QtdtOptions::associationModel = QTDT_GENERAL;
int            QtdtOptions::linkageModel = LINK_NONE;
int            QtdtOptions::flavour = FLAVOUR_ASSOC;
int            QtdtOptions::covariates = COVAR_USER;
int            QtdtOptions::minProbands = 30;
int            QtdtOptions::varComponentNull = 0;
int            QtdtOptions::varComponentFull = 0;
int            QtdtOptions::mcPermutations = 0;
int            QtdtOptions::nullSamples = 0;
int            QtdtOptions::numericMinimizer = 0;
int            QtdtOptions::transmissionScoring = 0;
int            QtdtOptions::randomSeed = 123456;
bool           QtdtOptions::excludeFounderPhenotypes = 0;
bool           QtdtOptions::multiAllelic = 0;
bool           QtdtOptions::dominance = 0;
bool           QtdtOptions::firstAlleleOnly = 0;
bool           QtdtOptions::composite = 0;
bool           QtdtOptions::diagnosticQ = 0;
bool           QtdtOptions::useDeviates = 0;
Random         QtdtOptions::rand;
IBDTable       QtdtOptions::ibdRepository;
FILE *         QtdtOptions::rawDataFile = NULL;
char *         QtdtOptions::whichStatistic = "Rsq";
char *         QtdtOptions::whichTest = "F";

bool   QtdtOptions::traitWise = true;
bool   QtdtOptions::markerWise = true;
bool   QtdtOptions::alleleWise = true;

QtdtOptions::QtdtOptions()
   { };

QtdtOptions::~QtdtOptions()
   { };

void QtdtOptions::ConsistencyCheck()
   {
   varComponentFull |= varComponentNull;

   if (imprinting && !(associationModel & QTDT_ALLOW_IMPRINTING))
      error("The requested association model cannot account for imprinting");
   if (dominance &&  !(associationModel & QTDT_ALLOW_DOMINANCE))
      error("The requested association model does not allow dominance");
   if (varComponentNull && associationModel & (RABINOWITZ_FAMILY | ZTEST_FAMILY))
       error("The association model is incompatible with variance components");
   if (multiAllelic && associationModel && associationModel & ~QTDT_ALLOW_MULTIALLELIC)
      error("A multi-allelic version of the requested test is not available");

#ifdef _UNSAFE_
   if (flavour == FLAVOUR_ASSOC && !associationModel)
      error("To perform an association test a model must be selected");
   if (flavour == FLAVOUR_LINKAGE && !linkageModel)
      error("To perform a linkage test a model must be selected");
   if (flavour == FLAVOUR_VARIANCE && varComponentFull == varComponentNull)
      error("No alternative hypothesis for variance components test");
   if (composite && (linkageModel || covariates & COVAR_PARENTS ||
                     varComponentFull || associationModel & RABINOWITZ_FAMILY))
      error("Composite phenotype construction is currently only compatible\n"
            "with simple association analysis");
#else
      flavour = varComponentFull != varComponentNull ?
                FLAVOUR_VARIANCE : FLAVOUR_ASSOC;

      if (flavour == FLAVOUR_ASSOC && !associationModel)
         error("No hypothesis specified");
#endif

   if ((mcPermutations || nullSamples) && (flavour != FLAVOUR_ASSOC ||
        associationModel & (QTDT_ALLISON | QTDT_TOTAL)))
      error("Monte-Carlo permutation not available with current selections");

   if (varComponentFull && varComponentNull == 0)
      error("A model for variances under the null must be provided\n");

   if (varComponentFull)
      {
      whichStatistic = "-LnLk";
      whichTest = "ChiSq";
      }

   if (associationModel & (RABINOWITZ_FAMILY | ZTEST_FAMILY))
      {
      whichStatistic = "Rab";
      whichTest = "T";
      }

   if (associationModel & ~(RABINOWITZ_FAMILY) && useDeviates)
      error("The --deviates option is only available for the Monks and Rabinowitz tests");

   if (useDeviates)
      covariates = COVAR_NONE;

   traitWise  = composite == false;
   markerWise = linkageModel || associationModel ||
                varComponentFull & (VAR_ADDITIVE | VAR_DOMINANCE);
   alleleWise = associationModel != QTDT_NONE && !multiAllelic;
   };
 
