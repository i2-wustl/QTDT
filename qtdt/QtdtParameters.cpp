////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtParameters.cpp 
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
 
#include "QtdtParameters.h"
#include "Genetics.h"
#include "Constant.h"
#include "MathNormal.h"
#include "Random.h"

BEGIN_OPTION_LIST(QtdtParameters::varianceOptions)
   {'e', "NON SHARED", VAR_NON_SHARED },
   {'c', "COMMON ENVIRONMENT", VAR_COMMON_ENV},
   {'g', "POLYGENIC", VAR_POLYGENIC },
   {'n', "NUCLEAR UNIT", VAR_NUCLEAR },
   {'t', "TWIN ENVIRONMENT", VAR_TWIN },
   {'a', "ADDITIVE", VAR_ADDITIVE },
   {'d', "DOMINANCE", VAR_DOMINANCE},
   {'-', "NOT MODELLED", VAR_OFF }
END_OPTION_LIST("NOT MODELLED");

BEGIN_OPTION_LIST(QtdtParameters::covariateOptions)
   {'p', "PARENTAL TRAITS", COVAR_PARENTS},
   {'s', "SEX", COVAR_SEX },
   {'u', "USER SPECIFIED", COVAR_USER },
   {'-', "NONE", COVAR_NONE }
END_OPTION_LIST("NONE");

BEGIN_OPTION_LIST(QtdtParameters::minimizerOptions)
   {'f', "FLETCHER-REEVES", NORMAL_FLETCHER_MIN},
   {'n', "NELDER AND MEAD", NORMAL_AMOEBA_MIN},
   {'p', "POWELL", NORMAL_POWELL_MIN}
END_OPTION_LIST("NONE");

BEGIN_OPTION_LIST(QtdtParameters::associationOptions)
   {'a', "ALLISON", QTDT_ALLISON},
   {'d', "DISCRETE TRAIT", AFF_TDT_ZTEST},
   {'f', "FULKER", QTDT_FULKER},
#ifdef _UNSAFE_
   {'l', "LINKAGE TDT", AFF_TDT_LINKAGE},
#endif
   {'m', "MONKS", QTDT_MONKS},
   {'o', "ORTHOGONAL", QTDT_GENERAL},
   {'p', "STRATIFICATION", QTDT_STRATIFICATION},
   {'r', "RABINOWITZ", QTDT_RABINOWITZ},
   {'t', "TOTAL", QTDT_TOTAL},
   {'w', "WITHIN ONLY", QTDT_WITHIN},
   {'-', "NONE", QTDT_NONE}
END_OPTION_LIST("NONE");

BEGIN_OPTION_LIST(QtdtParameters::imprintingOptions)
   {'f', "FULLY MODELLED", I_FULL},
   {'t', "TEST SIGNIFICANCE", I_IMPRINTING},
   {'m', "MATERNAL TEST", I_MATERNAL},
   {'p', "PATERNAL TEST", I_PATERNAL},
   {'-', "NONE", I_NONE}
END_OPTION_LIST("NONE");

BEGIN_OPTION_LIST(QtdtParameters::flavourOptions)
   {'a', "ASSOCIATION", FLAVOUR_ASSOC},
   {'l', "LINKAGE", FLAVOUR_LINKAGE},
   {'v', "VARIANCES", FLAVOUR_VARIANCE}
END_OPTION_LIST("ASSOCIATION")

BEGIN_OPTION_LIST(QtdtParameters::linkOptions)
   {'s', "SIMPLE", LINK_SIMPLE},
   {'-', "NONE", LINK_NONE}
END_OPTION_LIST("NONE")

BEGIN_OPTION_LIST(QtdtParameters::transmissionOptions)
   {'n', "NUCLEAR FAMILY", TRANSMISSION_NUCLEAR},
   {'p', "FULL PEDIGREE", TRANSMISSION_PEDIGREE}
END_OPTION_LIST("NONE")

BEGIN_LONG_PARAMETERS(QtdtParameters::longOptions)
   LONG_PARAMETER("dominance", &dominance)
   LONG_PARAMETER("snp", &firstAlleleOnly)
   LONG_PARAMETER("multi-allelic", &multiAllelic)
   LONG_PARAMETER("deviates", &QtdtParameters::useDeviates)
   LONG_PARAMETER("references", &QtdtParameters::showReferences)
   LONG_PARAMETER("exclude-founder-phenotypes", &QtdtParameters::excludeFounderPhenotypes)
   LONG_PARAMETER("p-values", &QtdtParameters::allProbabilities)
   LONG_PARAMETER("no-regress-tbl", &QtdtParameters::skipRawFile)
END_LONG_PARAMETERS()

QtdtParameters::QtdtParameters() :
   dataFile(QTDTDATA), pedFile(QTDTPED), ibdFile(QTDTIBD),
   rawFile(QTDTRAW), summaryFile("summary.tbl")
   {
   nullSamples = bestAssocs = 0;
   firstAlleleOnly = false;

   Parameter::SetNameLen(25);
   Parameter::SetStatusLen(20);

   Add (new StringParameter('d', "QTDT Data File", dataFile));
   Add (new StringParameter('p', "QTDT Pedigree File", pedFile));
   Add (new StringParameter('i', "QTDT IBD Status File", ibdFile));
   Add (new StringParameter('x', "Missing Value Code", Pedigree::missing));
#ifdef _UNSAFE_
   Add (new StringParameter('r', "Regression Tables File", rawFile));
   Add (new StringParameter('u', "Summary File", summaryFile));
   Add (new IntParameter('b', "Best Results in Summary", bestAssocs));
   Add (new ListParameter('f', "Flavour of Test", flavour, flavourOptions));
   Add (new ListParameter('l', "Linkage Model", linkageModel, linkOptions));
#endif
   Add (new SetParameter('c', "Covariates", covariates, covariateOptions));
   Add (new ListParameter('a', "Association Model", associationModel, associationOptions));
   Add (new SetParameter('v', "Full Model Variances", varComponentFull, varianceOptions));
   Add (new SetParameter('w', "Null Model Variances", varComponentNull, varianceOptions));
   Add (new HiddenSwitch('g', "Genetic Dominance", dominance));
   Add (new ListParameter('o', "Parent of Origin Effects", imprinting, imprintingOptions));
   Add (new IntParameter('m', "Monte-Carlo Permutations", mcPermutations));
   Add (new IntParameter('r', "Random Seed", randomSeed));
   Add (new ListParameter('n', "Numeric Minimizer", numericMinimizer, minimizerOptions));
   Add (new ListParameter('t', "Transmission Scoring", transmissionScoring, transmissionOptions));
   Add (new HiddenSwitch('1', "First Allele Only", firstAlleleOnly));
   Add (new HiddenSwitch('9', "Multi-Allelic Model", multiAllelic));
   Add (new LongParameters("Additional Options", longOptions));
#ifdef _UNSAFE_
   Add (new SwitchParameter('t', "Composite Trait", composite));
   Add (new SwitchParameter('q', "Q diagnostics", diagnosticQ));
   Add (new IntParameter('n', "Null Distribution Samples", nullSamples));
#endif
   }

void QtdtParameters::Read(int argc, char ** argv, int start)
   {
   ParameterList::Read(argc, argv, start);

   if ((associationModel & ~QTDT_ALLOW_PEDIGREE) || imprinting)
      transmissionScoring = TRANSMISSION_NUCLEAR;

   if (associationModel == QTDT_TOTAL)
      transmissionScoring = 2;

   globalRandom.Reset(randomSeed);
   }

bool QtdtParameters::showReferences = 0;
bool QtdtParameters::allProbabilities = false;
bool QtdtParameters::skipRawFile = false;

 
