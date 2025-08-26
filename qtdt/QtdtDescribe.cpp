////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtDescribe.cpp 
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
 
#include "QtdtDescribe.h"
#include "QtdtOptions.h"
#include "QtdtLinear.h"
#include "QtdtVar.h"
#include "MathMatrix.h"

void DescribeModels()
   {
   QtdtOptions qtdt;
   QtdtVariances variances;

   printf("The following models will be evaluated...\n");

   if (qtdt.associationModel & ZTEST_FAMILY)
      {
      printf("  NULL MODEL\n"
             "    Sum of W == 0 [among AFFECTEDS]\n\n"
             "  ALTERNATIVE HYPOTHESIS\n"
             "    Sum of W <> 0 [among AFFECTEDS]\n\n"
             "    W is a measure of allelic transmission (Abecasis et al, EJHG, 2000)\n"
             "    Variance estimate for W assumes no linkage%s\n\n"
             "  Test statistic T_affected in Abecasis et al. (EJHG, 2000)\n\n",
             qtdt.associationModel == AFF_TDT_ZTEST ? " disequilibrium" : "");
      return;
      }

   bool rab = (bool) (qtdt.associationModel & RABINOWITZ_FAMILY);

   if (rab)
      printf("  Expected Genotype  = Conditional on %s Genotypes\n"
             "  Expected Phenotype =",
         qtdt.associationModel == QTDT_RABINOWITZ ? "Parental" :
         qtdt.transmissionScoring == TRANSMISSION_NUCLEAR ? "Parental and Sibling" :
         "Ancestral and Sibling");
   else
      printf("  NULL MODEL\n     Means =");

   if (qtdt.useDeviates)
      printf(" 0\n\n");
   else
      {
      Matrix             nullMatrix("nullMatrix");
      QtdtLinearCounter  nullCounter;
      QtdtLinearLabeller nullLabeller(&nullMatrix);

      nullCounter.LinearModel(false);
      nullMatrix.Dimension(1, nullCounter.count);
      nullLabeller.LinearModel(false);

      for (int i = 0; i < nullCounter.count; i++)
         printf("%s %s", i ? " +" : "", nullMatrix.GetColumnLabel(i));
      printf("\n");

      int flag = 0;
      for (int i = 1; i <= qtdt.varComponentNull; i <<=1 )
         if (qtdt.varComponentNull & i)
            printf("%s %s", flag++ ? " +" : " Variances =", variances.GetLabel(i));
      printf(flag ? "\n\n" : "\n");
      }

   if (rab)
      {
      printf("  NULL HYPOTHESIS\n   (Observed - Expected Genotype) *"
             " (Observed - Expected Phenotype) = 0\n\n"
             "  ALTERNATIVE HYPOTHESIS\n   (Observed - Expected Genotype) *"
             " (Observed - Expected Phenotype) <> 0\n\n");

      if (qtdt.associationModel == QTDT_MONKS)
         printf("  Test statistic T_all in Abecasis et al. (EJHG, 2000)\n"
                "  equivalent to T_QPS in Monks and Kaplan (AJHG, 2000)\n\n");
      else
         printf("  Test statistic from Rabinowitz (1997)\n");

      return;
      }

   printf("  FULL MODEL\n");

   printf("     Means =");

   Matrix             fullMatrix("fullMatrix");
   QtdtLinearCounter  fullCounter;
   QtdtLinearLabeller fullLabeller(&fullMatrix);

   fullCounter.LinearModel(true);
   fullMatrix.Dimension(1, fullCounter.count);
   fullLabeller.LinearModel(true);

   for (int i = 0; i < fullCounter.count; i++)
      printf("%s %s", i ? " +" : "", fullMatrix.GetColumnLabel(i));
   printf("\n");

   int flag = 0;
   for (int i = 1; i <= qtdt.varComponentFull; i <<=1 )
      if (qtdt.varComponentFull & i)
        printf("%s %s", flag++ ? " +" : " Variances =", variances.GetLabel(i));
   printf(flag ? "\n\n" : "\n");

   if (qtdt.associationModel == QTDT_GENERAL &&
       qtdt.varComponentFull == (VAR_ADDITIVE | VAR_POLYGENIC | VAR_NON_SHARED))
     printf("  Likelihood ratio statistic from Abecasis et al (AJHG, 2000)\n\n");

   if (qtdt.flavour == FLAVOUR_LINKAGE && qtdt.linkageModel == LINK_SIMPLE)
      printf("** Only positive estimates for Pi*T considered significant\n\n");
   }

 
