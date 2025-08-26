////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtUi.cpp 
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
 
#include "QtdtUi.h"
#include "QtdtDescribe.h"
#include "Error.h"
#include "Sort.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

int main(int argc, char * argv[])
   {
   QtdtParameters p;

   p.Read(argc, argv);
   ShowBanner(p.showReferences);
   p.Status();
   Documentation();
   p.ConsistencyCheck();

   Pedigree ped;

   ReadInputFiles(p, ped);
   DescribeModels();
   RunTests(p, ped);

   return 0;
   }

void ShowBanner(bool showReferences)
   {
#ifndef VERSION
#define VERSION   "PRIVATE RELEASE"
#endif

   printf("QTDT - Quantitative TDT " VERSION "\n"
          "(c) 1998-2007 Goncalo Abecasis (goncalo@umich.edu)\n\n"
          "This program implements tests described by\n"
          "Abecasis et al, Am J Hum Genet 66:279-292 (2000)\n"
          "Abecasis et al, Eur J Hum Genet 8:545-551 (2000)\n"
          "%s",
          !showReferences ? "and others\n" :
          "Allison, Am J Hum Genet 60:676-690 (1997) [TDTQ5]\n"
          "Fulker et al, Am J Hum Genet 64:259-267 (1999)\n"
          "Monks et al, Am J Hum Genet 66:576-592 (2000)\n"
          "Rabinowitz, Hum Hered 47:342-350 (1997)\n");
   }

void Documentation()
   {
   printf("Online documentation http://www.sph.umich.edu/csg/abecasis/QTDT\n"
          "Comments, bugs: goncalo@umich.edu\n\n");
   }

void ReadInputFiles(QtdtParameters & p, Pedigree & ped)
   {
   ped.Prepare(p.dataFile);
   ped.Load(p.pedFile);
   ped.InheritanceCheck();
   p.ibdRepository.Load(ped, p.ibdFile);
   }

void RunTests(QtdtParameters & p, Pedigree & ped)
   {
   Qtdt engine;

   bool rab = (bool) (p.associationModel & (RABINOWITZ_FAMILY | ZTEST_FAMILY));

   char nullTag[20], fullTag[20], dfFullTag[20], dfNullTag[20], alleleTag[20];

   sprintf(nullTag, rab ? "Sum" : "%s(0)", p.whichStatistic);
   sprintf(fullTag, rab ? "StDev" : "%s(%s)", p.whichStatistic,
           engine.ModelAbbreviation());
   sprintf(dfNullTag, rab ? "df" : "df(0)");
   sprintf(dfFullTag, rab ? "-" : "df(%s)", engine.ModelAbbreviation());

   QtdtOptions::rawDataFile = p.skipRawFile ? NULL : fopen(p.rawFile, "wt");

   QtdtResult * best = NULL, overallBest;

   if (p.bestAssocs) best = new QtdtResult[p.bestAssocs];
   if (p.mcPermutations) engine.SetupPermutations(ped);

   int nullSampleCount = 0;

   int totalPhenotypes = QtdtBroker::CountPhenotypes();

   if (totalPhenotypes == 0)
      error("Pedigree includes no phenotypic data");

   if (ped.markerCount == 0 && p.markerWise)
      error("Pedigree includes no marker data");

   if (p.multiAllelic)           // Rare alleles are grouped to preserve power
      ped.LumpAlleles(0.05);
   else if (p.firstAlleleOnly)   // Sort alleles so we test the most common allele
      ped.LumpAlleles(0.0);
   else
      ped.LumpAlleles(0.0, false);

   long tests = 0;

   do {
      if (p.nullSamples)
         printf("Sampling from the null distribution (#%d)\n"
                "*****************************************\n",
                ++nullSampleCount);

      for (int t = 0; t < totalPhenotypes && (p.traitWise || t == 0); t++)
         {
         printf("Testing trait: %30s\n"
                "=============================================\n\n",
                p.traitWise ? (const char *) QtdtBroker::PhenotypeLabel(t) : "Composite");

         for (int m = 0; m < (p.markerWise ? Person::markerCount : 1); m++)
            {
            if (p.markerWise)
               printf("Testing marker: %29s\n"
                      "---------------------------------------------\n\n",
                      (const char *) Person::markerNames[m]);

            int alleleCount = p.markerWise ?
                              p.firstAlleleOnly ? 1 : ped.CountAlleles(m) : 0;

            printf("%7s %7s %8s %7s %8s %7s %7s\n",
                   "Allele", dfNullTag, nullTag, dfFullTag, fullTag,
                   p.whichTest, "p");

            for (int al = 1; al <= (p.alleleWise ? alleleCount : 1); al++)
               {
               sprintf(alleleTag,
                       p.markerWise ? (p.alleleWise ? "%s" : "All") : "N/A",
                       p.alleleWise ? (const char *) ped.GetMarkerInfo(m)->GetAlleleLabel(al) : "");

               engine.TestAllele(ped, t, m, al);

               if (engine.success)
                  {
                  if (p.composite) printf("%s :\n", engine.multiString);
                  printf("%7s %7d %8.2f %7d %8.2f %7.2f ",
                         alleleTag,
                         engine.dfNull,
                         engine.statisticNull,
                         engine.dfFull,
                         engine.statisticFull,
                         engine.testValue);

                  if (engine.p <= 0.000095)
                     printf("%7.0e", engine.p);
                  else if (engine.p <= 0.10 || p.allProbabilities)
                     printf("%7.4f", engine.p);
                  else
                     printf("%7s", "");

                  if (engine.probandCount == engine.informativeProbands)
                     printf("  (%d probands)\n", engine.probandCount);
                  else
                     printf("  (%4d/%d probands)\n", engine.informativeProbands,
                            engine.probandCount);

                  if (p.bestAssocs &&
                      best[p.bestAssocs-1].p >= engine.p)
                     {
                     best[p.bestAssocs-1] = engine;
                     QuickSort(best, p.bestAssocs, sizeof(QtdtResult),
                        COMPAREFUNC(QtdtResult::CompareSignificance));
                     }

                  if (overallBest.p >= engine.p)
                     overallBest = engine;

                  tests++;
                  }
               else
                  if (engine.informativeProbands)
                     {
                     if (engine.probandCount == engine.informativeProbands)
                        printf("%7s   *** not tested ***    "
                               "                           (%d probands)\n",
                               alleleTag, engine.probandCount);
                     else
                        printf("%7s   *** not tested ***    "
                               "                           (%4d/%d probands)\n",
                               alleleTag, engine.informativeProbands,
                               engine.probandCount);
                     }

               }
            printf("\n");
            }
         printf("\n");
         }
   } while (p.nullSamples != 0 && p.nullSamples > nullSampleCount);

   if (QtdtOptions::rawDataFile != NULL)
      fclose(QtdtOptions::rawDataFile);

   if (p.bestAssocs)
      {
      OutputSummary(p, best);
      delete [] best;
      }

   time_t t;
   time(&t);

   printf("Run completed on %s", ctime(&t));

   if (tests)
      {
      printf("%d tests carried out\n\n", tests);
      printf("The most significant result refers to:\n");
      printf("    Trait: %s\n", (const char *) QtdtBroker::PhenotypeLabel(overallBest.t));
      if (p.markerWise)
         printf("   Marker: %s\n", (const char *) ped.markerNames[overallBest.m]);
      if (p.alleleWise)
         printf("   Allele: %s\n", (const char *) ped.GetMarkerInfo(overallBest.m)->GetAlleleLabel(overallBest.al));
      printf("%9s: %.3f\n", p.whichTest, overallBest.testValue);
      printf("  p-value: "); PrintPValue(overallBest.p); printf("\n\n");

      if (overallBest.p > 1e-20 && overallBest.p < 1.0)
         printf("Overall Bonferroni significance level: "),
         PrintPValue(1.0 - exp(log(1.0 - overallBest.p) * tests)),
         printf("\n");

      if (p.mcPermutations)
         {
         printf("Overall empirical significance level:  ");
         PrintPValue(engine.GlobalSignificance(overallBest.p));
         printf("\n\nEmpirical thresholds:");
         printf("\n     .05 global significance level requires single test p-value of ");
         PrintPValue(engine.GlobalThreshold(0.05));
         printf("\n     .01 global significance level requires single test p-value of ");
         PrintPValue(engine.GlobalThreshold(0.01));
         printf("\n    .001 global significance level requires single test p-value of ");
         PrintPValue(engine.GlobalThreshold(0.001));
         printf("\n\nEmpirical significance based on McIntyre et al, Genet Epidemiol 19:18-29 (2001)\n");
         }

      printf("\n");
      }
   else
      printf("No tests carried out\n\n");
   }

void OutputSummary(QtdtParameters & p, QtdtResult * best)
   {
   if (p.summaryFile[0] == 0) return;

   QuickSort(best, p.bestAssocs, sizeof(QtdtResult),
            COMPAREFUNC(QtdtResult::CompareAssay));

   FILE * summary = fopen (p.summaryFile, "rt");
   int isNew = summary == NULL;
   if (!isNew) fclose(summary);

   summary = fopen(p.summaryFile, "at");
   if (summary == NULL) error("Can't open summary file");

   if (isNew)
      {
      bool  rab = (p.associationModel & (RABINOWITZ_FAMILY | ZTEST_FAMILY));
      char  nullTag[20], fullTag[20], dfFullTag[20], dfNullTag[20];
      sprintf(nullTag, rab ? "Sum" : "%s(0)", p.whichStatistic);
      sprintf(fullTag, rab ? "StDev" : "%s(1)", p.whichStatistic);
      sprintf(dfNullTag, rab ? "df" : "df(0)");
      sprintf(dfFullTag, rab ? "-" : "df(1)");

      fprintf(summary, "%15s %15s %10s %10s %10s "
         "%10s %10s %10s %10s %10s %s\n",
         "Trait", "Marker", "Allele", "Probands",
         dfNullTag, nullTag, dfFullTag, fullTag,
         p.whichTest, "p", "CommandLine");
      }

   printf("The %d most significant associations are:\n", p.bestAssocs);

   printf("%15s %15s %6s %6s %6s %6s\n",
         "Trait", "Marker", "Allele", "N", p.whichTest, "p");

   for (int b=0; b<p.bestAssocs; b++)
      {
      int precisionNull = 0;
      int precisionFull = 0;

      if (best[b].statisticNull > 1)
         precisionNull = (int) log10(best[b].statisticNull);
      precisionNull = precisionNull > 6 ? 0 : 6 - precisionNull;

      if (best[b].statisticFull > 1)
         precisionFull = (int) log10(best[b].statisticFull);
      precisionFull = precisionFull > 6 ? 0 : 6 - precisionFull;

      if (best[b].al)
      fprintf(summary, "%15s %15s %10d %10d %10d "
         "%10.*f %10d %10.*f %10.6f %10.4e %s\n",
         (const char *) QtdtBroker::PhenotypeLabel(best[b].t),
         (const char *) Person::markerNames[best[b].m],
         (const char *) Pedigree::GetMarkerInfo(best[b].m)->GetAlleleLabel(best[b].al),
         best[b].probandCount,
         best[b].dfNull, precisionNull, best[b].statisticNull,
         best[b].dfFull, precisionFull, best[b].statisticFull,
         best[b].testValue, best[b].p, p.string);

      if (best[b].p <= 0.10)
         {
         printf("%15s %15s %6s %6d %6.2f %6.4f\n",
            (const char *) QtdtBroker::PhenotypeLabel(best[b].t),
            (const char *) Person::markerNames[best[b].m],
            (const char *) Pedigree::GetMarkerInfo(best[b].m)->GetAlleleLabel(best[b].al),
            best[b].probandCount,
            best[b].testValue, best[b].p);
         }
      }

   fclose(summary);
   printf("\n");
   }

void PrintPValue(double p)
   {
   printf(p <= 0.000095 ? "%.0e" : "%.4f", p);
   }

 
