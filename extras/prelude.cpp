////////////////////////////////////////////////////////////////////// 
// extras/prelude.cpp 
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
 
#include "Constant.h"
#include "Parameters.h"
#include "Pedigree.h"
#include "MathConstant.h"
#include "Error.h"
#include "MapFunction.h"

#include <stdio.h>

class AlleleFrequencyParameter : public Parameter
   {
   public:
   AlleleFrequencyParameter(char c, const char * desc, int & how, String & file);

   virtual void Status();

   protected:
      String * filename;
      String status;

      virtual void Translate(const char * value);
   };

void PrepareRecombinationFractions(Pedigree & ped, double defaultTheta);

void WriteSimwalkDataFile(Pedigree & ped, int freqs);
void WriteSimwalkPedigreeFile(Pedigree & ped);
void WriteSimwalkBatchFile();

void WriteLinkageDataFile(Pedigree & ped, int freqs);
void WriteLinkagePedigreeFile(Pedigree & ped);
void WriteGenehunterBatchFile();

int main(int argc, char * argv[])
   {
   printf("Prelude - (c) 1999-2002 Goncalo Abecasis\n"
          "Preprocessor for external IBD sources\n");

   String pedName = QTDTPED;
   String dataName = QTDTDATA;
   String mapName = "qtdt.map";
   String freqName;
   double theta = 0.01;
   int    freqs = 0;

   ParameterList pl;

   pl.Add(new StringParameter('p', "Pedigree File", pedName));
   pl.Add(new StringParameter('d', "Data File", dataName));
   pl.Add(new StringParameter('m', "Map File (optional)", mapName));
   pl.Add(new DoubleParameter('t', "Theta between markers", theta));
   pl.Add(new AlleleFrequencyParameter('a', "Allele Frequencies", freqs, freqName));

   pl.Read(argc, argv);
   pl.Status();

   Pedigree ped;

   // Load data file
   ped.Prepare(dataName);

   // Load pedigree file
   ped.Load(pedName);
   ped.LumpAlleles(0.0, false);

   // Load map file
   ped.LoadMarkerMap(mapName);

   // Check for inheritance problems
   ped.InheritanceCheck();

   // Combine each MZ twin set into a single individual
   ped.MergeTwins();

   // Load allele frequencies (if file is available) and estimate them if required
   printf("\nIf an allele frequency file is provided, or a linkage format data\n"
          "file is used, prelude will use the frequencies specified there.\n\n");
   ped.LoadAlleleFrequencies(freqName);
   ped.EstimateFrequencies(freqs);

   // Calculate inter-marker recombination fractions
   PrepareRecombinationFractions(ped, theta);

   WriteSimwalkBatchFile();
   WriteSimwalkDataFile(ped, freqs);
   WriteSimwalkPedigreeFile(ped);

   printf("\n");
   WriteGenehunterBatchFile();
   WriteLinkageDataFile(ped, freqs);
   WriteLinkagePedigreeFile(ped);

   printf("\n\nYou may now run your IBD analysis.\n\n"
          "   To run SimWalk2 issue the command:    simwalk2\n"
          "   To run Genehunter2 issue the command: gh2 < genehunter.in\n\n"
          "When you are finished, run finale\n\n");

   printf("REMEMBER:\n"
          "   If your data consists of small and medium sized pedigrees,\n"
          "   you can use Merlin to estimate IBD in a single step.\n\n"
          "   Use the following command:\n\n"
          "   merlin -d %s -p %s -m %s --markerNames --ibd\n\n",
          (const char *) dataName, (const char *) pedName,
          (const char *) mapName);

   return 0;
   };

Vector theta;

void PrepareRecombinationFractions(Pedigree & ped, double defaultTheta)
   {
   IntArray markerOrder;
   int chromosome = -2;

   printf("If a map file is provided, or a linkage format data file is used,\n"
          "prelude will use the recombination fractions specified there.\n\n");

   if (ped.MarkerPositionsAvailable())
      ped.SortMarkersInMapOrder(markerOrder, chromosome);

   if (markerOrder.Length() == 0)
      {
      printf("Assuming markers are separated by %.3f recombination fraction\n\n",
             defaultTheta);
      theta.Dimension(ped.markerCount - 1);
      theta.Set(defaultTheta);
      return;
      }

   if (markerOrder.Length() != ped.markerCount)
      error("Not all markers could be placed");

   if (!markerOrder.isAscending())
      error("Marker orders in pedigree and genetic map do not match");

   double lastPosition = ped.GetMarkerInfo(markerOrder[0])->position;
   for (int i = 1; i < markerOrder.Length(); i++)
      {
      double position = ped.GetMarkerInfo(markerOrder[i])->position;

      theta.Push(DistanceToRecombination(position - lastPosition));

      lastPosition = position;
      }
   }

void WriteSimwalkBatchFile()
   {
   printf("Preparing SimWalk2 Run Specification ...\n");

   FILE * f = fopen("BATCH2.DAT", "wt");

   // Specify that we want to perform IBD analysises
   fprintf(f, "\n000001\n4\n");

   // Specify a descriptive run title
   fprintf(f, "\n000003\nCalculation of IBD probabilities for QTDT\n");

   // Recombination frequencies between markers -- assumed to be equal!
   fprintf(f, "\n000015    "
              "*** Edit these lines for unequal marker spacing ***\n"
              "%d\n", Person::markerCount);
   for (int i = 0; i < Person::markerCount - 1; i++)
      fprintf(f, "%8.6f\n", theta[i]);

   // label for 'affected' individuals
   fprintf(f, "\n000016\nCALC_IBD\n");

   // number of pedigrees to sample
   fprintf(f, "\n000020    "
              "*** Increase the following number for greater accuracy ***\n"
              "1000\n");

   // do not calculate IBD between markers
   fprintf(f, "\n000049\n1\n");

   // end of batch file
   fprintf(f, "\n000050\n\n");

   fclose(f);
   };

void WriteSimwalkDataFile(Pedigree & ped, int freqs)
   {
   printf("Preparing SimWalk2 Data Description ...\n");

   FILE * f = fopen("LOCUS.DAT", "wt");

   // Dummy trait locus where everyone is affected
   fprintf(f, "CALC_IBDAUTOSOME 2 2\n"
              "     ONE 0.99999\n"
              "     TWO 0.00001\n"
              "CALC_IBD1\n"
              "ONE/ONE\n"
              "IGNORED 2\n"
              "ONE/TWO\n"
              "TWO/TWO\n");

   // The actual markers, in the same order as in the data file
   for (int m = 0; m < Person::markerCount; m++)
      {
      MarkerInfo * info = ped.GetMarkerInfo(m);

      int alleleCount = info->freq.Length() - 1;

      // Simwalk2 requires at least two alleles per marker
      if (alleleCount <= 1)
         {
         fprintf(f, "%-8.8sAUTOSOME 2 0\n"
                    "       10.999999\n"
                    "       20.000001\n",
                    (const char *) info->name);
         continue;
         }

      // The actual allele frequencies
      fprintf(f, "%-8.8sAUTOSOME%2d 0\n",
               (const char *) Person::markerNames[m], max(alleleCount, 2));

      // If one allele has frequency one and others have frequency zero,
      // we perform minor adjustment to avoid Simwalk2 complaints
      double min_freq = info->freq.CountIfGreater(0.0) == 1 ? 0.000001 : 0.0;

      for (int i = 1; i <= alleleCount; i++)
         fprintf(f, "%8d%8.6f\n", i, max(min(info->freq[i], 0.999999), min_freq));
      }

   fclose(f);
   };

void WriteSimwalkPedigreeFile(Pedigree & ped)
   {
   printf("Preparing SimWalk2 Pedigree Data ...\n");

   FILE * f = fopen("PEDIGREE.DAT", "wt");

   // Specify the output format
   fprintf(f, "(I5,1X,A8)\n(3A8,1X,2A1,A8,(T36,3A8,:))\n");

   // The major hazards are exact column widths and
   // converting missing data to empty space
   for (int i = 0; i < ped.familyCount; i++)
      {
      IntArray mzTwins;
      char *   key = " MF";

      fprintf(f, "%5d %8.8s\n",
         ped.families[i]->count, (const char *) ped.families[i]->famid);

      for (int j = ped.families[i]->first; j <= ped.families[i]->last; j++)
         {
         if (ped[j].pid.Length() > 8)
            printf("WARNING - Simwalk2 only supports 8 character names, but\n"
                   "          [%s] in family [%s] has an %d character name\n",
                   (const char *) ped[j].famid,
                   (const char *) ped[j].pid,
                   ped[j].pid.Length());

         if (ped[j].isFounder())
            fprintf(f, "%-8.8s%-8.8s%-8.8s %c CALC_IBD",
                    (const char *) ped[j].pid, "", "", key[ped[j].sex]);
         else
            {
            char twin_status = ' ';

            if (ped[j].isMzTwin(ped[j]))
               {
               int k;

               for (k = 0; k < mzTwins.Length(); k++)
                  if (ped[j].isMzTwin(ped[mzTwins[k]]))
                     break;

               if (k == mzTwins.Length())
                  mzTwins.Push(j);

               twin_status = (char) ('A' + k);
               }

            fprintf(f, "%-8.8s%-8.8s%-8.8s %c%cCALC_IBD",
                    (const char *) ped[j].pid,
                    (const char *) ped[j].fatid, (const char *)  ped[j].motid,
                    key[ped[j].sex], twin_status);
            }

         for (int m = 0; m < Person::markerCount; m++)
            {
            if ((m != 0) && (m % 3 == 0))
               fprintf(f,"\n%35s", "");
            if (ped[j].markers[m].isKnown())
              fprintf(f," %3d/%3d", ped[j].markers[m][0], ped[j].markers[m][1]);
            else
              fprintf(f,"        ");
            }
         fprintf(f, "\n");
         }
      }

   fclose(f);
   }

void WriteLinkageDataFile(Pedigree & ped, int freqs)
   {
   printf("Preparing Linkage Data Description ...\n");

   FILE * f = fopen("linkage.dat", "wt");

   // Header lines with number of markers
   // and other information we don't need
   fprintf(f, "%d 0 0 0\n0 0 0 0\n", Person::markerCount + 2);

   // Marker order
   for (int i = 0; i <= Person::markerCount + 1; i++)
      fprintf(f, "%d ", i + 1);

   // Dummy affection locus where everyone is affected
   fprintf(f, "\n1  2  # dummy\n0.0001 0.9999\n1\n0.0 0.0 1.0\n");

   // The actual markers, in the same order as in the data file
   for (int m = 0; m < Person::markerCount; m++)
      {
      MarkerInfo * info = ped.GetMarkerInfo(m);

      int alleleCount = info->freq.Length() - 1;

      fprintf(f, "3 %2d  # %s\n", alleleCount,
                 (const char *) Person::markerNames[m]);

      for (int i = 1; i <= alleleCount; i++)
         fprintf(f, "%8.6f ", info->freq[i]);

      fprintf(f, "\n");
      }

   // Dummy marker locus where everyone is typed
   fprintf(f, "3  2  # undrop  << dummy locus to stop pedigree dropping\n"
              "0.5 0.5\n");

   // Footer -- no sex difference, no interference
   // Default recombination fractions for dummy locus and undrop marker are 0.5
   fprintf(f, "0 0\n0.50 ");

   for (int i = 0; i < Person::markerCount - 1; i++)
      fprintf(f, "%f ", theta[i]);
   fprintf(f, " 0.50\n");

   fclose(f);
   };

void WriteLinkagePedigreeFile(Pedigree & ped)
   {
   printf("Preparing Linkage Pedigree Data ...\n");

   FILE * f = fopen("linkage.ped", "wt");

   for (int i = 0; i < ped.familyCount; i++)
      {
      for (int j = ped.families[i]->first; j <= ped.families[i]->last; j++)
         {
         if (ped[j].isFounder())
            fprintf(f, "%4s %4s    0    0 %d  1 ",
                    (const char *) ped.families[i]->famid,
                    (const char *) ped[j].pid, ped[j].sex);
         else
            fprintf(f, "%4s %4s %4s %4s %d  1 ",
                    (const char *) ped.families[i]->famid,
                    (const char *) ped[j].pid,
                    (const char *) ped[j].father->pid,
                    (const char *) ped[j].mother->pid,
                    ped[j].sex);

         for (int m = 0; m < Person::markerCount; m++)
            if (ped[j].markers[m].isKnown())
              fprintf(f," %2d %2d", ped[j].markers[m][0], ped[j].markers[m][1]);
            else
              fprintf(f,"  0  0");
         fprintf(f, "  2 2\n");
         }
      }

   fclose(f);
   }

void WriteGenehunterBatchFile()
   {
   printf("Preparing Genehunter Run Specification ...\n");

   FILE * f = fopen("genehunter.in", "wt");

   // specify that we need a log file for error-checking
   fprintf(f, "photo genehunter.out\n");

   // increase max bits
   fprintf(f, "max bits 20\n");

   // one increment between markers
   fprintf(f, "increment step 1\n");

   // don't output linkage statistics
   fprintf(f, "display scores off\n");

   // read data file
   fprintf(f, "load markers linkage.dat\n");

   // setup distances
   fprintf(f, "system # To change marker order and/or spacing\n"
              "system # edit the next line, but do not remove\n"
              "system # the USE command, as it is required by\n"
              "system # by finale.\nuse ");

   // Recombination frequencies between markers -- assumed to be equal!
   for (int i = 0; i < Person::markerCount - 1; i++)
      fprintf(f, "%s %8.6f ",
              (const char *) Person::markerNames[i], theta[i]);
   fprintf(f,  "%s 0.5 undrop\n",
              (const char *) Person::markerNames[Person::markerCount-1]);

   // analyse pedigrees
   fprintf(f, "scan pedigree linkage.ped\n");

   // dump ibd information
   fprintf(f, "system rm -f dump.ibd\ndump ibd\ndump.ibd\n");

   fprintf(f, "quit\n");


   fclose(f);
   };

AlleleFrequencyParameter::AlleleFrequencyParameter
   (char c, const char * desc, int & how, String & file) :
   Parameter(c, desc, &how), status("ALL INDIVIDUALS")
   {
   how = FREQ_ALL;
   filename = &file;
   }

void AlleleFrequencyParameter::Translate(const char * str)
   {
   String value(str);
   int * estimator = (int *) var;

   if (value == "a" || value == "A")
      {
      * estimator = FREQ_ALL;
      * filename = "";
      status = "ALL INDIVIDUALS";
      }
   else if (value == "e" || value == "E")
      {
      * estimator = FREQ_EQUAL;
      * filename = "";
      status = "ASSUMED EQUAL";
      }
   else if (value == "f" || value == "F")
      {
      * estimator = FREQ_FOUNDERS;
      * filename = "";
      status = "FOUNDERS ONLY";
      }
   else
      {
      * estimator = FREQ_ALL;
      * filename = value;
      status = value;
      }
   }

void AlleleFrequencyParameter::Status()
   {
   printf("%*s : %*s (-%c[a|e|f|file])\n", nameCol, description,
          statusCol, (const char *) status, ch);
   }

 
