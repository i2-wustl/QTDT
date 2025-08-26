////////////////////////////////////////////////////////////////////// 
// extras/finale.cpp 
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
 
#include    "StringArray.h"
#include    "Constant.h"
#include    "Error.h"

#include    <stdio.h>
#include    <string.h>

int main(int argc, char * argv[])
   {
   printf("Finale - (c) 1999-2002 Goncalo Abecasis\n"
          "Postprocessor for external IBD sources\n\n"
          "Typical usage:\n"
          "==============\n\n"
          "SIMWALK2 output    --- finale IBD-01.*\n"
          "Genehunter2 output --- finale genehunter.in\n\n");

   if (argc <= 1)
      error("Nothing to do - no input files specified\n\n");

   printf("Preparing qtdt.ibd file ...\n\n");
   FILE * output = fopen("qtdt.ibd", "wt");

   bool negative_probabilities = false, twins_mode = false;

   if (argc != 2 || strcmp(argv[1], "genehunter.in"))
      {
      for (int i = 1; i < argc; i++)
         {
         FILE * f = fopen(argv[i], "rt");
         if (f == NULL) error("Opening File %s", argv[i]);

         char buffer[1024] = "";

         // Find the family name ...
         char famid[100] = "famid";
         do if (fscanf(f, " %s ", buffer) == EOF)
            error("Parsing SimWalk2 output file %s\n",
                  (const char *) argv[i]);
         while (strcmp(buffer, "named:") != 0);
         fscanf(f, " %s ", famid);

         // Skip the rest of the blurb
         do if (fscanf(f, " %s ", buffer) == EOF)
            error("Parsing SimWalk2 output file %s\n",
                  (const char *) argv[i]);
         while (buffer[0] != '_');

         // Find out if this file was created by a new version of simwalk2
         fscanf(f, " %s ", buffer);
         bool v286 = strcmp(buffer, "Pair-of-Individuals") == 0;

         bool done = false;
         while (!done)
            {
            do if (fscanf(f, " %s ", buffer) == EOF)
               error("Parsing SimWalk2 output file %s\n",
                     (const char *) argv[i]);
            while (strcmp(buffer, "COEF.") != 0);

            // Retrieve information on pair of relatives for the first marker
            char pid1[100] = "pid1", pid2[100] = "pid2";
            char marker[100] = "marker";
            char ibd0[100] = "ibd0", ibd1[100] = "ibd1", ibd2[100] = "ibd2";

            if (v286)
               fscanf(f, " %s %s %*s",
                      pid1, pid2);
            else
               fscanf(f, " %s %s %s %*s %s %s %s %*s",
                      pid1, pid2, marker, ibd0, ibd1, ibd2);

            // Find out if we have MZ twins twin1$twin2
            StringArray twins1, twins2;

            twins1.AddTokens(pid1, '$');
            twins2.AddTokens(pid2, '$');

            if (!twins_mode && (twins1.Length() > 1 || twins2.Length() > 1))
               printf("MZ twins compatibility mode activated ...\n\n",
                      twins_mode = true);

            if (!v286)
               for (int id1 = 0; id1 < twins1.Length(); id1++)
                  for (int id2 = 0; id2 < twins2.Length(); id2++)
                     // And output it to the QTDT ibd table
                     fprintf(output, "%8s %8s %8s %8s %8s %8s %8s\n",
                             famid,
                             (const char *) twins1[id1],
                             (const char *) twins2[id2],
                             marker, ibd0, ibd1, ibd2);

            // Repeat the process for additional markers (if any!)
            while (true)
               {
               if (fscanf(f, " %s %*s %s %s %s %*s", marker, ibd0, ibd1, ibd2)!= 4)
                  {
                  done = true;
                  break;
                  }

               if (strcmp(marker, v286 ? "Pair-of-Individuals" : "PAIR-OF-AFFECTEDS") == 0)
                  break;

               if (strcmp(marker, "--") == 0)
                  continue;

               for (int id1 = 0; id1 < twins1.Length(); id1++)
                  for (int id2 = 0; id2 < twins2.Length(); id2++)
                     fprintf(output, "%8s %8s %8s %8s %8s %8s %8s\n",
                         famid,
                         (const char *) twins1[id1],
                         (const char *) twins2[id2],
                         marker, ibd0, ibd1, ibd2);
               }
            };
         fclose(f);
         }
      printf("Processed %d Simwalk2 IBD files ...\n\n", argc - 1);
      }
   else
      {
      printf("Parsing genehunter commands -- EXPERIMENTAL\n\n");

      StringArray input;
      input.Read("genehunter.in");

      StringArray markers;

      for (int i = 0; i < input.Length(); i++)
         {
         StringArray tokens;
         tokens.AddTokens(input[i], WHITESPACE);

         if (tokens.Length() == 0) continue;

         if (SlowCompare(tokens[0], "use") == 0)
            {
            printf("Found USE command\n"
                   "-----------------\n\n"
                   "Marker Order\n");
            markers.Dimension(tokens.Length() >> 1);
            for (int j = 1; j < tokens.Length(); j += 2)
               {
               markers[j >> 1] = tokens[j];
               if (markers[j >> 1].SlowCompare("undrop"))
                  printf("Marker %d - %s\n", j / 2 + 1,
                         (const char *) markers[j >> 1]);
               }
            printf("\n");
            }

         if (SlowCompare(tokens[0], "dump") == 0 && tokens.Length() > 1 &&
             SlowCompare(tokens[1], "ibd") == 0)
            {
            printf("Found DUMP IBD command\n"
                   "----------------------\n\n");

            i++;

            if (i < input.Length())
               {
               if (markers.Length() == 0)
                  error("Did not find USE command to determine marker order");

               if (input[i] == "")
                  input[i] = "ibd_dist.out";

               printf("Reading IBD status from file %s ...\n",
                      (const char *) input[i]);

               StringArray dump;
               dump.Read(input[i]);

               if (dump.Length() == 0)
                  error("IBD file seems to be empty\n");

               StringArray previous;
               previous.AddTokens(dump[0], " \t");

               if (previous.Length() != 12)
                  error("File header doesn't look like a Genehunter2 header");

               for (int j = 1; j < dump.Length(); j++)
                  {
                  StringArray columns;
                  columns.AddTokens(dump[j], WHITESPACE ",");

                  if (columns.Length() == 0)
                     continue;

                  if (columns.Length() != 10)
                     error("Error parsing line %d of IBD file\n", j + 1);

                  bool check = true;

                  if ((j - 1) % markers.Length() == 0)
                     if (previous[2]==columns[2] && previous[3]==columns[3])
                        check = false;

                  if ((j - 1) % markers.Length())
                      if (previous[2]!=columns[2] || previous[3]!=columns[3])
                        check = false;

                  if (!check)
                   error("Prelude is having trouble parsing this file ...\n\n"
                         "Expecting the following layout:\n\n"
                         " * Single Header Line\n"
                         " * IBD at each of %d markers for pair 1\n"
                         " * IBD at each of %d markers for pair 2\n"
                         " * IBD at each of %d markers for ...\n",
                         markers.Length(), markers.Length(), markers.Length());

                  if ((j - 1) % markers.Length() &&
                      (double) previous[0] > (double) columns[0])
                     error("Prelude is having trouble parsing this file ...\n\n"
                           "For each pair of individuals, IBD results at each\n"
                           "marker should appear in map order.\n\n"
                           "However, for pair %s and %s in family %s\n"
                           "IBD at position %s precedes IBD at position %s\n",
                           (const char *) columns[2],
                           (const char *) columns[3],
                           (const char *) columns[1],
                           (const char *) previous[0],
                           (const char *) columns[0]);

                  if (columns[1][0] == '-' || columns[2][0] == '-' ||
                      columns[3][0] == '-' && !negative_probabilities)
                      {
                      negative_probabilities = true;
                      printf("Finale has detected negative IBD probabilities\n"
                             "in your genehunter output ...\n\n");
                      }

                  StringArray twins1, twins2;

                  twins1.AddTokens(columns[2], '$');
                  twins2.AddTokens(columns[3], '$');

                  if (!twins_mode && (twins1.Length()>1 || twins2.Length()>1))
                     printf("MZ twins compatibility mode activated ...\n\n",
                            twins_mode = true);

                  if (markers[(j - 1) % markers.Length()].SlowCompare("undrop"))
                     for (int id1 = 0; id1 < twins1.Length(); id1++)
                        for (int id2 = 0; id2 < twins2.Length(); id2++)
                           fprintf(output, "%s %s %s %s %s %s %s\n",
                             (const char *) columns[1],
                             (const char *) twins1[id1],
                             (const char *) twins2[id2],
                             (const char *) markers[(j - 1) % markers.Length()],
                             (const char *) columns[7],
                             (const char *) columns[8],
                             (const char *) columns[9]);

                  previous.Clear();
                  previous.AddTokens(dump[j], WHITESPACE ",");
                  }
               }
            }
         }
      }

   fclose(output);

   return 0;
   }
 
