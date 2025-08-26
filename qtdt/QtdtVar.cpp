////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtVar.cpp 
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
 
#include "QtdtVar.h"
#include "QtdtLinear.h"
#include "IBD.h"

#include <string.h>
#include <math.h>

QtdtVariances::QtdtVariances() : kin()
    {
    marker = trait = allele = -1;
    n = NULL;
    f = NULL;
    probands = NULL;
    vCount = pCount = matrix = 0;
    };

int QtdtVariances::CountVarianceComponents(bool full)
   {
   int temp = full ? varComponentFull : varComponentNull;
   int count = 0;

   while (temp)
      {
      count += temp & 1;
      temp >>= 1;
      }

   return count;
   }

void QtdtVariances::FillMatrices(bool full)
   {
   int temp = full ? varComponentFull : varComponentNull;
   if (temp & VAR_NON_SHARED) FillResiduals();
   if (temp & VAR_COMMON_ENV) FillCommonEnvironment();
   if (temp & VAR_NUCLEAR)    FillNuclearUnit();
   if (temp & VAR_TWIN)       FillTwinEnvironment();
   if (temp & VAR_POLYGENIC)  FillPolygenic();
   if (temp & VAR_ADDITIVE)   FillAdditive();
   if (temp & VAR_DOMINANCE)  FillDominance();
   };

void QtdtVariances::FillResiduals()
   {
   LabelMatrix("Ve");

   n->varComponents[matrix].Zero();

   for (int i = 0; i < pCount; i++)
      n->varComponents[matrix][i][i] = 1.0;

   matrix++;
   };

void QtdtVariances::FillCommonEnvironment()
   {
   LabelMatrix("Vc");

   n->varComponents[matrix].Set(1.0);

   matrix++;
   };

void QtdtVariances::FillPolygenic()
   {
   LabelMatrix("Vg");

   kin.Setup(*f);

   // if (kin.isInbred())
   //   error("Family %s is inbred\n"
   //         "Please break the loop and try again\n", (const char *) f->famid);

   for (int i = 0; i < pCount; i++)
      for (int j = i; j < pCount; j++)
         n->varComponents[matrix][i][j] =
         n->varComponents[matrix][j][i] =
            kin(*(probands[i]), *(probands[j])) * 2;

   matrix++;
   };

void QtdtVariances::FillAdditive()
   {
   LabelMatrix("Va");

   for (int i = 0; i < pCount; i++)
      for (int j = i; j < pCount; j++)
         n->varComponents[matrix][i][j] =
         n->varComponents[matrix][j][i] =
            ibdRepository.Lookup(marker, *probands[i], *probands[j])
            ->expected();

   matrix++;
   };

void QtdtVariances::FillDominance()
   {
   LabelMatrix("Vd");

   for (int i = 0; i < pCount; i++)
      for (int j = i; j < pCount; j++)
         n->varComponents[matrix][i][j] =
         n->varComponents[matrix][j][i] =
            ibdRepository.Lookup(marker, *probands[i], *probands[j])
            -> p2;

   matrix++;
   };

void QtdtVariances::FillTwinEnvironment()
   {
   LabelMatrix("Vd");

   for (int i = 0; i < pCount; i++)
      {
      n->varComponents[matrix][i][i] = 1.0;

      for (int j = i + 1; j < pCount; j++)
         n->varComponents[matrix][i][j] =
         n->varComponents[matrix][j][i] =
            probands[i]->isTwin(*probands[j]);
      }

   matrix++;
   };

void QtdtVariances::FillNuclearUnit()
   {
   LabelMatrix("Vnuclear");

   for (int i = 0; i < pCount; i++)
      for (int j = i; j < pCount; j++)
         {
         double group = 0.0;

         if (i == j) group = 2.0;
         else if (probands[i]->isSib(*(probands[j]))) group = 1.0;
         else if (probands[i]->father == probands[j] ||
                  probands[i]->mother == probands[j] ||
                  probands[j]->father == probands[i] ||
                  probands[j]->mother == probands[i])
             group = 1.0;
         else for (int k = f->founders; k < f->count; k++)
               {
               if (f->ped[f->path[k]].father == probands[i] &&
                   f->ped[f->path[k]].mother == probands[j])
                   group = 1.0;
               if (f->ped[f->path[k]].mother == probands[i] &&
                   f->ped[f->path[k]].father == probands[j])
                   group = 1.0;
               }

         n->varComponents[matrix][i][j] =
         n->varComponents[matrix][j][i] =
            group;
         }

   matrix++;
   };

void QtdtVariances::LabelMatrix(char * name)
   {
   n -> varComponents[matrix].SetLabel(name);
   n -> varComponents[matrix].Dimension(pCount, pCount);

   for (int i = 0; i < pCount; i++)
      {
      n->varComponents[matrix][i].SetLabel(probands[i]->pid);
      n->varComponents[matrix].SetColumnLabel(i, probands[i]->pid);
      }
   }

char * vc_labels[] = {"Ve", "Vc", "Vnuclear", "Vtwin", "Vg", "Va", "Vd" };

char * QtdtVariances::GetLabel(int vc)
   {
   int index = 0;
   while (vc /= 2)
      index++;
   return vc_labels[index];
   }

int QtdtVariances::GetCode(char * label)
   {
   int index = 0;

   while (vc_labels[index] != NULL)
      if (stricmp(label, vc_labels[index]) == 0)
         return 1 << index;
      else
         index++;

   return -1;
   }

int QtdtVariances::GetMatrixIndex(int vc, bool full)
   {
   int temp = full ? varComponentFull : varComponentNull;

   if ((temp & vc) == 0) return -1;

   int index = 0;
   while (vc /= 2)
      if (vc & temp) index++;
   return index;
   }
 
