////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtVar.h 
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
 
#ifndef __QTDTVAR_H__
#define __QTDTVAR_H__

#include "QtdtOptions.h"
#include "MathNormal.h"
#include "Pedigree.h"
#include "Kinship.h"

class QtdtVariances : protected QtdtOptions
   {
   public:
      NormalEquations * n;
      Family          * f;
      Person         ** probands;
      int               vCount, pCount;
      int               matrix;
      int               trait, marker, allele;

      QtdtVariances();

      int  CountVarianceComponents(bool full);

      void Select(int t, int m, int al)
         {
         trait = t;
         marker = m;
         allele = al;
         }

      void Select(NormalEquations * ne,
                  Family * fam,
                  Person ** people, int peopleCount)
         {
         n = ne;
         f = fam;
         probands = people;
         pCount = peopleCount;
         matrix = 0;
         }

      void FillMatrices(bool full);

      void FillResiduals();
      void FillCommonEnvironment();
      void FillPolygenic();
      void FillNuclearUnit();
      void FillAdditive();
      void FillDominance();
      void FillGenotypeInteraction();
      void FillTwinEnvironment();

      void LabelMatrix(char * name);

      // Utility functions to allow identification of named
      // variance components
      static char * GetLabel(int vc_code);
      static int    GetCode(char * vc_label);
      static int    GetMatrixIndex(int vc, bool full);

   private:
      Kinship kin;
   };

#endif 
