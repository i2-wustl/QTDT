////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtLinear.h 
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
 
#ifndef __QTDTLINEAR_H__
#define __QTDTLINEAR_H__

#include "MathMatrix.h"
#include "Pedigree.h"
#include "Constant.h"
#include "QtdtOptions.h"

class QtdtLinear : protected QtdtOptions
   {
   public:
      int  firstAllele;
      int  lastAllele;
      int  t;
      int  m;

      QtdtLinear()
         { firstAllele = lastAllele = 0;
           t = m = -1; }

      // The labeller, parameter counter and modeller classes
      // override these virtual functions
      // They are as specific and simple as possible to make
      // their implementation unambiguous
      virtual void GrandMean() = 0;

      virtual void Covariates() = 0;
      virtual void Sex() = 0;
      virtual void PaternalPhenotype() = 0;
      virtual void MaternalPhenotype() = 0;

      virtual void W(int al) = 0;
      virtual void B(int al) = 0;
      virtual void Wd(int al) = 0;
      virtual void Bd(int al) = 0;
      virtual void Wi(int al) = 0;
      virtual void Bi(int al) = 0;
      virtual void AllisonMatingTypes(int al) = 0;
      virtual void ImprintingMatingTypes(int al) = 0;
      virtual void X(int al) = 0;
      virtual void Xsq(int al) = 0;
      virtual void Xi(int al) = 0;

      virtual void SibIBD() = 0;
      virtual void SibTrait() = 0;
      virtual void SibIBDTrait() = 0;

      // These functions are composites of the virtual functions,
      // and ensure they are always called in the same order
      void ParentalPhenotypes();

      void AllisonControl();
      void FulkerControl();
      void LinkageControl();
      void AssociationControl();

      void AllisonTest();
      void FulkerTest();
      void LinkageTest();
      void AssociationTest();

      void LinearModel(bool full = true);
   };

class QtdtLinearCounter : public QtdtLinear
   {
   public:
      int count;

      QtdtLinearCounter() : QtdtLinear()
         { count = 0; }

      virtual void GrandMean();

      virtual void Covariates();
      virtual void Sex();
      virtual void PaternalPhenotype();
      virtual void MaternalPhenotype();

      virtual void W(int al);
      virtual void B(int al);
      virtual void Wd(int al);
      virtual void Bd(int al);
      virtual void Wi(int al);
      virtual void Bi(int al);
      virtual void AllisonMatingTypes(int al);
      virtual void ImprintingMatingTypes(int al);
      virtual void X(int al);
      virtual void Xsq(int al);
      virtual void Xi(int al);

      virtual void SibIBD();
      virtual void SibTrait();
      virtual void SibIBDTrait();
   };

class QtdtLinearLabeller : public QtdtLinear
   {
   public:
      char buffer[50];
      int count;
      Matrix * M;

      QtdtLinearLabeller(Matrix * m) : QtdtLinear()
         {
         count = 0;
         M = m;
         }

      virtual void GrandMean();

      virtual void Covariates();
      virtual void Sex();
      virtual void PaternalPhenotype();
      virtual void MaternalPhenotype();

      virtual void W(int al);
      virtual void B(int al);
      virtual void Wd(int al);
      virtual void Bd(int al);
      virtual void Wi(int al);
      virtual void Bi(int al);
      virtual void AllisonMatingTypes(int al);
      virtual void ImprintingMatingTypes(int al);
      virtual void X(int al);
      virtual void Xsq(int al);
      virtual void Xi(int al);

      virtual void SibIBD();
      virtual void SibTrait();
      virtual void SibIBDTrait();
   };

class QtdtLinearModeller : public QtdtLinear
   {
   private:
      bool     cacheFulker, cacheImprinting, cacheIBD;
      double   cacheXi, cachePi;
      char     buffer[BUFSIZE];
      int      column;
      Matrix   cacheB;
      Vector   * row;
      Person   * proband, * father, * mother, * sib;
      Pedigree * ped;

      void IndividualB(Person * p, int al);

   public:
      QtdtLinearModeller() : QtdtLinear(), cacheB("B")
         {
         column = -1 ;
         row = NULL;
         proband = NULL;
         }

      void Select(Pedigree * pedigree, Person * p, Vector * r);
      void SelectSib(Person * p)
         { sib = p; };

      virtual void GrandMean();

      virtual void Covariates();
      virtual void Sex();
      virtual void PaternalPhenotype();
      virtual void MaternalPhenotype();

      virtual void W(int al);
      virtual void B(int al);
      virtual void Wd(int al);
      virtual void Bd(int al);
      virtual void Wi(int al);
      virtual void Bi(int al);
      virtual void AllisonMatingTypes(int al);
      virtual void ImprintingMatingTypes(int al);
      virtual void X(int al);
      virtual void Xsq(int al);
      virtual void Xi(int al);

      virtual void SibIBD();
      virtual void SibTrait();
      virtual void SibIBDTrait();
   };


#endif 
