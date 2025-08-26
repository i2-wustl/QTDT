////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtLinear.cpp 
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
 
#include "QtdtLinear.h"
#include "QtdtBroker.h"

#include "Genetics.h"
#include "IBD.h"
#include "MathVector.h"

#include <string.h>
#include <math.h>

#define B_ADD     0
#define B_DOM     1
#define B_IMP     2

void QtdtLinear::ParentalPhenotypes()
   {
   PaternalPhenotype();
   MaternalPhenotype();
   }

void QtdtLinear::AllisonControl()
   {
   if (!associationModel) return;
   if (imprinting)
      {
      ImprintingMatingTypes(firstAllele);
      if (imprinting == I_IMPRINTING)
         {
         X(firstAllele);
         if (dominance) Xsq(firstAllele);
         }
      }
   else
      AllisonMatingTypes(firstAllele);
   }

void QtdtLinear::AllisonTest()
   {
   for (int i = firstAllele; i <= lastAllele; i++)
      {
      if (imprinting != I_PATERNAL &&
          imprinting != I_MATERNAL &&
          imprinting != I_IMPRINTING)
         {
         X(i);
         if (dominance) Xsq(i);
         }
      if (imprinting) Xi(i);
      }
   }

void QtdtLinear::FulkerTest()
   {
   for (int a = firstAllele; a <= lastAllele; a++)
      {
      if (imprinting != I_PATERNAL &&
          imprinting != I_MATERNAL &&
          imprinting != I_IMPRINTING)
          {
          W(a);
          if (dominance) Wd(a);
          }
      if (imprinting) Wi(a);
      }
   }

void QtdtLinear::FulkerControl()
   {
   for (int a = firstAllele; a <= lastAllele; a++)
      {
      B(a);
      if (dominance) Bd(a);
      if (imprinting) Bi(a);
      if (imprinting == I_IMPRINTING)
          {
          W(a);
          if (dominance) Wd(a);
          }
      }
   }

void QtdtLinear::LinkageControl()
   {
   if (!linkageModel) return;
   SibIBD();
   SibTrait();
   }

void QtdtLinear::LinkageTest()
   {
   if (linkageModel == LINK_SIMPLE)
      SibIBDTrait();
   }

void QtdtLinear::AssociationControl()
   {
   switch (associationModel)
      {
      case QTDT_ALLISON :
         AllisonControl();
         break;
      case QTDT_FULKER :
      case QTDT_ORTHOGONAL :
      case QTDT_GENERAL :
         FulkerControl();
         break;
      case QTDT_STRATIFICATION :
         AllisonTest();
         break;
      case QTDT_TOTAL :
         if (imprinting == I_IMPRINTING) X(firstAllele);
      case QTDT_WITHIN :
      case QTDT_MONKS :
      case QTDT_RABINOWITZ :
         break;
      }
   }

void QtdtLinear::AssociationTest()
   {
   switch (associationModel)
      {
      case QTDT_TOTAL :
      case QTDT_ALLISON :
         AllisonTest();
         break;
      case QTDT_FULKER :
      case QTDT_ORTHOGONAL :
      case QTDT_GENERAL :
      case QTDT_WITHIN :
      case QTDT_STRATIFICATION :
         FulkerTest();
         break;
      case QTDT_MONKS :
      case QTDT_RABINOWITZ :
         break;
      }
   }

void QtdtLinear::LinearModel(bool full)
   {
   if (associationModel != QTDT_ALLISON)
      GrandMean();

   if (covariates & COVAR_USER)    Covariates();
   if (covariates & COVAR_PARENTS) ParentalPhenotypes();
   if (covariates & COVAR_SEX)     Sex();

   AssociationControl();
   LinkageControl();

   switch(flavour) {
      case FLAVOUR_ASSOC :
         LinkageTest();
         if (!full) return;
         AssociationTest();
         break;

      case FLAVOUR_LINKAGE :
         AssociationTest();
         if (!full) return;
         LinkageTest();
         break;

      case FLAVOUR_VARIANCE :
         AssociationTest();
         LinkageTest();
      }
   }

// QTDT Linear Counter class
//

void QtdtLinearCounter::GrandMean()
   { count ++; }

void QtdtLinearCounter::Covariates()
   { count += Person::covariateCount; }

void QtdtLinearCounter::Sex()
   { count ++; }

void QtdtLinearCounter::PaternalPhenotype()
   { count ++; }

void QtdtLinearCounter::MaternalPhenotype()
   { count ++; }

void QtdtLinearCounter::W(int )
   { count ++; }

void QtdtLinearCounter::B(int )
   { count ++; }

void QtdtLinearCounter::Wd(int )
   { count ++; }

void QtdtLinearCounter::Bd(int )
   { count ++; }

void QtdtLinearCounter::Wi(int )
   { count ++; }

void QtdtLinearCounter::Bi(int )
   { count ++; }

void QtdtLinearCounter::AllisonMatingTypes(int )
   { count += 3; }

void QtdtLinearCounter::ImprintingMatingTypes(int )
   { count += 5; }

void QtdtLinearCounter::X(int )
   { count ++; }

void QtdtLinearCounter::Xsq(int )
   { count ++; }

void QtdtLinearCounter::Xi(int )
   { count ++; }

void QtdtLinearCounter::SibIBD()
   { count ++; }

void QtdtLinearCounter::SibTrait()
   { count ++; }

void QtdtLinearCounter::SibIBDTrait()
   { count ++; }

// QTDT Linear Labeller class
//

void QtdtLinearLabeller::GrandMean()
   {
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "Mu");
   }

void QtdtLinearLabeller::Covariates()
   {
   for (int i = 0; i < Person::covariateCount; i++)
      {
      M->SetColumnLabel(count + i, Person::covariateNames[i]);
      M->SetColWidth(i, 7);
      M->SetColPrecision(i, 2);
      }
   count += Person::covariateCount;
   }

void QtdtLinearLabeller::Sex()
   {
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "sex");
   }

void QtdtLinearLabeller::PaternalPhenotype()
   {
   M->SetColWidth(count, 7);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, "father_y");
   }

void QtdtLinearLabeller::MaternalPhenotype()
   {
   M->SetColWidth(count, 7);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, "mother_y");
   }

void QtdtLinearLabeller::W(int al)
   {
   if (firstAllele != lastAllele)
      sprintf(buffer, "W(%s)", (const char *)
              Pedigree::GetMarkerInfo(m)->GetAlleleLabel(al));
   else
      sprintf(buffer, "W");

   M->SetColWidth(count, 5);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, buffer);
   }

void QtdtLinearLabeller::B(int al)
   {
   if (firstAllele != lastAllele)
      sprintf(buffer, "B(%s)", (const char *)
              Pedigree::GetMarkerInfo(m)->GetAlleleLabel(al));
   else
      sprintf(buffer, "B");

   M->SetColWidth(count, 5);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, buffer);
   }

void QtdtLinearLabeller::Wd(int )
   {
   M->SetColWidth(count, 5);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, "Wd");
   }

void QtdtLinearLabeller::Bd(int )
   {
   M->SetColWidth(count, 5);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, "Bd");
   }

void QtdtLinearLabeller::Wi(int al)
   {
   if (firstAllele != lastAllele)
      sprintf(buffer, imprinting == I_PATERNAL ? "Wpat(%s)" : "Wmat(%s)",
              (const char *) Pedigree::GetMarkerInfo(m)->GetAlleleLabel(al));
   else
      sprintf(buffer, imprinting == I_PATERNAL ? "Wpat" : "Wmat" );

   M->SetColWidth(count, 5);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, buffer);
   }

void QtdtLinearLabeller::Bi(int al)
   {
   if (firstAllele != lastAllele)
      sprintf(buffer, imprinting == I_PATERNAL ? "Bpat(%s)" : "Bmat(%s)",
              (const char *) Pedigree::GetMarkerInfo(m)->GetAlleleLabel(al));
   else
      sprintf(buffer, imprinting == I_PATERNAL ? "Bpat" : "Bmat");

   M->SetColWidth(count, 5);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, buffer);
   }

void QtdtLinearLabeller::AllisonMatingTypes(int )
   {
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "Aa*AA");
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "Aa*Aa");
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "Aa*aa");
   }

void QtdtLinearLabeller::ImprintingMatingTypes(int )
   {
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "mot_aA*AA");
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "mot_aA*aA");
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "mot_aA*aa");
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "mot_AA*Aa");
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "mot_aa*Aa");
   }

void QtdtLinearLabeller::X(int al )
   {
   if (firstAllele != lastAllele)
      sprintf(buffer, "X(%s)",
              (const char *) Pedigree::GetMarkerInfo(m)->GetAlleleLabel(al));
   else
      sprintf(buffer, "X");

   M->SetColWidth(count, 3);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, buffer);
   }

void QtdtLinearLabeller::Xsq(int )
   {
   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, "X*X");
   }

void QtdtLinearLabeller::Xi(int al )
   {
   if (firstAllele != lastAllele)
      sprintf(buffer, imprinting == I_PATERNAL ? "Xpat(%s)" : "Xmat(%s)",
              (const char *) Pedigree::GetMarkerInfo(m)->GetAlleleLabel(al));
   else
      sprintf(buffer, imprinting == I_PATERNAL ? "Xpat" : "Xmat");

   M->SetColWidth(count, 1);
   M->SetColPrecision(count, 0);
   M->SetColumnLabel(count++, buffer);
   }

void QtdtLinearLabeller::SibIBD()
   {
   M->SetColWidth(count, 4);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, "Pi");
   }

void QtdtLinearLabeller::SibTrait()
   {
   M->SetColWidth(count, 7);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, "T");
   }

void QtdtLinearLabeller::SibIBDTrait()
   {
   M->SetColWidth(count, 7);
   M->SetColPrecision(count, 2);
   M->SetColumnLabel(count++, "Pi*T");
   }

// QTDT linear modeller
//

void QtdtLinearModeller::Select(Pedigree * pedigree, Person * p, Vector * r)
   {
   if (proband == NULL || (strcmp(proband->famid, p->famid) != 0))
      cacheFulker = false;

   row = r;
   proband = p;
   ped = pedigree;
   column = 0;

   father = proband->father;
   mother = proband->mother;
   sib = NULL;

   cacheImprinting = cacheIBD = false;

   row->SetLabel((const char *) (proband->famid + '.' + proband->pid));
   }

void QtdtLinearModeller::GrandMean()
   {
   (*row)[column ++] = 1.0;
   }

void QtdtLinearModeller::Covariates()
   {
   for (int i = 0; i < Person::covariateCount; i++)
      (*row)[column ++] = proband->covariates[i];
   }

void QtdtLinearModeller::Sex()
   {
   (*row)[column ++] = proband->sex;
   }

void QtdtLinearModeller::PaternalPhenotype()
   {
   (*row)[column ++] = QtdtBroker::GetPhenotype(*father, t);
   }

void QtdtLinearModeller::MaternalPhenotype()
   {
   (*row)[column ++] = QtdtBroker::GetPhenotype(*mother, t);
   }

void QtdtLinearModeller::IndividualB(Person * p, int al)
   {
   if (!p->isFounder())
      {
      if (!imprinting && transmissionScoring == TRANSMISSION_PEDIGREE)
         {
         // use all available information in the pedigree
         double fatadd = cacheB[p->father->traverse][B_ADD];
         double motadd = cacheB[p->mother->traverse][B_ADD];
         double fatdom = cacheB[p->father->traverse][B_DOM];
         double motdom = cacheB[p->mother->traverse][B_DOM];

         if ((motadd != _NAN_) && (fatadd != _NAN_))
            {
            double fat11 = (1 - fatdom - fatadd) * 0.5, fat22 = fatadd + fat11;
            double mot11 = (1 - motdom - motadd) * 0.5, mot22 = motadd + mot11;

            cacheB[p->traverse][B_ADD] = (fatadd + motadd) * 0.5;
            cacheB[p->traverse][B_DOM] = fat11*mot22 + fat22*mot11 + 0.5 *
               (fatdom*(1 - motdom) + (1 - fatdom)*motdom + fatdom*motdom);
            cacheB[p->traverse][B_IMP] = _NAN_;
            return;
            }
         }
      else if (associationModel != QTDT_FULKER && p->isGenotyped(m) &&
               p->mother->isGenotyped(m) && p->father->isGenotyped(m))
         {
         // Use the parental average
         bool a = p->father->markers[m].one == al,
              b = p->father->markers[m].two == al,
              c = p->mother->markers[m].one == al,
              d = p->mother->markers[m].two == al;
         cacheB[p->traverse][B_ADD] =((a+c) + (a+d) + (b+c) + (b+d))*0.25 - 1;
         cacheB[p->traverse][B_DOM] =((a^c) + (a^d) + (b^c) + (b^d))*0.25;
         cacheB[p->traverse][B_IMP] =(imprinting == I_PATERNAL) ?
                                       ((a|0) + (a|0) + (b|0) + (b|0))*0.25 :
                                       ((0|c) + (0|d) + (0|c) + (0|d))*0.25;
         return;
         }

      // Use the sibling average
      double add = 0, dom = 0;
      int    sibs = 0;

      for (int i = 0; i < p->sibCount; i++)
         {
         if (p->sibs[i]->isGenotyped(m))
            {
            bool a = p->sibs[i]->markers[m].one == al,
                 b = p->sibs[i]->markers[m].two == al;

            add += a + b;
            dom += a ^ b;
            sibs++;
            }
         }

      if (sibs)
         {
         add = add / sibs - 1;
         dom = dom / sibs;

         cacheB[p->traverse][B_ADD] = add;
         cacheB[p->traverse][B_DOM] = dom;
         cacheB[p->traverse][B_IMP] = _NAN_;
         return;
         }
      }

   // Use self genotype
   if (p->isGenotyped(m))
      {
      bool a = p->markers[m].one == al,
           b = p->markers[m].two == al;

      cacheB[p->traverse][B_ADD] = a + b - 1;
      cacheB[p->traverse][B_DOM] = a ^ b;
      cacheB[p->traverse][B_IMP] = _NAN_;
      return;
      }

   // Bail-out
   cacheB[p->traverse].Set(_NAN_);
   }

void QtdtLinearModeller::B(int al)
   {
   if (!cacheFulker)
      {
      Family * f = ped->FindFamily(proband->famid);
      cacheB.Dimension(f->count, 3);
      for (int i = 0; i < f->count; i++)
         IndividualB(ped->persons[f->path[i]], al);
      cacheFulker = firstAllele == lastAllele;
      }
   (*row)[column ++] = cacheB[proband->traverse][B_ADD];
   }

void QtdtLinearModeller::W(int al )
   {
   // after B() calculates the between means, and
   // stores them in the cache variables it is easy
   B(al);
   (*row)[column - 1] = proband->markers[m].countAlleles(al) - 1
                        - cacheB[proband->traverse][B_ADD];
   }

void QtdtLinearModeller::Bd(int al )
   {
   // after B() calculates the between means, and
   // stores them in the cache variables it is easy
   B(al);
   (*row)[column - 1] = cacheB[proband->traverse][B_DOM];
   }

void QtdtLinearModeller::Wd(int al )
   {
   // after B() calculates the between means, and
   // stores them in the cache variables it is easy
   B(al);
   (*row)[column - 1] = (proband->markers[m].countAlleles(al) == 1)
                      - cacheB[proband->traverse][B_DOM];
   }

void QtdtLinearModeller::Bi(int al)
   {
   // after B() calculates the between means, and
   // stores them in the cache variables it is easy
   B(al);
   (*row)[column - 1] = cacheB[proband->traverse][B_IMP];
   }

void QtdtLinearModeller::Wi(int al)
   {
   // after B() calculates the between means, and
   // Xi() calculates the probands genotype it is easy
   B(al);
   column--;
   Xi(al);
   (*row)[column - 1] = cacheXi - cacheB[proband->traverse][B_IMP];
   }

void QtdtLinearModeller::AllisonMatingTypes(int al)
      {
      // the next three columns encode unordered mating type,
      // as follows (little a is the test allele):
      // Col     1    2    3
      // aAxAA   1    0    0
      // aAxaA   0    1    0
      // aAxaa   0    0    1

      int aCount = mother->markers[m].countAlleles(al) +
                   father->markers[m].countAlleles(al);

      (*row)[column ++] = aCount == 1;
      (*row)[column ++] = aCount == 2;
      (*row)[column ++] = aCount == 3;
      }

void QtdtLinearModeller::ImprintingMatingTypes(int al)
      {
      // the next three columns encode unordered mating type,
      // as follows (little a is the test allele):
      // Col     1    2    3
      // aAxAA   1    0    0
      // aAxaA   0    1    0
      // aAxaa   0    0    1

      int mCount = mother->markers[m].countAlleles(al);
      int fCount = father->markers[m].countAlleles(al);

      (*row)[column ++] = mCount == 1 && fCount == 0;
      (*row)[column ++] = mCount == 1 && fCount == 1;
      (*row)[column ++] = mCount == 1 && fCount == 2;
      (*row)[column ++] = mCount == 0 && fCount == 1;
      (*row)[column ++] = mCount == 2 && fCount == 1;
      }

void QtdtLinearModeller::X(int al)
   {
   double X = proband->markers[m].countAlleles(al) - 1;
   (*row)[column++] = X;
   }

void QtdtLinearModeller::Xsq(int al )
   {
   double X = proband->markers[m].countAlleles(al) - 1;
   (*row)[column++] = X * X;
   }

void QtdtLinearModeller::Xi(int al )
   {
   if (father->markers[m] == mother->markers[m] && father->markers[m].isHeterozygousFor(al))
      printf("PROBLEM");

   // For homozygotes, maternal allele is ...
   if (proband->markers[m].isHomozygousFor(al))
      cacheXi = proband->markers[m].one == al ? 1 : 0;
   // For individuals with an homozygous mother, maternal allele is ...
   else if (mother->markers[m].isKnown() && mother->markers[m].isHomozygousFor(al))
      cacheXi = mother->markers[m].one == al;
   // For individuals with an homozygous father, maternal allele is ...
   else if (father->markers[m].isKnown() && father->markers[m].isHomozygousFor(al))
      cacheXi = proband->markers[m].countAlleles(al) -
           (father->markers[m].one == al);
   // For all other individuals, we assume that parents have different
   // heterozygous genotypes (enforced in the isProband() function)
   else
      cacheXi = father->markers[m].hasAllele(
               proband->markers[m].otherAllele(al));

   // Reverse maternal vs. paternal imprinting
   if (imprinting == I_PATERNAL)
      cacheXi = proband->markers[m].countAlleles(al) - cacheXi;

   // Finito!
   (*row)[column ++] = cacheXi;
   }

void QtdtLinearModeller::SibIBD()
   {
   if (!cacheIBD)
      {
      IBD ibd;
      cachePi = ibdRepository.Lookup(m, *proband, *sib)->expected();
      }
   (*row)[column ++] = cachePi;
   }

void QtdtLinearModeller::SibTrait()
   {
   (*row)[column ++] = QtdtBroker::GetPhenotype(*sib, t);
   }

void QtdtLinearModeller::SibIBDTrait()
   {
   SibIBD();
   (*row)[column - 1] = cachePi * QtdtBroker::GetPhenotype(*sib, t);
   }

 
