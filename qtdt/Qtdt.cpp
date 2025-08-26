////////////////////////////////////////////////////////////////////// 
// qtdt/Qtdt.cpp 
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
 
#include "Qtdt.h"
#include "QtdtVar.h"
#include "Constant.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

Qtdt::Qtdt() : M("Qtdt.linear"), y("Qtdt.y"), means("Qtdt.means"),
               variances("Qtdt.variances")
   { Wipe(); }

void Qtdt::SetupPermutations(Pedigree & p)
   {
   // Flips array must be large enough to accomodate one sign per family
   flips.Dimension(p.familyCount);

   // Best scores array must accomodate the highest score overall
   bestScores.Dimension(mcPermutations);
   bestScores.Set(mcPermutations);

   // pScores array must accomodate the scores for this crop of replicates
   pScores.Dimension(mcPermutations);
   }

void Qtdt::TestAllele(Pedigree & p, int trait, int marker, int allele)
   {
   // Wipe results
   Wipe();

   t = trait;
   m = marker;
   al = allele;
   ped = &p;

   // Reseting the random number generation, ensures that we permute
   // families in exactly the same manner for all traits and markers,
   // and allows us to compute global p-values
   if (mcPermutations) rand.Reset(randomSeed);

   if (varComponentFull)
      vcTestAllele();
   else if (associationModel & RABINOWITZ_FAMILY)
      rabTestAllele();
   else if (associationModel & ZTEST_FAMILY)
      zTestAllele();
   else
      linearTestAllele();
   }

Person ** Qtdt::ListProbands()
   {
   // Allocate an array big enough to include all informative offspring
   Person ** probands = new Person *[ped->count];

   // Families array helps us group informative individuals by family
   families.Dimension(ped->count);

   for (int f=0; f < ped->familyCount; f++)
      for (int i=ped->families[f]->first; i <= ped->families[f]->last; i++)
         if (isProband((*ped)[i]))
            {
            // add to list of probands
            families[probandCount] = f;
            probands[probandCount++] = &((*ped)[i]);
            }

   return probands;
   }

void Qtdt::linearTestAllele()
   {
   Person ** probands = ListProbands();

   if (informativeProbands < minProbands)
      {
      delete [] probands;
      return;
      }

   // Build matrices for regression
   int parametersNull = CountParametersNull();
   int parametersFull = CountParametersFull() - parametersNull;

   M.Dimension(probandCount, parametersFull + parametersNull);
   BuildMatrix(probands, probandCount);

   parametersNull -= EditDegenerates(0, parametersNull);
   parametersFull -= EditDegenerates(0, parametersNull + parametersFull);

   dfNull = probandCount - parametersNull;
   dfFull = dfNull - parametersFull;

   SVD nullModel;
   SVD fullModel;

   nullModel.Decompose(M, probandCount, parametersNull);
   nullModel.Edit();

   M.Dimension(probandCount, M.cols);
   if (nullSamples) Permute(M, parametersNull);

   fullModel.Decompose(M, probandCount, parametersNull + parametersFull);
   fullModel.Edit();

   /* MultivariateY(probands, nullModel, fullModel); */
   BuildY(probands, probandCount);
   nullModel.BackSubst(y);
   fullModel.BackSubst(y);

   variance = y.Var();
   double nullVar = nullModel.RSS(M, y) / (probandCount - 1);
   statisticNull = 1.0 - (nullVar  / (variance + ZEPS)) + ZEPS;

   double RSS = fullModel.RSS(M, y);
   double fullVar = RSS / (probandCount - 1);
   statisticFull = 1.0 - (fullVar  / (variance + ZEPS)) + ZEPS;

   statisticFull = statisticFull > statisticNull ? // eliminate roundoff error,
                   statisticFull : statisticNull;  // numbers as negative -1e-14
                                                   // can cause fdist to abort
                                                   // below

   PrintLinearSummary(nullModel, y, nullVar*(probandCount-1), variance, probandCount - dfNull, false);
   PrintLinearSummary(fullModel, y, fullVar*(probandCount-1), variance, probandCount - dfFull, true);

   variances.Dimension(1);
   variances[0] = fullVar;
   means = fullModel.x;

   if (dfNull != dfFull)
      {
      testValue = (statisticFull - statisticNull) / (dfNull - dfFull);
      testValue /= (1 - statisticFull)/dfFull;

      if (!mcPermutations)
         {
         p = fdist(testValue, dfNull - dfFull, dfFull);
         if (flavour == FLAVOUR_LINKAGE && linkageModel == LINK_SIMPLE)
            p = (fullModel.x[fullModel.x.dim-1] > 0) ? p * 0.5 : 1.0 - p * 0.5;
         }
      else
         {
         int tailCount = 0;

         for (int i = 0; i < mcPermutations; i++)
            {
            Permute(M, parametersNull);
            fullModel.Decompose(M, probandCount, parametersNull+parametersFull);
            fullModel.Edit();
            fullModel.BackSubst(y);

            double score = fullModel.RSS(M, y);
            tailCount += RSS - score > -1e-7;
            pScores[i] = -score;
            }

         p = (double) tailCount / (double) mcPermutations;

         UpdatePermutationScores();
         }
      }
   else
      {
      testValue = 0.0;
      p = 1.0;
      }

   delete [] probands;
   success = 1;
   }

void Qtdt::zTestAllele()
   {
   Person ** probands = ListProbands();

   if (informativeProbands < minProbands)
      {
      delete [] probands;
      return;
      }

   // LinearModeller calculates W scores
   // Terms vector holds the score for each family (proband)
   QtdtLinearModeller helper;
   Vector v(1), terms(probandCount);

   helper.m = m;

   int    lastFamily = families[0];
   int    n = 0, termCount = 0;
   double U = 0;

   dfNull = 0;

   for (int i = 0; i < probandCount; i++)
      {
      if ((associationModel == AFF_TDT_LINKAGE || families[i] != lastFamily)
          && n > 0)
         {
         families[termCount] = lastFamily;
         terms[termCount++] = U / n;
         variance += (U / n) * (U / n);
         dfNull += U ? 1 : 0;
         // N++;

         U = n = 0;
         lastFamily = families[i];
         }
      n++;

      helper.Select(ped, probands[i], &v);
      helper.W(al);

      U += v[0] * GetPhenotype(*(probands[i]), t);
      }

   terms[termCount++] = U / n;
   variance += (U / n) * (U / n);
   // N++;

   double & mean = statisticNull = 0;
   double & stdev = statisticFull;

   terms.Dimension(termCount);
   if (nullSamples) Permute(terms);

   for (int i = 0; i < termCount; i++)
      mean += terms[i];

   stdev = sqrt(variance);

   testValue = fabs(mean) < ZEPS ? 0.0 : mean / stdev;

   if (dfNull <= 0)
      p = 1.0;
   else if (!mcPermutations)
      p = tdist(testValue, dfNull);
   else
      {
      double tailCount = 0;
      double fmean = fabs(mean);

      for (int i = 0; i < mcPermutations; i++)
         {
         Permute(terms);

         double permuted = 0.0;
         for (int j = 0; j < termCount; j++)
            permuted += terms[j];

         double score = fabs(permuted);

         if (fmean - score > -1e-7)
            tailCount ++;

         pScores[i] = score;
         }

      p = (double) tailCount / (double) mcPermutations;

      UpdatePermutationScores();
      }

   delete [] probands;
   success = 1;
   }

void Qtdt::rabTestAllele()
   {
   Person ** probands = ListProbands();

   if (informativeProbands < minProbands)
      {
      delete [] probands;
      return;
      }

   // Score phenotype
   BuildY(probands, probandCount);

   // Calculate after regressing out trait
   Vector residuals;

   // Number of model parameters, used to adjust df
   int parameters = 0;

   // Used to estimate covariate effects and sample mean
   SVD model;

   if (useDeviates)
      // Phenotypes are already deviates [e.g., in some selected samples]
      {
      residuals = y;
      dfNull = probandCount;
      }
   else
      {
      // Build matrices for regression
      parameters = CountParametersFull();

      M.Dimension(probandCount, parameters);
      BuildMatrix(probands, probandCount);

      parameters -= EditDegenerates(0, parameters);

      model.Decompose(M, probandCount, parameters);
      model.Edit();

      model.BackSubst(y);

      // Calculate residuals * w ...
      model.Residuals(M, y, residuals);

      dfNull = probandCount - parameters;
      }

   dfFull = 0;

   Vector terms(probandCount);
   int    termCount = 0;
   double variance = 0;

   if (associationModel == QTDT_RABINOWITZ)
      {
      for (int i = 0; i < probandCount; i++)
         {
         Person * p = probands[i];

         int mother = p->mother->markers[m].countAlleles(al);
         int father = p->father->markers[m].countAlleles(al);
         int child  = p->markers[m].countAlleles(al);

         double delta = (double) child - (mother + father) * 0.5;

         terms[termCount++] = residuals[i] * delta;
         variance += residuals[i] * residuals[i] *
                 ((mother == 1) + (father == 1)) * 0.25;
         }
      }
   else
   // if associationModel is QTDT_MONKS
      {
      QtdtLinearModeller helper;
      Vector v(1);

      helper.m = m;

      int    lastFamily = families[0];
      int    n = 0, N = 0;
      double U = 0;

      dfNull = -parameters;

      for (int i = 0; i < probandCount; i++)
         {
         if (families[i] != lastFamily && n > 0)
            {
            families[termCount] = lastFamily;
            terms[termCount++] = U / n;
            variance += (U / n) * (U / n);
            N++;
            dfNull += U ? 1 : 0;

            U = n = 0;
            lastFamily = families[i];
            }
         n++;

         helper.Select(ped, probands[i], &v);
         helper.W(al);

         U += residuals[i] * v[0];
         }

      terms[termCount++] = U / n;
      variance += (U / n) * (U / n);
      N++;
      }

   double & mean = statisticNull = 0;
   double & stdev = statisticFull;

   terms.Dimension(termCount);
   if (nullSamples) Permute(terms);

   for (int i = 0; i < termCount; i++)
      mean += terms[i];

   stdev = sqrt(variance);

   testValue = fabs(mean) < ZEPS ? 0.0 : mean / stdev;

   if (dfNull <= 0)
      p = 1.0;
   else if (!mcPermutations)
      p = tdist(testValue, dfNull);
   else
      {
      double fmean = fabs(mean);
      int tailCount = 0;

      for (int i = 0; i < mcPermutations; i++)
         {
         Permute(terms);

         double permuted = 0.0;
         for (int j = 0; j < termCount; j++)
            permuted += terms[j];

         double score = fabs(permuted);
         if (fmean - score > -1e-7) tailCount++;
         pScores[i] = score;
         }

      p = (double) tailCount / (double) mcPermutations;

      UpdatePermutationScores();
      }

   if (!useDeviates)
      PrintLinearSummary(model, y, 0, variance, parameters, false);

   delete [] probands;
   success = 1;
   }

void Qtdt::vcTestAllele()
   {
   // The NormalSet class handles the multivariate normal data
   // Optional threading support is off by default
   // To turn on set parameter to NormalSet constructor
   QtdtVariances vc;
   NormalSet * normalNull, * normalFull;

   normalNull = new NormalSet;
   normalFull = new NormalSet;

   int vcNull = vc.CountVarianceComponents(false);
   int vcFull = vc.CountVarianceComponents(true);

   normalNull->Dimension(ped->familyCount, vcNull);
   normalFull->Dimension(ped->familyCount, vcFull);

   // List all informative offspring in array
   Person ** probands = new Person *[ped->count];

   // Index informative families for permutation
   families.Dimension(ped->familyCount);

   // linear parameters
   int linearNull = CountParametersNull();
   int linearFull = CountParametersFull();

   vc.Select(t, m, al);

   // iterate through families
   int index = 0;
   for (int family = 0; family < ped->familyCount; family++)
      {
      int selected = 0;

      for(int i=ped->families[family]->first; i<=ped->families[family]->last; i++)
         if (isProband((*ped)[i]))
            // add to list of probands
            probands[selected++] = &((*ped)[i]);

      if (!selected) continue;

      probandCount += selected;

      families[index] = family;

      M.Dimension(selected, linearFull);
      BuildMatrix(probands, selected);
      (*normalFull)[index].linearModel = M;
      (*normalFull)[index].linearModel.CopyLabels(M);

      M.Dimension(selected, linearNull);
      (*normalNull)[index].linearModel = M;
      (*normalNull)[index].linearModel.CopyLabels(M);

      (*normalNull)[index].scores.Dimension(selected);
      BuildY(probands, selected);
      (*normalFull)[index].scores = (*normalNull)[index].scores = y;

      vc.Select(&(*normalNull)[index], ped->families[family], probands, selected);
      vc.FillMatrices(false);
      vc.Select(&(*normalFull)[index], ped->families[family], probands, selected);
      vc.FillMatrices(true);

      index++;
      }

   delete [] probands;

   if (informativeProbands < minProbands)
      {
      delete normalNull;
      delete normalFull;
      return;
      }

   normalNull->Dimension(index, vcNull);
   normalFull->Dimension(index, vcFull);

   if (nullSamples) Permute(*normalNull, *normalFull);

   normalFull->numericMinimizer = numericMinimizer;
   normalNull->numericMinimizer = numericMinimizer;

   normalNull->DisableConstant();
   normalFull->DisableConstant();

   normalNull->Solve();
   normalFull->Solve();

   normalNull->EnableConstant();
   normalFull->EnableConstant();

   statisticNull = normalNull->Evaluate();
   dfNull = probandCount - normalNull->CountParameters();

   statisticFull = normalFull->Evaluate();
   dfFull = probandCount - normalFull->CountParameters();

   testValue = 2 * (statisticNull - statisticFull);
   if (testValue < 0) testValue = 0;

   PrintNormalSummary(*normalNull, false);
   PrintNormalSummary(*normalFull, true);

   means = normalFull->means;
   variances = normalFull->variances;

   if (mcPermutations == 0)
      p = dfNull > dfFull ? chidist(testValue, dfNull - dfFull) : 1.0;
   else
      {
      int tailCount = 0;

      for (int i = 0; i < mcPermutations; i++)
         {
         Permute(*normalNull, *normalFull);
         normalFull->Solve();

         double score = normalFull->Evaluate();
         tailCount += (statisticFull - score > -1e-7);
         pScores[i] = -score;
         }

      p = (double) tailCount / (double) mcPermutations;

      UpdatePermutationScores();
      }

   delete normalNull;
   delete normalFull;

   success = 1;
   }

void Qtdt::Wipe()
   {
   Clear();
   multiString[0] = 0;
   }

void Qtdt::LinearDriver(QtdtLinear & engine, bool full)
   {
   if (multiAllelic)
      {
      engine.firstAllele = 1;
      engine.lastAllele = ped->CountAlleles(m);
      }
   else
      engine.firstAllele = engine.lastAllele = al;
   engine.t = t;
   engine.m = m;
   engine.LinearModel(full);
   }

int Qtdt::CountParametersNull()
   {
   QtdtLinearCounter counter;
   LinearDriver(counter, false);
   return counter.count;
   };

int Qtdt::CountParametersFull()
   {
   QtdtLinearCounter counter;
   LinearDriver(counter);
   return counter.count;
   }

bool Qtdt::isProband(Person & proband)
   {
   // if this is not a child return
   // the linear models cannot be fitted
   if (associationModel & ~QTDT_ALLOW_FOUNDERS ||
       excludeFounderPhenotypes || linkageModel ||
       covariates & COVAR_PARENTS || imprinting)
       if (proband.isFounder())
        return false;

   Person * father = proband.father;
   Person * mother = proband.mother;

   // have we got the indispensable phenotype information
   if (composite ? !isFullyPhenotyped(proband) : !isPhenotyped(proband, t))
      return false;

   // Genotypes for association
   if (associationModel && !proband.isGenotyped(m)) return false;

   // IBD information for linkage
   if (varComponentFull & (VAR_ADDITIVE | VAR_DOMINANCE) &&
       !ibdRepository.HaveFamily(m, ped->FindFamily(proband.famid)))
       return false;

   // if the current model allows for covariates
   if (~associationModel & ZTEST_FAMILY && !useDeviates)
      // check that we have all necessary covariates
      if ( covariates & COVAR_USER && !proband.isFullyControlled() ||
           covariates & COVAR_SEX && !proband.isSexed() ||
           covariates & COVAR_PARENTS &&
           (!isPhenotyped(*mother, t) || !isPhenotyped(*father, t)))
         return false;

#define ISHETEROZYGOUS(parent) (multiAllelic ? \
                               (parent->markers[m].isHeterozygous()) : \
                               (parent->markers[m].isHeterozygousFor(al)))

   // No imprinting test is possible when both parents are heterozygous
   // and have identical genotypes. If carrying out a disequilibrium test,
   // we need genotypes for both parents. If carrying a total association
   // test, then we can use families where one parent is missing but the
   // parent of interest is homozygous.
   if (imprinting &&
       // Discard if one parental genotype is missing or
       // both parents have identical heterozygous genotypes
       (!father->isGenotyped(m) || !mother->isGenotyped(m) ||
        (ISHETEROZYGOUS(father) && father->markers[m] == mother->markers[m])
       ) &&
       // For population level tests we can recover information
       // from families where one parent is homozygous
       !(associationModel == QTDT_TOTAL &&
         (father->isGenotyped(m) && !ISHETEROZYGOUS(father) ||
          mother->isGenotyped(m) && !ISHETEROZYGOUS(mother))
        )
       )
       return false;

   bool isInformative = flavour == FLAVOUR_VARIANCE;

   switch (associationModel)
      {
      case QTDT_TOTAL :
         // Unless the model includes imprinting, no
         // special conditions need be met...
         isInformative = true;
         break;
      case QTDT_RABINOWITZ :
      case QTDT_ALLISON :
         // Genotyped offspring and parents, one parent must
         // be heterozygous
         if (!mother->isGenotyped(m) ||
             !father->isGenotyped(m) ||
             !mother->markers[m].isHeterozygousFor(al) &&
             !father->markers[m].isHeterozygousFor(al))
             return false;
         isInformative = true;
         break;
      case QTDT_MONKS :
      case AFF_TDT_LINKAGE :
      case AFF_TDT_ZTEST :
      case QTDT_WITHIN :
      case QTDT_STRATIFICATION :
      case QTDT_GENERAL :
         // Can include everyone in the model
         if (transmissionScoring == TRANSMISSION_PEDIGREE)
            {
            isInformative = isInformativeTDT(proband);
            break;
            }
      case QTDT_ORTHOGONAL :
         if (proband.isFounder())
            break;
         if (multiAllelic ?
          (father->markers[m].isHeterozygous() && mother->isGenotyped(m) ||
           mother->markers[m].isHeterozygous() && father->isGenotyped(m)) :
          (father->markers[m].isHeterozygousFor(al) && mother->isGenotyped(m) ||
           mother->markers[m].isHeterozygousFor(al) && father->isGenotyped(m)))
             isInformative = true;
      case QTDT_FULKER :
         if (!isInformative)
            {
            int sibCount = 0;
            for (int i = 0; i < proband.sibCount; i++)
               if (proband.sibs[i]->isGenotyped(m))
                  {
                  sibCount ++;
                  if (multiAllelic ?
                      (proband.markers[m] != proband.sibs[i]->markers[m]) :
                      (proband.markers[m].countAlleles(al) !=
                       proband.sibs[i]->markers[m].countAlleles(al)))
                     isInformative = true;
                  }
            }
         break;
      case QTDT_NONE :
         break;
      default :
         error ("Internal error - invalid association model");
      }

   switch (linkageModel)
      {
      case LINK_SIMPLE :
         // we must have a control sibling
         // and two phenotyped parents -> use external ibds
         if (// !father->isGenotyped(m) ||
             // !mother->isGenotyped(m) ||
             (GetSib(&proband) == NULL))
            return false;
         // if (mother->markers[m].isHeterozygous() ||
         //     father->markers[m].isHeterozygous())
            isInformative = true;
         break;
      case LINK_NONE :
         break;
      default :
         error("Internal error - invalid linkage model");
      }

   // If we are trying to score parental origin of alleles
   // only consider matings where it is unambigous:
   //       * the parental genotypes are different
   //         eg: 12 * 22, 12 * 11, 12 * 13
   if (imprinting && (associationModel & ~QTDT_TOTAL) && (multiAllelic ?
        (imprinting == I_PATERNAL && father->markers[m].isHomozygous() ||
         imprinting == I_MATERNAL && mother->markers[m].isHomozygous()) :
        (imprinting == I_PATERNAL && father->markers[m].isHomozygousFor(al) ||
         imprinting == I_MATERNAL && mother->markers[m].isHomozygousFor(al))))
      return false;

   informativeProbands += isInformative;

   return true;
   }

bool Qtdt::isInformativeTDT(Person & proband)
   {
   if (proband.isFounder()) return false;

   int    * pedinfo = new int [proband.traverse + 1];
   Family * f = ped->FindFamily(proband.famid);

   for (int i = 0; i <= proband.traverse; i++)
      {
      // Untyped individuals as -1
      // Heterozygotes are marked as 1
      // Informative individuals as 3
      // Homozygotes are labelled with even number describing allele
      int fat, mot;
      Person * p = ped->persons[f->path[i]];

      if (!p->isFounder() && (fat = pedinfo[p->father->traverse]) != -1 &&
                             (mot = pedinfo[p->mother->traverse]) != -1 )
         pedinfo[i] = (fat & 1 || mot & 1) ? 3 : (fat == mot ? fat : 1);
      else if (p->isGenotyped(m))
         {
         pedinfo[i] = !multiAllelic ? p->markers[m].countAlleles(al) :
                      p->markers[m].isHeterozygous() ? 1 : 2*p->markers[m].one;
         for (int j = 0; j < p->sibCount; j++)
            if (p->sibs[j]->markers[m].isKnown())
               if (!multiAllelic ?
                    (pedinfo[i] != p->sibs[j]->markers[m].countAlleles(al)) :
                    (p->markers[m] == p->sibs[j]->markers[m]) )
                  { pedinfo[i] = 3; break; }
         }
      else
         pedinfo[i] = -1;
      }

   bool answer = pedinfo[proband.traverse] == 3;
   delete [] pedinfo;

   return answer;
   }

/*
void Qtdt::MultivariateY(Person ** probands, SVD & null, SVD & full)
   {
   if (!composite)
      BuildY(probands, probandCount);
   else
      {
      y.Dimension(probandCount);
      MultiVariateHelper mvh;
      mvh.Setup(&null, &full, M, probands);
      mvh.Optimize();
      mvh.RetrieveY(y);
      mvh.DescribeY(multiString);
      }
   }
*/

void Qtdt::BuildY(Person ** probands, int headCount)
   {
   y.Dimension(headCount);

   for (int i = 0; i < headCount; i++)
      y[i] = GetPhenotype(*(probands[i]), t);;
   }

void Qtdt::BuildMatrix(Person ** probands, int headCount)
   {
   QtdtLinearLabeller label(&M);
   LinearDriver(label);

   QtdtLinearModeller engine;

   for (int i = 0; i < headCount; i++)
      {
      Person * p = probands[i];
      Person * sib = (linkageModel == LINK_SIMPLE) ? GetSib(p) : NULL;

      engine.Select(ped, p, &(M[i]));
      engine.SelectSib(sib);

      LinearDriver(engine);
      }
   };

int Qtdt::EditDegenerates(int first, int last)
   {
   SVD engine;

   int i, col;
   for (i = col = first; i < last; i++)
      {
      engine.Decompose(M, probandCount, col + 1);
      engine.Edit();

      bool redundant = false;
      for (int j = 0; j <= col; j++)
         if (engine.w[j] == 0.0)
            redundant = true;

      if (redundant)
         M.DeleteColumn(col);
      else
         col++;
      }
   return last - col;
   };

char * Qtdt::ModelAbbreviation()
   {
   if (flavour == FLAVOUR_VARIANCE) return "V";
   if (flavour == FLAVOUR_LINKAGE)  return "L";
   if (associationModel == QTDT_STRATIFICATION) return "S";
   if (associationModel == QTDT_TOTAL) return "X";
   if (imprinting == I_FULL) return "TI";
   if (imprinting == I_IMPRINTING) return "I";
   return "T";
   }

Person * Qtdt::GetSib(Person * p)
   {
   for (int i = 0; i < p->sibCount; i++)
      {
      Person * sib = p->sibs[i];

      if (sib->isGenotyped(m) && ( composite ? isFullyPhenotyped(*sib) :
                                               isPhenotyped(*sib, t)))
         return (sib == p) ? NULL : sib;
      }
   return NULL;
   }

void Qtdt::PrintLinearSummary(SVD & model, Vector & y, double RSS, double varY,
                              int parameters, bool full)
   {
   if (rawDataFile == NULL) return;

   int df = full ? dfFull : dfNull;

   if (!full)
      {
      fprintf(rawDataFile, "Trait: %-15s ",
              traitWise ? (const char *) PhenotypeLabel(t) : "Composite");
      if (markerWise)
         fprintf(rawDataFile, "Marker: %-15s ",
                 (const char *) Person::markerNames[m]);
      if (alleleWise)
         fprintf(rawDataFile, "Allele: %s",
                 (const char *) Person::GetMarkerInfo(m)->GetAlleleLabel(al));
      fprintf(rawDataFile, "\n=============================="
                           "=============================\n");
      fprintf(rawDataFile, "Total Probands: %d\n", probandCount);
      fprintf(rawDataFile, "Var(Y):         %.2f\n", varY);
      fprintf(rawDataFile, "First %d Y values - ", y.dim > 10 ? 10 : y.dim);
      y.Print(rawDataFile, 10);
      fprintf(rawDataFile, "\n");
      }

   fprintf(rawDataFile, "%s HYPOTHESIS\n---------------\n\n",
                        full ? "FULL" : "NULL");

   fprintf(rawDataFile, "First %d lines of regression matrix...",
           M.rows > 10 ? 10 : M.rows);
   M.Print(rawDataFile, 10, parameters);

   fprintf(rawDataFile, "\nSome useful information...\n"
                        "               df : %d\n"
                        "              RSS : %.2f\n", df, RSS);

   fprintf(rawDataFile, "Estimates - ");
   model.x.Print(rawDataFile);

   fprintf(rawDataFile, full ? "\n\n" : "\n");
   }

void Qtdt::PrintNormalSummary(NormalSet & normal, bool full)
   {
   if (rawDataFile == NULL) return;

   int df  = full ? dfFull : dfNull;

   Vector & y = normal[0].scores;
   Matrix & M = normal[0].linearModel;

   if (!full)
      {
      fprintf(rawDataFile, "Trait: %-15s",
              traitWise ? (const char *) PhenotypeLabel(t) : "Composite");
      if (markerWise)
         fprintf(rawDataFile, "Marker: %-15s ",
                              (const char *) Person::markerNames[m]);
      if (alleleWise)
         fprintf(rawDataFile, "Allele: %s",
                 (const char *) Person::GetMarkerInfo(m)->GetAlleleLabel(al));
      fprintf(rawDataFile, "\n================================"
                           "================================\n");
      fprintf(rawDataFile, "Total Probands: %5d\n", probandCount);
      fprintf(rawDataFile, "Family #1 Phenotypes - ");
      y.Print(rawDataFile, 10);
      fprintf(rawDataFile, "\n");
      }

   fprintf(rawDataFile, "%s HYPOTHESIS\n---------------\n\n",
           full ? "FULL" : "NULL");

   fprintf(rawDataFile, "Family #1 var-covar matrix terms [%d]...",
           normal.variances.dim);
   for(int i = 0; i < normal.variances.dim; i++)
#ifdef _UNSAFE_
      normal[0].varComponents[i].Print(rawDataFile, 10, 10);
#else
      fprintf(rawDataFile, "[%s]", (const char *) normal[0].varComponents[i].label);
#endif

   fprintf(rawDataFile, "\nFamily #1 regression matrix...",
           M.rows > 10 ? 10 : M.rows);
   M.Print(rawDataFile, 10, probandCount - df);

   fprintf(rawDataFile, "\nSome useful information...\n"
                        "               df : %d\n"
                        "  log(likelihood) : %.2f\n", df, normal.likelihood);

   fprintf(rawDataFile, "        ");
   normal.variances.Print(rawDataFile);
   fprintf(rawDataFile, "            ");
   normal.means.Print(rawDataFile);

   if (diagnosticQ)
      {
      fprintf(rawDataFile, "Diagnostic Statistics for all families\n"
                           "======================================\n\n"
                           "Individual  Residual  Q~ChiSq   "
                           "FQ~ChiSq  FQ~N(0,1)\n");

      double sum = 0.0;

      for (int j = 0; j < normal.count; j++)
         {
         int i;

         normal[j].Diagnostics();

         sum += normal[j].rawQ;

         for (i = 1; i < normal[j].residuals.dim; i++)
            fprintf(rawDataFile, "%10s  %8.3lf  %8.3lf\n",
               (const char *) normal[j].linearModel[i].label,
               normal[j].residuals[i],
               normal[j].Qi[i]);

         fprintf(rawDataFile, "%10s  %8.3lf  %8.3lf  %8.3lf  %8.3lf\n",
               (const char *) normal[j].linearModel[i].label,
               normal[j].residuals[i],
               normal[j].Qi[i],
               normal[j].rawQ,
               normal[j].Q);
         }
      fprintf(rawDataFile, "%30s  %8.3lf", "Grand Total", sum);
      }

   fprintf(rawDataFile, full ? "\n\n" : "\n");
   };

void Qtdt::Permute(Vector & scores)
   {
   // TO DO -- Use Mersenne Twister for Prime Generation

   // Select a set of random signs (one per family)
   for (int i = 0; i < flips.Length(); i++)
      flips[i] = rand.Binary();

   // Permute scores for each family as appropriate
   for (int i = 0; i < scores.dim; i++)
      if (flips[families[i]])
         scores[i] *= -1;
   };

void Qtdt::Permute(NormalSet & normalNull, NormalSet & normalFull)
   {
   int firstW = normalNull[0].linearModel.cols;
   int lastW = normalFull[0].linearModel.cols;

   // Select a set of random signs (one per family)
   for (int i = 0; i < flips.Length(); i++)
      flips[i] = rand.Binary();

   // permute the w scores in each family
   for (int j = 0; j < normalFull.count; j++)
      if (flips[families[j]])
         for (int k = 0; k < normalFull[j].linearModel.rows; k++)
            for (int l = firstW; l < lastW; l++)
               normalFull[j].linearModel[k][l] *= -1;
   };

void Qtdt::Permute(Matrix & M, int firstW)
   {
   int lastW = M.cols;

   // Select a set of random signs (one per family)
   for (int i = 0; i < flips.Length(); i++)
      flips[i] = rand.Binary();

   // Permute scores for each family as appropriate
   for (int i = 0; i < M.rows; i++)
      if (flips[families[i]])
         for (int j = firstW; j < lastW; j++)
            M[i][j] *= -1;
   };

void Qtdt::UpdatePermutationScores()
   {
   pIndex.Index(pScores);

   int rank = mcPermutations;
   double last = pScores[pIndex[0]];

   for (int i = 0; i < mcPermutations; i++)
      {
      if (pScores[pIndex[i]] > last)
         {
         last = pScores[pIndex[i]];
         rank = mcPermutations - i;
         }
      pScores[pIndex[i]] = rank;
      }

   for (int i = 0; i < mcPermutations; i++)
      if (pScores[i] < bestScores[i])
         bestScores[i] = pScores[i];
   }

double Qtdt::GlobalSignificance(double p)
   {
   return (mcPermutations - bestScores.CountIfGreater(p * mcPermutations + 0.5))
        / (double) mcPermutations;
   }

double Qtdt::GlobalThreshold(double p)
   {
   pIndex.Index(bestScores);

   return bestScores[pIndex[p]] / (double) mcPermutations;
   }
 
