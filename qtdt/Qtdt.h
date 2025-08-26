////////////////////////////////////////////////////////////////////// 
// qtdt/Qtdt.h 
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
 
#ifndef _QTDT_H_
#define _QTDT_H_

#include "QtdtResult.h"
#include "QtdtLinear.h"
#include "QtdtBroker.h"

#include "Pedigree.h"
#include "Genetics.h"
#include "MathSVD.h"
#include "MathStats.h"
#include "MathNormal.h"
#include "QuickIndex.h"

class Qtdt : public QtdtResult, public QtdtBroker
   {
   public:
      Matrix      M;
      Vector      y, means, variances;
      char        multiString[500];

      Qtdt();

      // Test a specific trait, marker and allele combination for association
      void TestAllele(Pedigree & ped, int trait, int marker, int al);

      // Get a one letter abbreviation for this model
      char * ModelAbbreviation();

      // Allocate memory required for storing intermediate results
      // that allow us to calculate global permutation p-value
      void SetupPermutations(Pedigree & ped);

      // Use permutation data to assess a significance level globally
      double GlobalSignificance(double p);

      // Calculate an empirical p-value threshold corresponding to adjusted p
      double GlobalThreshold(double p);

   private:
      // These routines implement individual tests
      void linearTestAllele();
      void vcTestAllele();
      void rabTestAllele();
      void zTestAllele();

      // Helpers for TestAllele function
      // Is p an informative child?
      bool isProband(Person & p);
      bool isInformativeTDT(Person & p);

      // Count the degrees of freedom, give the current parameters
      int CountParametersNull();
      int CountParametersFull();

      // Fill in the regression matrix
      void BuildMatrix(Person ** probands, int headCount);
      void BuildY(Person ** probands, int headCount);

      // Wrapper around calls to QtdtLinear which fills in means model
      void LinearDriver(QtdtLinear & engine, bool full = true);

      // Zero degenerate columns, and return their number
      int  EditDegenerates(int first, int last);

      // Clear the TDT engine, before new test
      void Wipe();

      // Find the first sibling that is phenotyped for t and m
      // Returns only independent sib-pairs
      Person * GetSib(Person * p);

      void PrintLinearSummary(SVD & model, Vector & y, double RSS, double varY,
                              int parameters, bool full);
      void PrintNormalSummary(NormalSet & v, bool full);

      void Permute(NormalSet & null, NormalSet & full);
      void Permute(Matrix & M, int firstW);
      void Permute(Vector & scores);

      // This array is used to group informative individuals by family
      IntArray families;

      // This array stores signs for each family
      IntArray flips;

      // Lists all informative probands in an array
      Person ** ListProbands();

      // Array of maximum scores (over all traits and markers) for each permutation
      Vector     bestScores;

      // This array temporarily scores permutated statistics for each analysis
      Vector     pScores;
      QuickIndex pIndex;

      // Summarize permutations
      void UpdatePermutationScores();
   };

#endif
 
