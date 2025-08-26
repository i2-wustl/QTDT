////////////////////////////////////////////////////////////////////// 
// qtdt/QtdtBroker.cpp 
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
 
#include "QtdtBroker.h"
#include "Error.h"

int QtdtBroker::CountPhenotypes()
   {
   if (associationModel & ZTEST_FAMILY)
      return Person::affectionCount;

   if (associationModel & RABINOWITZ_FAMILY)
      return Person::traitCount + Person::affectionCount;

   return Person::traitCount;
   }

String QtdtBroker::PhenotypeLabel(int t)
   {
   if (associationModel & ZTEST_FAMILY)
      return Person::affectionNames[t] + " [Affecteds]";

   if (t < Person::traitCount)
      return Person::traitNames[t];

   t -= Person::traitCount;

   return Person::affectionNames[t] + " [Binary]";
   }

bool QtdtBroker::isPhenotyped(Person & p, int t)
   {
   if (associationModel & ZTEST_FAMILY)
      return p.affections[t] == 2;

   if (t < Person::traitCount)
      return p.isPhenotyped(t);

   t -= Person::traitCount;

   return p.isDiagnosed(t);
   }

bool QtdtBroker::isFullyPhenotyped(Person & p)
   {
   if (associationModel & ZTEST_FAMILY)
      error("QTDT Phenotype Broker detected inconsistent state");

   if (CountPhenotypes() == Person::traitCount)
      return p.isFullyPhenotyped();

   return p.isFullyPhenotyped() && p.isFullyDiagnosed();
   }

double QtdtBroker::GetPhenotype(Person & p, int t)
   {
   if (associationModel & ZTEST_FAMILY)
      return p.affections[t] == 2;

   if (t < Person::traitCount)
      return p.traits[t];

   t -= Person::traitCount;

   switch (p.affections[t])
      {
      case 1 : return 0.0;
      case 2 : return 1.0;
      default: return _NAN_;
      }
   }

bool QtdtBroker::isBinaryPhenotype(int t)
   {
   return (t >= Person::traitCount || associationModel & ZTEST_FAMILY);
   } 
