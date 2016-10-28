/*
  by Didem Unat 
  */

#ifndef _READWRITEPROP_H
#define _READWRITEPROP_H

//#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>

#include "../datatypes.h"


using namespace std;
using namespace SageInterface;

class ReadWriteProp
{
public:
 void run(SgNode* node);
 
 ReadWriteProp(){}

 void getReadOnlyVars(set<NameComp_t>&); //{ return readOnlyVars;}
 void getWriteOnlyVars(set<NameComp_t>&); //{ return writeOnlyVars;}
 void getBothReadWrittenVars(set<NameComp_t> &); //{ return bothReadWrittenVars;}
 //set<NameComp_t> getLoopIndependentVars(){ return loopIndepVars;} //returns arrays whose's indices are not dependent on loop indices

 private:

 set<NameComp_t> readVars; 
 set<NameComp_t> readOnlyVars; 
 set<NameComp_t> writeVars; 
 set<NameComp_t> writeOnlyVars; 
 set<NameComp_t> bothReadWrittenVars; 
 set<NameComp_t> allVars; 


 bool collectReadWriteVariables(SgStatement* stmt, 
				set<NameComp_t>& readVars, 
				set<NameComp_t>& writeVars);
 bool collectReadWriteRefs(SgStatement* stmt, 
			   std::vector<SgNode*>& readRefs, 
			   std::vector<SgNode*>& writeRefs);
  

};

#endif //_READWRITEPROP_H
