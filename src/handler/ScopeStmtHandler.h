/*
  by Didem Unat 
  */

#ifndef _SCOPESTMTHANDLER_H
#define _SCOPESTMTHANDLER_H

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>

#include "../datatypes.h"
#include "./NodeHandler.h"

using namespace std;

class ScopeStmtHandler : public NodeHandler
{
 protected:

  //data structures
  Flops_t flops; 

  std::set<SgVarRefExp*> scalarVars;
  set<NameComp_t> readOnlyVars;
  set<NameComp_t> writeOnlyVars;
  set<NameComp_t> bothReadWrittenVars;
  ArraysAccessInfo_t accessInfo;
  FirstAccessMap_t firstAccess;

  //helper functions for dump
  void dumpVectorInfo(set<NameComp_t> vars,
		      string accesstype);

  void dumpScalarInfo( std::set<SgInitializedName*>& scalarsInames, 
		       ScalarRWMap_t& rw,
		       FirstAccessMap_t& , 
		       string xml_tag);

  void computeFlops(SgNode*);
  void computeArrayAccesses(SgNode* );

 private:

  //helper functions for remove duplicates
  void removeDuplicatesInScalars(ScopeStmtHandler* parent, ScopeStmtHandler* child);
  void removeDuplicatesInFlops(ScopeStmtHandler* parent, ScopeStmtHandler* child);
  void removeDuplicatesInArrayAccesses(ScopeStmtHandler* parent, ScopeStmtHandler* child);
  

  //helper function for dump vector 
  string fixAccessType(AccessListPerArray_t accList,
		       FirstAccess_t fa,
		       string rose_accesstype);

  void getScalarsOfParentNest(ScopeStmtHandler* cur,
			      vector< set<SgVarRefExp*>* >& scalars);

  void getLoopIndicesOfParentNest(ScopeStmtHandler* cur,
				  vector<SgInitializedName*>& loopIndices);


 public:

  void removeDuplicates();

  //get this and that
  Flops_t* getFlops()                  {return &flops; }
  ArraysAccessInfo_t* getAccessInfo()  {return &accessInfo;}
  std::set<SgVarRefExp*>* getScalars() {return &scalarVars; }

  virtual SgInitializedName* getLoopIndex()    {return NULL;}
};

#endif //_SCOPESTMTHANDLER_H
