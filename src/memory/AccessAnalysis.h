/*
  by Didem Unat 
  */

#ifndef _ACCESSANALYSIS_H
#define _ACCESSANALYSIS_H


#include <iostream>
#include <string>
#include <sstream>
#include "../datatypes.h"

using namespace std;
typedef struct
{
  bool loopVar; 
  bool sign;//true is positive
  SgExpression* lhs;
  SgExpression* rhs;
} clause;

class AccessAnalysis
{
 private:

  int loopNo;
  std::vector<SgInitializedName*>* loopIndices;
  vector< std::set<SgVarRefExp*>*>* scalarVars;
//  vector< std::list<SgVarRefExp*>*>* sortedVars;
  list<clause*> clauseList;
  list<SgNode*> flatArrayRefList;
  list<SgNode*> hArrayRefList;

  bool isPartofFuncCallExp(SgNode* node);

  SgInitializedName* getNonNullLoopIndex(int pos);

  void computeOffsets(AccessPattern_t& this_access,
		      vector<SgExpression*> subscripts);
  void factorizeTerms(AccessPattern_t& this_access, SgExpression*);
  void normalizeTerms(SgNode*, bool);

  string getAccessString(vector<SgExpression*> subs);

  string getComponent(vector<SgExpression*> subs);

  void computeDependentLoopVar(AccessPattern_t& this_access,
			       vector<SgExpression*> &subscripts);
  void eraseClauseList(){
    while(!clauseList.empty()){
      clause *c= clauseList.front();
      delete c;
      clauseList.pop_front();
    }
  }

 public:

  AccessAnalysis(){}

  void setLoopNo(int loopnum){ loopNo = loopnum ;}
  void setLoopIndices(std::vector<SgInitializedName*>* loopInd) {loopIndices = loopInd ;}

  void setScalarVars(std::vector< std::set<SgVarRefExp*>* >* scalarV) { scalarVars = scalarV ;}
  
  void arrayAccessAnalysis(SgNode* node, //loop body 
			   ArraysAccessInfo_t& arrList, //array-comp list mapped to access patterns 
			   FirstAccessMap_t& ); //array-comp list mapped to first access info

  bool isRead(SgExpression* exp); //checks if the expression is a read or write

  void firstWrittenOrReadAnalysis(FirstAccessMap_t& firstAccess,
				  SgExpression* exp, 
				  NameComp_t namecomp);

 
};

#endif //_ACCESSANALYSIS_H
