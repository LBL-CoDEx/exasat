/*
  by Didem Unat 
  */

#ifndef _FUNCTIONHANDLER_H
#define _FUNCTIONHANDLER_H 1

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>
#include "../datatypes.h"
#include "ASTtools.hh"
#include "DefUseAnalysis.h"
#include "LivenessAnalysis.h"
#include "../handler/NodeHandler.h"
#include "../handler/ScopeStmtHandler.h"

using namespace std;
using namespace SageInterface;

class FunctionHandler 
{
 public:

  void process(SgNode* node);

  FunctionHandler();
  ~FunctionHandler();
  void setFortranLanguage(bool value){isFortranLanguage = value;}

 private:

  void dump(std::vector<NodeHandler*>* selectednodes);

  vector<NodeHandler*> selected; 

  void setParentChild(stack<NodeHandler*>& handlerStack);
  void buildParentToChildrenMap(SgNode*, std::stack<NodeHandler*>& nodeStack, set<SgNode*>& visited);
  void collectLocalsNonLocals(SgNode* node);
  int isFortranLanguage;

  
};

#endif //_FUNCTIONHANDLER_H
