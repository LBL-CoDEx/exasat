/*
  by Didem Unat 
  */

#ifndef _IFHANDLER_H
#define _IFHANDLER_H 1

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>
#include "../datatypes.h"
#include "ASTtools.hh"
#include "../handler/NodeHandler.h"
#include "../handler/ScopeStmtHandler.h"

using namespace std;
using namespace SageInterface;

class IfHandler : public ScopeStmtHandler
{
  
 private:

  IfAttributes_t ifAttr;
 
  void processIfBody(SgIfStmt* ifbody);

 public:

  bool isTrueBody; //true means true body, false means else's body

  void dump();
  void process();
  void closeXMLTag();

  IfHandler(SgNode* node, bool iforelse);
  ~IfHandler();

};

#endif //_IFHANDLER_H
