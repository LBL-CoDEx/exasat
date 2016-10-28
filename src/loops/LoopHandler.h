/*
  by Didem Unat 
  */

#ifndef _LOOPHANDLER_H
#define _LOOPHANDLER_H 1

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>
#include "../datatypes.h"
#include "ASTtools.hh"
#include "../handler/ScopeStmtHandler.h"

using namespace std;
using namespace SageInterface;

class LoopHandler : public ScopeStmtHandler 
{
  
 private:

  LoopAttributes_t loopAttr;

 public:

  void dump();
  void process();
  void closeXMLTag();

  //get this and that
  SgInitializedName* getLoopIndex()    {return loopAttr.indexVar ; }

  LoopHandler(SgNode* node);
  ~LoopHandler();

};

#endif //_LOOPHANDLER_H
