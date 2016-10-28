#ifndef _SPECIALBB_H
#define _SPECIALBB_H 

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

class SpecialBBHandler : public ScopeStmtHandler
{
  
 private:

  SpecialBBAttributes_t spBBAttr;
  SpecialBBAttributes_t specialBBAttr;
 
 public:

  void dump();
  void process();
  void closeXMLTag();

  SpecialBBHandler(SgNode* node);
  ~SpecialBBHandler();

};

#endif 
