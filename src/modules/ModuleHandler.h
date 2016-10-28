
/*
  by Didem Unat 
  */

#ifndef _MODULEHANDLER_H
#define _MODULEHANDLER_H

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>

#include "ASTtools.hh"
#include "DefUseAnalysis.h"
#include "LivenessAnalysis.h"

using namespace std;
using namespace SageInterface;

class ModuleHandler
{
 public:

  void processModule(SgNode* node);
  static void getDependentModuleList(SgNode* node);

  ModuleHandler(){}

 private:
  
};


#endif //_MODULEHANDLER_H
