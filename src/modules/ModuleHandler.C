/*
  Didem Unat
 */

#include <rose.h>
#include "VarSym.hh"
#include "./ModuleHandler.h"
#include "../functions/FunctionHandler.h"
#include "../utils/OutputFormatting.h"

using namespace std;
using namespace SageInterface;
using namespace OutputFormat;


void ModuleHandler::processModule(SgNode* node)
{
  SgModuleStatement* modNode = isSgModuleStatement(node);
  ROSE_ASSERT(modNode);

  //get the list of modules that this modules uses 
  // use module_name ...
  ModuleHandler::getDependentModuleList(node);

  //find all the fortran subroutines in the module and process them one by one
  Rose_STL_Container<SgNode*> funcList = NodeQuery::querySubTree(node, V_SgFunctionDefinition);
  Rose_STL_Container<SgNode*>::iterator funcListItr = funcList.begin() ;
  
  for(; funcListItr != funcList.end(); funcListItr++)
    {
      SgFunctionDefinition* func= isSgFunctionDefinition(*funcListItr);
      ROSE_ASSERT(func);

      SgFunctionDeclaration* func_dec = func -> get_declaration();
      ROSE_ASSERT(func_dec);
      
      SgName func_name = func_dec->get_name();

      outputBegin("function", func_name.str()); 

      //not all the modules are declared at the beginning of a module
      //some modules specific to function can be declared at the beginning of a function 
      getDependentModuleList(isSgNode(func)); //of this function

      FunctionHandler* funcHandler = new FunctionHandler();
      funcHandler -> process(isSgNode(func));
      delete funcHandler; 

      outputEnd("function");
      
    }
}


void ModuleHandler::getDependentModuleList(SgNode* node)
{  
  const string xml_tag = "use";

  Rose_STL_Container<SgNode*> modules  = NodeQuery::querySubTree(node, V_SgUseStatement);
  Rose_STL_Container<SgNode*>::iterator m = modules.begin();

  for( ; m != modules.end(); m++)
    {
      SgUseStatement* useStmt = isSgUseStatement(*m);
      ROSE_ASSERT(useStmt);

      //use module statement appears in a module or in a function definition
      if(getEnclosingFunctionDefinition(useStmt) == NULL || isSgFunctionDefinition(node))
	{
	  SgName useName = useStmt->get_name();

	  std::vector < pair < string, string> >attr;

	  attr.push_back(make_pair("name", useName.str()));
	  outputBeginEnd(xml_tag, attr);
	}
    }
}

// eof


