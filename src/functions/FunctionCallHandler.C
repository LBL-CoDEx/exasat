/*
  Didem Unat
 */

#include <rose.h>

#include "VarSym.hh"
#include "FunctionCallHandler.h"

#include "../constants.h"
#include "../datatypes.h"
#include "../utils/OutputFormatting.h"

using namespace std;
using namespace SageInterface;
using namespace OutputFormat;

FunctionCallHandler::FunctionCallHandler(SgNode* node)
{
  myNode = node; 
  visited = false;
}

FunctionCallHandler::~FunctionCallHandler()
{
  (fCallInfo.arg_map).clear();
  children.clear();
}

bool FunctionCallHandler::isFortranMathFunction(SgNode* node)
{

  Rose_STL_Container<SgNode*> funcCalls = NodeQuery::querySubTree(node, V_SgFunctionRefExp); 
  Rose_STL_Container<SgNode*>::iterator funcItr = funcCalls.begin();
  
  for( ; funcItr != funcCalls.end(); funcItr++)
    {
      SgFunctionRefExp* exp = isSgFunctionRefExp(*funcItr);
      SgFunctionSymbol* sym = exp->get_symbol();
      ROSE_ASSERT(sym);
      
      string func_name = sym->get_name().str();

      if(math_func_list.find(func_name) != math_func_list.end())
	  return true;

      else if(math_func_low_lat_list.find(func_name)!= math_func_low_lat_list.end())
	return true;

      else if(fortran_intrinsic_func_list.find(func_name) != fortran_intrinsic_func_list.end())
	return true;
    }
  return false;
}

void FunctionCallHandler::process()
{
  //Step 1: Get the function refences and its symbol, set the default values
  //Step 2: Get the declaration and function scope
  //Step 3: There might be more than function definitions with the same interface, create a candidate list
  //Step 4: Parse the function call expressions and its parameters
  //Step 5: Match the function arguments with the call arguments. 
  //Step 6: They might be optional arguments with explicit naming (ActualArgumentExpression),
  //match the argument name with the name on the call site
  //Step 7: If all the arguments match with the parameters, that's the right function, return. Otherwise 
  //check the next candidate function symbol

  /* Didem Unat, Jan 17, 2014
     On the function call arguments 

     There might be multiple interfaces for a function
     For example: 
        interface simple_func
	  procedure simple_func_1
	  procedure simple_func_2
	end interface
    
     ROSE gives the references to the first one on the list in the interface, 
     Thus argument list might be wrong, 
     Simple_func_1 might have no params and simple_func_2 might have 1 param 
     Even the call is simple_func(A), ROSE returns the pointer to the simple_func_1 

     So, we get the function symbol, its declaration and then scope and then the symbol table 
     and then go through all the function symbols that have the same name (string).
     SgRenameSymbol has the interface name (e.g. simple_func) 
     but original function symbol has the actual name for example simple_func_1.
 
     Then We go through the candidateFunc symbols and see if any of them has a matching 
     parameter list with the func call arguments. 

   */

  SgNode* node = this-> myNode; 

  //Step 1
  SgFunctionCallExp* func_call= isSgFunctionCallExp(node);
  ROSE_ASSERT(func_call);

  //set the default values 
  fCallInfo.lineno = func_call -> get_file_info() -> get_line();
  fCallInfo.func_name = "";
  fCallInfo.func_rename = "";
  fCallInfo.module_name = "unknown";

  SgFunctionRefExp * func_ref = isSgFunctionRefExp(func_call->get_function());

  if (func_ref != NULL) {

    std::set<SgFunctionSymbol*> candidateFuncs;

    //Step 1: Get the function declaration and its arguments 
    SgFunctionSymbol* func_call_sym = func_ref->get_symbol_i();
    ROSE_ASSERT(func_call_sym);
    string func_call_name = func_call_sym->get_name().str();
    fCallInfo.func_name = func_call_name; 

    candidateFuncs.insert(func_call_sym);

    //Step 2: get declaration and scope
    SgFunctionDeclaration* func_decl = func_call_sym->get_declaration();
    ROSE_ASSERT(func_decl);

    SgScopeStatement* func_scope= func_decl-> get_scope();
    ROSE_ASSERT(func_scope);

    //Step 3: Get the candidate function symbols 
    SgSymbolTable* sym_table = func_scope->get_symbol_table();
    SgNodeSet nodeset = sym_table->get_symbolSet();

    for (SgNodeSet::iterator sym=nodeset.begin(); sym!=nodeset.end();sym++)
      {
	SgNode* symbol = *sym;  

	if(isSgRenameSymbol(symbol))//if function has an interface with multiple names
	  {
	    SgRenameSymbol* renamesymbol = isSgRenameSymbol(symbol);

	    //check if the function name matches with the call name 
	    if(renamesymbol->get_name().str() == func_call_name){
		candidateFuncs.insert(isSgFunctionSymbol(renamesymbol));
	    }
	  }
      }
    
    //Step 4, 5, 6, and 7 go through the candidate func symbols and check any of them 
    //has a matching parameter list with the func call arguments
    for(std::set<SgFunctionSymbol*>::iterator sym = candidateFuncs.begin(); 
	sym!= candidateFuncs.end(); sym++)
      {
	SgFunctionSymbol* func_sym = *sym; 
	ROSE_ASSERT(func_sym);

	//interface symbol 
	SgRenameSymbol* rename_sym = isSgRenameSymbol(func_sym);

	if(rename_sym){
	  //interface name 
	  fCallInfo.func_rename = rename_sym->get_name().str();

	  //original function symbol 
	  func_sym = isSgFunctionSymbol(rename_sym -> get_original_symbol());
	  ROSE_ASSERT(func_sym);
	}

	//orgininal name
	fCallInfo.func_name = func_sym->get_name().str();

	SgFunctionDeclaration* func_decl = func_sym->get_declaration();
	ROSE_ASSERT(func_decl);

	//Get the function parameters
	SgInitializedNamePtrList& funcParamList = func_decl->get_parameterList()->get_args();
	
	//Get the arguments in a function call exp 
	SgExpressionPtrList& callArgList = func_call->get_args()->get_expressions();
	
	if(matchArgumentsWithParameters(callArgList, funcParamList))
	  {
	    //Getting module name that the function is defined
	    SgScopeStatement* func_scope= func_decl-> get_scope();
	    ROSE_ASSERT(func_scope);
	    
	    if(isSgClassDefinition(func_scope))
	      {
		SgClassDeclaration* class_decl = isSgClassDefinition(func_scope)-> get_declaration();
		if(isSgModuleStatement(class_decl))
		  {
		    SgModuleStatement* module_stmt = isSgModuleStatement(class_decl);
		    SgName module_name = module_stmt->get_name();
		    fCallInfo.module_name = module_name.str();
		  }
	      }
	    break;
	  }//end of match 
      }
  }
}


bool FunctionCallHandler::matchArgumentsWithParameters(SgExpressionPtrList callArgList, 
						       SgInitializedNamePtrList funcParamList)
{
  if(callArgList.size() > funcParamList.size())
    return false; 

  //can we exhaust all the non-optional parameters?
  //get the number of non-optional parameters
  //if that is lower than the number of arguments
  int nonoptional = 0; 
  SgInitializedNamePtrList::iterator nameItr;
  for(nameItr = funcParamList.begin(); nameItr!= funcParamList.end(); nameItr++)
    {
      SgInitializedName* iname = *nameItr; 
      ROSE_ASSERT(iname);
      
      SgDeclarationStatement* var_decl = iname ->get_declaration();;
      ROSE_ASSERT(var_decl);
      
      //Fortran variables can have intent modifier, we skip those because they are not local
      if(!(var_decl -> get_declarationModifier().get_typeModifier().isOptional()))
	nonoptional++;
    }
  //we don't have enough arguments, return false
  if(callArgList.size() < nonoptional)
    return false; 

  //Match the function arguments with the call arguments 
  SgExpressionPtrList::iterator callItr= callArgList.begin();
  SgExpressionPtrList::iterator actualArgItr = callArgList.end();
  SgInitializedNamePtrList::const_iterator paramItr = funcParamList.begin();
  
  for(callItr = callArgList.begin(); callItr != callArgList.end(); callItr++, paramItr++)
    {
      SgExpression* argExp = *callItr; 
      ROSE_ASSERT(argExp);
      
      if(isSgActualArgumentExpression(argExp)){
	
	//this is an optional parameter and explicit naming is used
	//for example: call func(TIME=2)
	//once we see one argument of this type, 
	//the rest of the arguments are in this type, break
	actualArgItr = callItr ; 
	break; 
      }
      else 
	{
	  SgInitializedName* iname_param  = *paramItr;
	  ROSE_ASSERT(iname_param);

	  string arg_str = argExp-> unparseToString(); 
	  
	  fCallInfo.arg_map.insert(make_pair(iname_param, arg_str ));
	}
    }

  //Handle actual argument expressions for ex. call func(1, TIME=2)
  //continue where the iterator was at
  for(; actualArgItr != callArgList.end(); actualArgItr++)
    {
      SgExpression* paramCallExp = *actualArgItr; 
      ROSE_ASSERT(paramCallExp);
      
      SgActualArgumentExpression *actualArgExp = isSgActualArgumentExpression(paramCallExp);
      if(actualArgExp != NULL)
	{
	  string param_str =  actualArgExp->get_argument_name().str() ;
	  string arg_str = actualArgExp->get_expression() -> unparseToString();
	  
	  SgInitializedName* iname_param = NULL;
	  
	  //search the initializedname in the function parameter list
	  SgInitializedNamePtrList::iterator nameItr ;
	  for(nameItr = funcParamList.begin(); nameItr!= funcParamList.end(); nameItr++)
	    {
	      SgInitializedName* iname = *nameItr; 
	      ROSE_ASSERT(iname);
	      
	      string param_name_str = iname-> get_name().str();
	      if(param_str == param_name_str){ //found the initializedname 
		iname_param = iname; 
		break;
	      }
	    }
	  
	  ROSE_ASSERT(iname_param);
	  
	  fCallInfo.arg_map.insert(make_pair(iname_param, arg_str));
	  
	}
    }
  //everything matched, return true
  return true;
}

inline
string find_replace(string &originalStr, string replacedSubStr, string replacementStr){
    size_t position= originalStr.find(replacedSubStr);
    while(position!= std::string::npos){
      originalStr.replace(position, replacedSubStr.length(), replacementStr);
      position= originalStr.find(replacedSubStr, position + replacementStr.size());
    }
    return originalStr;
}

void FunctionCallHandler::dump()
{
  string func_name = fCallInfo.func_name; 
  string module_name = fCallInfo.module_name;
  string func_rename = fCallInfo.func_rename;

  std::map<SgInitializedName*, string>* arg_map = &(fCallInfo.arg_map);

  ostringstream linenum_str ;
  linenum_str << fCallInfo.lineno ; 

  std::vector< pair<string, string> > attributes;
  attributes.push_back(make_pair("linenum", linenum_str.str()));
  attributes.push_back(make_pair("name", func_rename));
  if(func_name != "")
    attributes.push_back(make_pair("origname", func_name));
  attributes.push_back(make_pair("module", module_name));
  attributes.push_back(make_pair("flops", "unknown"));

  outputBegin("funccall", attributes);

  //going through argument list 
  std::map<SgInitializedName*, string>::iterator argItr = arg_map-> begin();
  for(; argItr != arg_map-> end(); argItr++)
    {
      SgInitializedName* iname = argItr -> first;
      string arg_name = argItr-> second;

      string param_name = iname-> get_name().str(); 

      find_replace(arg_name, "&", "&amp;");
      find_replace(arg_name, "<", "&lt;");
      find_replace(arg_name, ">", "&gt;");
      find_replace(arg_name, "\"", "&quot;");
      find_replace(arg_name, "'", "&apos;");

      std::vector< pair<string, string> > arg_attr;
      arg_attr.push_back(make_pair("paramname", param_name));
      arg_attr.push_back(make_pair("argname", arg_name));

      outputBeginEnd("arg", arg_attr);
    }

}

void FunctionCallHandler::closeXMLTag()
{
  outputEnd("funccall");
}
// eof



