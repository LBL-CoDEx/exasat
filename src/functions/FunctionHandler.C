/*
  Didem Unat
 */

#include <rose.h>

#include "VarSym.hh"
#include "FunctionHandler.h"

#include "../flops/Flops.h"

#include "../handler/ScopeStmtHandler.h"
#include "../loops/LoopHandler.h"
#include "../conditionals/IfHandler.h"
#include "./FunctionCallHandler.h"
#include "../specialBasicBlocks/specialBBHandler.h"

#include "../wholearrays/WholeArrayExp.h"
#include "../utils/OutputFormatting.h"
#include "../utils/ArrayUtils.h"

using namespace std;
using namespace SageInterface;
using namespace OutputFormat;
using namespace ArrayUtils;


FunctionHandler::FunctionHandler()
{
  isFortranLanguage = true;
}

FunctionHandler::~FunctionHandler()
{
  vector<NodeHandler*>::iterator node;
  for(node = selected.begin(); node!= selected.end(); node++)
    {
      NodeHandler* handler = *node; 
      delete handler; // calls the destructor of the respective class
    }
  selected.clear();
}

void findLoopIndices1(SgBinaryOp* op, std::vector< std::pair<SgInitializedName*, SgExpression*> > &indexVars){
  SgExpression *lhs = op->get_lhs_operand();
  SgExpression *rhs = op->get_rhs_operand();
  if(isSgVarRefExp(lhs)){
    indexVars.push_back(make_pair(isSgVarRefExp(lhs)->get_symbol()->get_declaration(), rhs));
  }
  if(isSgBinaryOp(lhs))findLoopIndices1(isSgBinaryOp(lhs), indexVars);
  if(isSgBinaryOp(rhs))findLoopIndices1(isSgBinaryOp(rhs), indexVars);
}
SgExpression* findStride1(SgForStatement* forNode, SgInitializedName* indexVar){
  SgExpression* inc= forNode->get_increment();
  if(inc){
    //look for increment exp associated with this indexVar, e.g. i++
    Rose_STL_Container<SgNode*> unaryOps  = NodeQuery::querySubTree(isSgNode(inc), V_SgUnaryOp);
    for(Rose_STL_Container<SgNode*>::iterator it= unaryOps.begin(); it!= unaryOps.end(); it++){
      SgUnaryOp* op = isSgUnaryOp(*it);
      if(op){
        SgVarRefExp* varRef= isSgVarRefExp(op->get_operand());
        if(varRef->get_symbol()->get_declaration() == indexVar){
          if(isSgMinusMinusOp(op)) return buildIntVal(-1);
          if(isSgPlusPlusOp(op)) return buildIntVal(1);
        }
      }
    }
    //if we didn't find the index in a unary op form, we look for binary ops (e.g. i+=1)
    Rose_STL_Container<SgNode*> binaryOps  = NodeQuery::querySubTree(isSgNode(inc), V_SgBinaryOp);
    for(Rose_STL_Container<SgNode*>::iterator it= binaryOps.begin(); it!= binaryOps.end(); it++){
      SgBinaryOp* op = isSgBinaryOp(*it);
      if(op){
        SgExpression* lhs= op->get_lhs_operand();
        SgExpression* rhs= op->get_rhs_operand();
        SgVarRefExp* varRef= isSgVarRefExp(lhs);
        if(varRef)
          if(varRef->get_symbol()->get_declaration() == indexVar){
            if(isSgPlusAssignOp(op)) return rhs;
            if(isSgMinusAssignOp(op)) return buildMinusOp(rhs);
          }
      }
    }
    //if we didn't find the index in a binary op form, we look for assign ops (e.g. i=i+1)
    Rose_STL_Container<SgNode*> assignOps  = NodeQuery::querySubTree(isSgNode(inc), V_SgAssignOp);
    for(Rose_STL_Container<SgNode*>::iterator it= assignOps.begin(); it!= assignOps.end(); it++){
      SgBinaryOp* op = isSgAssignOp(*it);
      if(op){
        SgExpression* lhs= op->get_lhs_operand();
        SgExpression* rhs= op->get_rhs_operand();
        SgVarRefExp* varRef= isSgVarRefExp(lhs);
        if(varRef)
          if(varRef->get_symbol()->get_declaration() == indexVar){
              if(isSgBinaryOp(rhs)){
                if(isSgVarRefExp(isSgBinaryOp(rhs)->get_lhs_operand()))
                  if(isSgVarRefExp(isSgBinaryOp(rhs)->get_lhs_operand())->get_symbol()->get_declaration() == indexVar){
                    if(isSgAddOp(rhs)) return isSgBinaryOp(rhs)->get_rhs_operand();
                    if(isSgSubtractOp(rhs)) return buildMinusOp(isSgBinaryOp(rhs)->get_rhs_operand());
                  }
              }
          }
      }
    }
    //We couldn't find/detect the stride
    return buildNullExpression();
  }
}

void removeRestrictKeywords(SgNode* node){
  Rose_STL_Container<SgNode*> refs = NodeQuery::querySubTree(node, V_SgVarRefExp);
  Rose_STL_Container<SgNode*>::iterator it = refs.begin();
  for(; it != refs.end(); it++){
    SgVarRefExp *ref= isSgVarRefExp(*it);
    SgType* type= ref->get_type();
    if(type){
      if(type->unparseToString().size()>12)
      if(type->unparseToString().substr(type->unparseToString().size()-12, 12) == "__restrict__"){
        SgModifierNodes* modifier=  type->get_modifiers();
        if(modifier){
          SgTypeModifier *type_modifier= isSgTypeModifier(modifier);
          if(type_modifier->isRestrict()){ 
            type_modifier->unsetRestrict();
          }
        }else{
          SgInitializedName* initName= ref->get_symbol()->get_declaration();
          if(initName){
            string typeStr= type->unparseToString();
            string removeRestrictStr= typeStr.erase(type->unparseToString().size()-15, type->unparseToString().size()-1); //remove " * __restrict__"
            SgType* newBaseType= buildOpaqueType(typeStr.c_str(), initName->get_scope());
            SgPointerType* newType= buildPointerType(newBaseType);
            newType->set_base_type(newBaseType);
            initName->set_type(newType);
          }
        }
      }
    }
  }
} 

void FunctionHandler::process(SgNode* node)
{
  //Step 1: convert whole array expressions to fortran do loops
  //Step 2: get all the variable information for the function (local vs non local vars)
  //Step 3: Create the tree structure of selected nodes
  //Step 4: process all the do loops, if stmts and func call exp
  //Step 5: Dump information collected in prev step

  SgFunctionDefinition* func_def= isSgFunctionDefinition(node);
  ROSE_ASSERT(func_def);

  removeRestrictKeywords(node);

  SgStatement* func_body = func_def -> get_body();
  ROSE_ASSERT(func_body);

  //In order to treat whole array statements as same as fortran do loops, we convert all 
  //the whole array references into nested loop
    WholeArrayExp* wholeArr = new WholeArrayExp();
    wholeArr->convertWholeArraysToFortranDos(isSgNode(func_body));
    delete wholeArr; 

  //creates a tree like data structure and construct their objects for dumping the results 
  //selected objects are kept in selected list and the tree structure is kept in childList
  std::stack<NodeHandler*> nodeHandlerStack;
  std::set<SgNode*> visited;
  buildParentToChildrenMap(node, nodeHandlerStack, visited);

  //preprocess loops that have 2 or more loop variables
  vector<SgNode*> forLoops= NodeQuery::querySubTree(node, V_SgForStatement); 
  for(vector<SgNode*>::iterator it= forLoops.begin(); it!=forLoops.end(); it++){
    SgForStatement* forLoop= isSgForStatement(*it); 
    vector< std::pair<SgInitializedName*, SgExpression*> > indexVars;
    std::vector<SgStatement*> forInitStmtList= forLoop->get_init_stmt();
    for(int i=0 ; i<forInitStmtList.size(); i++){
      SgStatement* initStmtList= isSgExprStatement(forInitStmtList[i]);
      if(initStmtList){
        SgExpression* exp= isSgExprStatement(initStmtList)->get_expression();
        if(exp){
          if(isSgBinaryOp(exp)) findLoopIndices1(isSgBinaryOp(exp), indexVars);
        }
      }
    }
    if(indexVars.size()>1){
      SgExpression* testExp= forLoop->get_test_expr();
      if(isSgLessThanOp(testExp) || isSgLessOrEqualOp(testExp) || isSgGreaterThanOp(testExp) || isSgGreaterOrEqualOp(testExp) || isSgEqualityOp(testExp)|| isSgNotEqualOp(testExp))
      {  
        if(isSgVarRefExp(isSgBinaryOp(forLoop->get_test_expr())->get_lhs_operand()))
        {
          SgInitializedName *indexVar = isSgVarRefExp(isSgBinaryOp(forLoop->get_test_expr())->get_lhs_operand())->get_symbol()->get_declaration();
          SgExpression* stride0= findStride1(forLoop, indexVar);
          if(!isSgNullExpression(stride0)){
            std::vector< std::pair<SgInitializedName*, SgExpression*> >::iterator itr;
            for(itr= indexVars.begin(); itr != indexVars.end(); itr++)
              if(itr->first==indexVar) break;
            SgExpression *init0=itr->second;
            for(itr= indexVars.begin(); itr != indexVars.end(); itr++){
              SgExpression *init, *stride;
              if(itr->first != indexVar){
                stride = findStride1(forLoop, itr->first);
                init= itr->second;
                if(isSgIntVal(init0) && isSgIntVal(stride0) && isSgIntVal(init) && isSgIntVal(stride)){
                  SgExpression* replacementExp;
                  if(isSgIntVal(stride0)->get_value()==1){
                    int offset= isSgIntVal(init)->get_value()- isSgIntVal(init0)->get_value();
                    int speedDiff= isSgIntVal(stride)->get_value()- isSgIntVal(stride0)->get_value();
                    if(speedDiff==0)replacementExp= isSgBinaryOp(testExp)->get_lhs_operand();
                    else replacementExp= buildMultiplyOp(buildIntVal(1+speedDiff), isSgBinaryOp(testExp)->get_lhs_operand());
                    if(offset!=0) replacementExp= buildAddOp(replacementExp, buildIntVal(offset));
                  }else{
                      ;
                  }
                  Rose_STL_Container<SgNode*> refs = NodeQuery::querySubTree(forLoop->get_loop_body(), V_SgVarRefExp);
                  Rose_STL_Container<SgNode*>::iterator ref = refs.begin();
                  for(; ref != refs.end(); ref++){
                    if(isSgVarRefExp(*ref))
                    if(isSgVarRefExp(*ref)->get_symbol())
                    if(isSgVarRefExp(*ref)->get_symbol()->get_declaration()== itr->first)
                      replaceExpression(isSgVarRefExp(*ref), replacementExp);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //process selected nodes
  //selected nodes can be typed as function body, ifstmt, doloop and funccall exp
  //we could have processed them in the previous function but it is cleaner
  vector<NodeHandler*>::iterator nd;
  for(nd = selected.begin(); nd != selected.end(); nd++)
    {
      NodeHandler* handler = *nd ; 
      handler -> process();
    }

  //recursively dump the results by respecting tree structure 
  dump(&selected);
  //get the function level information about local and non local variables 
  collectLocalsNonLocals(isSgNode(func_body));
  outputEnd("function");

}


void FunctionHandler::buildParentToChildrenMap(SgNode* node, 
					       std::stack<NodeHandler*>& handlerStack, 
					       std::set<SgNode*>& visited)
{
  //This function creates a nested node structure and 
  //Create a respective object of its kind
  //Currently nodes can be SgIfStmt, SgFortranDo and SgFunctionCallExp
  //SfIfStmt or SgFortranDo may have nesting so we recursively call buildParentToChildrenMap from 
  //these nodes 

  //    Selected list preserves the order of XML output 
  //    Children List keeps the structure of child list for each parent 

  /*in a nested structure like this one
    L1
     L2 
       L3
       L4
     L5
       L6 
    L7 

   Selected list will be
    L1 L2 L3 L4 L5 L6 L7 
 
   Children will be 
    L1 -> L2, L5
    L2 -> L3, L4
    L5 -> L6
    the rest has an empty list

   Parent List
    L1 -> NULL
    L2 -> L1 
    L3 -> L2 etc 
   */
  if(!node) return;

  //Get all the loops in the nest including itself
  Rose_STL_Container<SgNode*> nodeList = NodeQuery::querySubTree(node, V_SgLocatedNode);

  if(nodeList.size() == 0) //leaf node, no children
    {
      return; 
    }    
  else
    {
      Rose_STL_Container<SgNode*>::iterator nodeItr = nodeList.begin();
      for(; nodeItr!= nodeList.end(); nodeItr++)
	{
	  SgNode* node = *nodeItr; 

	  //We are interested in the following node types 
	  if(!(isSgIfStmt(node)||
	       isSgFortranDo(node)||isSgForStatement(node)||isSgFunctionCallExp(node) || isSgBasicBlock(node)))
	    continue;

          if(isSgBasicBlock(node))
	    if(!isSgFunctionDefinition(node->get_parent())) continue;

	  //Node has to be unvisited
	  if(visited.find(node) == visited.end())
	    {
	      visited.insert(node);

	      if(isSgIfStmt(node)) //there true and false body
		{		  
		  SgIfStmt* ifstmt = isSgIfStmt(node);
		  SgBasicBlock* true_body = isSgBasicBlock(ifstmt -> get_true_body());
		  SgBasicBlock* false_body = isSgBasicBlock(ifstmt -> get_false_body());		  

		  IfHandler* ifTrueHandler = new IfHandler(node, true);
		  
		  handlerStack.push(ifTrueHandler);
		  selected.push_back(ifTrueHandler);
		  buildParentToChildrenMap(true_body, handlerStack, visited);
		  setParentChild(handlerStack);

		  if(false_body){ 
		    //this is not enough, false_body is never NULL
		    //ROSE creates a false body for all ifs, check if there are any statements in the body 
		    if(false_body -> get_statements().size() > 0){
		      
		      IfHandler* elseHandler = new IfHandler(node, false);

		      handlerStack.push(elseHandler);
		      selected.push_back(elseHandler);
		      buildParentToChildrenMap(false_body, handlerStack, visited);
		      setParentChild(handlerStack);
		    }//end of false body
		  }
		}
	      else if (isSgFortranDo(node))
		{		  
		  LoopHandler* loopHandler = new LoopHandler(node);
		  SgFortranDo* forstmt = isSgFortranDo(node);
		  SgBasicBlock* body = forstmt -> get_body();
		  
		  handlerStack.push(loopHandler);
		  selected.push_back(loopHandler);
		  buildParentToChildrenMap(body, handlerStack, visited);		  
		  setParentChild(handlerStack);

		}
              else if (isSgForStatement(node))
                {
                  LoopHandler* loopHandler = new LoopHandler(node);
                  SgForStatement* forstmt = isSgForStatement(node);
                  SgStatement* body_stmt = forstmt -> get_loop_body();
                  SgBasicBlock* body =NULL;
                  if(!isSgBasicBlock(body_stmt)){
		    body = buildBasicBlock(body_stmt);
                    body->set_parent(forstmt);
                    forstmt->set_loop_body(body);
                  }else{
		     body = isSgBasicBlock(body_stmt);
                  }
                  handlerStack.push(loopHandler);
                  selected.push_back(loopHandler);
                  buildParentToChildrenMap(body, handlerStack, visited);
                  setParentChild(handlerStack);

                }
	      else if (isSgFunctionCallExp(node))
		{
		  if(FunctionCallHandler::isFortranMathFunction(node))
		    continue;

		  FunctionCallHandler* funcCallHandler = new FunctionCallHandler(node); 
		  
		  handlerStack.push(funcCallHandler);
		  selected.push_back(funcCallHandler);
		  setParentChild(handlerStack);
		}
	      else if(isSgBasicBlock(node))
		{
		  SpecialBBHandler* specialBBHandler= new SpecialBBHandler(node);
		  handlerStack.push(specialBBHandler);
		  selected.push_back(specialBBHandler);
                  buildParentToChildrenMap(isSgBasicBlock(node), handlerStack, visited);
		  setParentChild(handlerStack);
		}
	    } // end of visited
	}// end of nodes loop
    }
}

void FunctionHandler::setParentChild(stack<NodeHandler*>& handlerStack)
{
  NodeHandler* cur_handler = NULL; 
  NodeHandler* parent = NULL;
  
  if(!handlerStack.empty()){
    cur_handler = handlerStack.top();
    handlerStack.pop();
  }
  if(!handlerStack.empty())
    parent = handlerStack.top(); 
  
  //set the parent handler 
  if(cur_handler != NULL ){
    cur_handler-> setParent(parent);
    
    //add yourself to your parent's children list 
    if(parent != NULL ){
      parent -> addChild(cur_handler);
    }
  }
}

void FunctionHandler::dump(vector<NodeHandler*>* selected)
{
  //Step 1: Go through the selected NodeHandlers. NodeHandlers can be DoFortran, IfStmt, FuncCallExp etc
  //Step 2: Check if the handler is visited before
  //Step 3: If it is unvisited, then dump its info 
  //Step 4: Check if it has any children 
  //Step 5: Recursively call dump on the children
  //Step 6: Dump the closing tag for the XML key

  for(vector<NodeHandler*>::iterator nd= selected->begin(); nd!= selected->end(); nd++)
    {
      //Node Handler can be DoFortran, IfStmt, FuncCallExp etc 
      NodeHandler* nodeHandler = *nd;

      if(nodeHandler -> getVisited())
	continue; 

      else {
	//mark this handler as visited
	nodeHandler -> setVisited(true);

	//removes the children's data and makes parent's data exclusive
	nodeHandler-> removeDuplicates();

	//Dump all the info about this handler
	nodeHandler -> dump();      

	//Get the list of children's for this node
	vector<NodeHandler*>* children = nodeHandler -> getChildren();

	if(children-> size()> 0)
	  dump(children);
	
	nodeHandler -> closeXMLTag();
      }
    }
}

void FunctionHandler::collectLocalsNonLocals(SgNode* node)
{
  /**************************************************************
                  Variable Info (Local vs Non Local vars)
  ***************************************************************/
  SgStatement* func_body = isSgStatement(node);
 
  //L = {symbols defined within 'func_body'}, local variables declared within function body  
  ASTtools::VarSymSet_t L, Locals;
  ArrayUtils::collectLocalVarSyms(func_body, L);
  ArrayUtils::getOnlyArraySymbols(L, Locals); 
  OutputFormat::dump (Locals, "Locals:", "local"); 

  // U = {symbols used within function body }
  ASTtools::VarSymSet_t U, Used;
  ASTtools::collectRefdVarSyms (func_body, U);
  ArrayUtils::getOnlyArraySymbols(U, Used); 
  
  // U - L = {symbols used within 'func_body' but not defined in 'func_body'}    
  // variable references to non-local-declared variables
  ASTtools::VarSymSet_t diff_U_L;
  set_difference (Used.begin (), Used.end (), Locals.begin (), Locals.end (),
                  inserter (diff_U_L, diff_U_L.begin ()));
  OutputFormat::dump (diff_U_L, "Non-Locals:", "nonlocal");          


  /******************************************************
                    eof  Variable Info
  ********************************************************/

}


// eof


