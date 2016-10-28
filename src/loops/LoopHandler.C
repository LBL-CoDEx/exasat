
#include <rose.h>
#include "LoopHandler.h"

#include "../scalars/Scalars.h"
#include "../utils/ArrayUtils.h"
#include "../utils/OutputFormatting.h"


#include "ArrayAnnot.h"
#include "ArrayInterface.h"
#include "LoopTransformInterface.h"
#include "AstInterface_ROSE.h"
#include "DepInfoAnal.h"

#define MAX_NUM_LOOPS 100

using namespace std;
using namespace SageInterface;
using namespace ArrayUtils;
using namespace OutputFormat; 

LoopHandler::LoopHandler(SgNode* node)
{
  myNode = node;
  visited = false;
  parent = NULL;

}

LoopHandler::~LoopHandler()
{
  //clean up 
  accessInfo.clear();
  readOnlyVars.clear();
  writeOnlyVars.clear();
  bothReadWrittenVars.clear();
  scalarVars.clear();
  firstAccess.clear();
  children.clear();
}

void LoopHandler::process()
{
  //Step 1: Get the loop attributes such as lower bound, upper bound etc 
  //Step 2: Process the loop body 
  SgNode* node = this-> myNode;
  SgBasicBlock* body;

  SgFortranDo* cur_loop= isSgFortranDo(node);
  if(cur_loop) body = cur_loop->get_body(); //fortran Do loop
  else{ //C for loop
    body = isSgBasicBlock(isSgForStatement(node)->get_loop_body()); //we already standardize for statements so they always come with a basic block
  }
  ROSE_ASSERT(body);

  SgInitializedName* orig_index;
  SgExpression *orig_lower, *orig_upper, *orig_stride;

  if(cur_loop) //fortran Do loop
  {
    bool result = doLoopNormalization(cur_loop);
    bool is_canonical = isCanonicalDoLoop (cur_loop, &orig_index, & orig_lower,
                                         &orig_upper, &orig_stride, NULL, NULL, NULL);
    loopAttr.lower = orig_lower;
    loopAttr.upper = orig_upper;
  }
  else //C for loop
  { 
    bool isIncremental;
    bool isInclusiveUpperBound;
    bool is_canonical = isCanonicalForLoop (isSgForStatement(node), &orig_index, & orig_lower,
                                         &orig_upper, &orig_stride, NULL, &isIncremental, &isInclusiveUpperBound);
    if(isInclusiveUpperBound)
      loopAttr.upper = orig_upper;
    else{
      if(isIncremental)
        loopAttr.upper = buildSubtractOp(orig_upper, orig_stride);
      else
        loopAttr.upper = buildAddOp(orig_upper, orig_stride);
    }
  }
  int lineno =(isSgNode(node))->get_file_info()->get_line();

  loopAttr.lineno = lineno;
  loopAttr.indexVar = orig_index;
  loopAttr.lower = orig_lower;
  loopAttr.stride = orig_stride;

  computeFlops(isSgNode(body));

#if 1
  Rose_STL_Container<SgNode*> prefs = NodeQuery::querySubTree(isSgNode(body), V_SgPointerDerefExp);
  Rose_STL_Container<SgNode*>::iterator pref = prefs.begin();
  for(; pref != prefs.end(); pref++)
  {
    if(isSgPointerDerefExp(*pref)){
      SgVarRefExp* pointerVarRef;
      Rose_STL_Container<SgNode*> varList = NodeQuery::querySubTree(*pref, V_SgVarRefExp);
      for(Rose_STL_Container<SgNode*>::iterator it= varList.begin(); it!= varList.end(); it++){
        SgType* type = isSgVarRefExp(*it)->get_type();
        if(isSgPointerType(type)){
          pointerVarRef= isSgVarRefExp(*it);
          break;
        }
      }
      //find the enclosing loop, there is at least one (i.e. node)
      SgNode* enclosingLoop=pointerVarRef;
      while(enclosingLoop!=NULL && !isSgForStatement(enclosingLoop) && !isSgFortranDo(enclosingLoop)){
        enclosingLoop= enclosingLoop->get_parent();
      }
      ROSE_ASSERT(enclosingLoop);
//run side effect analysis
      SgFunctionDefinition* funcDef = SageInterface::getEnclosingFunctionDefinition(enclosingLoop);
      ROSE_ASSERT(funcDef != NULL);
      SgBasicBlock* funcBody = funcDef->get_body();
      ROSE_ASSERT(funcBody!= NULL);

      AstInterfaceImpl faImpl(funcBody);
      AstInterface fa(&faImpl);
      DoublyLinkedListWrap<AstNodePtr> rRef1, wRef1;
      CollectDoublyLinkedList<AstNodePtr> crRef1(rRef1),cwRef1(wRef1);
      AstNodePtr s1 = AstNodePtrImpl(enclosingLoop);
      AnalyzeStmtRefs(fa, s1, cwRef1, crRef1);

      bool isModified=false;
      Rose_STL_Container<SgNode*> varList1 = NodeQuery::querySubTree(enclosingLoop, V_SgVarRefExp);
      for(Rose_STL_Container<SgNode*>::iterator it= varList1.begin(); it!= varList1.end(); it++){
        SgVariableSymbol* sym= isSgVarRefExp(*it)->get_symbol(); 
        if(sym == pointerVarRef->get_symbol()){ //we will find at least one
          for (DoublyLinkedEntryWrap<AstNodePtr>* p = wRef1.First(); p != 0; )
          {
            DoublyLinkedEntryWrap<AstNodePtr>* p1 = p;
            p = wRef1.Next(p);
            AstNodePtr cur = p1->GetEntry();
            SgNode* sgRef = AstNodePtrImpl(cur).get_ptr();
            ROSE_ASSERT(sgRef != NULL);
            if(isSgVarRefExp(*it) == sgRef) isModified=true;
          }
        }
      }
      if(isModified==false){
        //this is a simple transformation to replace dereference exp with array reference exp, we will improve its later
        SgNode *bb= enclosingLoop;
        while(!isSgBasicBlock(bb)){
           bb= bb->get_parent();
        }
        if(isSgBasicBlock(bb)){
          SgStatement* lastStmtFound=NULL;
          SgStatementPtrList stmtList = isSgBasicBlock(bb)->get_statements();
          for (SgStatementPtrList::iterator i = stmtList.begin(); i != stmtList.end(); i++)
          {
            if(isSgExprStatement(*i))
            { 
              SgExpression* expStmt =isSgExprStatement(*i)->get_expression();
              ROSE_ASSERT(expStmt);
              if(isSgAssignOp(expStmt))
              {
                SgExpression * lhs = isSgAssignOp(expStmt)->get_lhs_operand();
                if(isSgVarRefExp(lhs)){
                  SgSymbol *symbol= isSgVarRefExp(lhs)->get_symbol();
                  if(symbol == pointerVarRef->get_symbol()){ 
                    SgExpression *rhs= isSgAssignOp(expStmt)->get_rhs_operand();
                    if(isSgAddressOfOp(rhs))
                      if(isSgPntrArrRefExp(isSgAddressOfOp(rhs)->get_operand())){
                         SgExpression * op= isSgPointerDerefExp(*pref)->get_operand();
                         if(isSgAddOp(op) || isSgSubtractOp(op)){
                           if(isSgBinaryOp(op)->get_lhs_operand() == pointerVarRef || isSgBinaryOp(op)->get_rhs_operand() == pointerVarRef){
                              SgExpression* arrNameExp = NULL;
                              std::vector<SgExpression*> subscripts; 
                              SgExpression* replacementArrayRefExp= copyExpression(isSgPntrArrRefExp(isSgAddressOfOp(rhs)->get_operand()));
                              isArrayReference(replacementArrayRefExp, &arrNameExp, &subscripts);
                              SgExpression* copiedExp= copyExpression(*(subscripts.end()-1));
                              replaceExpression(isSgBinaryOp(op)->get_lhs_operand() == pointerVarRef?isSgBinaryOp(op)->get_lhs_operand():isSgBinaryOp(op)->get_rhs_operand(), copiedExp); 
                              replaceExpression(*(subscripts.end()-1), copyExpression(op));             
                              replaceExpression(isSgPointerDerefExp(*pref), replacementArrayRefExp); 
                           }
                         }
                      }
                  }
                }
              }
            }
            if(isSgVariableDeclaration(*i))
            {
              SgVariableDeclaration *declStmt = isSgVariableDeclaration(*i);
              SgInitializedName* initializedName = *(declStmt->get_variables().begin());
              if(initializedName){
                SgSymbol *symbol= initializedName->get_symbol_from_symbol_table();
                if(symbol == pointerVarRef->get_symbol()){
                  SgInitializer* initPtr = initializedName->get_initptr();
                  if(initPtr){
                    SgExpression* exp = ((SgAssignInitializer*)initPtr)->get_operand_i();
                    if(isSgAddressOfOp(exp))
                      if(isSgPntrArrRefExp(isSgAddressOfOp(exp)->get_operand())){
                         SgExpression * op= isSgPointerDerefExp(*pref)->get_operand();
                         if(isSgAddOp(op) || isSgSubtractOp(op)){
                           if(isSgBinaryOp(op)->get_lhs_operand() == pointerVarRef || isSgBinaryOp(op)->get_rhs_operand() == pointerVarRef){
                              SgExpression* arrNameExp = NULL;
                              std::vector<SgExpression*> subscripts;
                              SgExpression* replacementArrayRefExp= copyExpression(isSgPntrArrRefExp(isSgAddressOfOp(exp)->get_operand()));
                              isArrayReference(replacementArrayRefExp, &arrNameExp, &subscripts);
                              SgExpression* copiedExp= copyExpression(*(subscripts.end()-1));
                              SgExpression* pointerExp= copyExpression(op);
                              replaceExpression(isSgBinaryOp(op)->get_lhs_operand() == pointerVarRef?isSgBinaryOp(pointerExp)->get_lhs_operand():isSgBinaryOp(pointerExp)->get_rhs_operand(), copiedExp);
                              replaceExpression(*(subscripts.end()-1), pointerExp);
                              replaceExpression(isSgPointerDerefExp(*pref), replacementArrayRefExp);
                           }
                         }
                      }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
#endif
  //scalar variables in the loop
  //we include the scalars used in the loop increment and condition to the current loop
  if(cur_loop)//fortran Do loop
    Scalars::getActiveScalarVars(isSgNode(cur_loop), scalarVars);
  else //C for loop
    Scalars::getActiveScalarVars(isSgNode(isSgForStatement(node)), scalarVars);

  computeArrayAccesses(isSgNode(body));
}

void LoopHandler::closeXMLTag()
{
  outputEnd("loop");
}

void LoopHandler::dump()
{
  std::vector< pair<string, string> > attributes;
  
  ostringstream linenum_str ;
  linenum_str << loopAttr.lineno;

  SgInitializedName* index = loopAttr.indexVar;
  string index_str = index->get_name().str(); 

  //loop attributes 
  attributes.push_back(make_pair("linenum", linenum_str.str()));
  attributes.push_back(make_pair("loopvar", index_str));
  attributes.push_back(make_pair("lowerbound",  loopAttr.lower->unparseToString()));
  attributes.push_back(make_pair("upperbound", loopAttr.upper->unparseToString()));
  attributes.push_back(make_pair("stride", loopAttr.stride->unparseToString()));

  ostringstream count_add_str, count_mul_str, count_special_str, count_div_str; 
  //flops 
  count_add_str << flops.addOps;
  count_mul_str << flops.mulOps;
  count_div_str << flops.divOps;
  count_special_str << flops.specialOps;
  
  attributes.push_back(make_pair("adds", count_add_str.str()));
  attributes.push_back(make_pair("multiplies", count_mul_str.str() ));
  attributes.push_back(make_pair("divides", count_div_str.str() ));
  attributes.push_back(make_pair("specials", count_special_str.str()));

  //Output loop attributes 
  outputBegin("loop", attributes);
  
  //we need to remove the dublicates first for scalars because those are varref
  //when we convert them into iname, they are unique
  ScalarRWMap_t scalarRW;
  FirstAccessMap_t scalarFirstAccess;
  set<SgInitializedName*> scalarInames; 
  
  Scalars::getScalarReadWriteCounts(scalarRW, scalarFirstAccess, scalarVars);
  Scalars::getScalarInames(scalarInames, scalarVars);

  //dump scalar variables
  dumpScalarInfo(scalarInames, scalarRW, scalarFirstAccess, "scalar");
  
  dumpVectorInfo(readOnlyVars, "readonly");
  dumpVectorInfo(writeOnlyVars, "writeonly");
  dumpVectorInfo(bothReadWrittenVars, "readwrite");
  
}
//eof








