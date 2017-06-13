//Didem Unat
//Major changes in Dec 2013

#include <string>
#include <map>
#include <boost/regex.hpp>

#include <iostream>


#include <rose.h>
#include "sage3basic.h"
#include "../datatypes.h"
#include "../functions/FunctionCallHandler.h"

#include "AccessAnalysis.h"
#include "../utils/ArrayUtils.h"
#include "../utils/Utils.h"
#include "../utils/ExaSatOptions.h"
#include "../utils/OutputFormatting.h"

#include "ArrayAnnot.h"
#include "ArrayInterface.h"
#include "LoopTransformInterface.h"
#include "AstInterface_ROSE.h"
#include "DepInfoAnal.h"

#define DEF_NUM_DIM 3

using namespace std;
using namespace SageInterface;

using namespace boost; 
using namespace ArrayUtils;
using namespace Utils; 
using namespace OutputFormat;

string AccessAnalysis::getAccessString(vector<SgExpression*> subs)
{

  //converts the subscript expression into a string for using it as a key for the access list map

  ostringstream access_str; 
  access_str << "("; 
  int numIndices = subs.size(); 

  vector<SgExpression*>::const_iterator it ;
  for(it= subs.begin(); it != subs.end(); it++)
    {
      SgExpression* index = (*it); 
      ROSE_ASSERT(index);

      access_str << index-> unparseToString() ; 

      //for proper formatting 
      if(it != subs.end() - 1)
	access_str << ",";
      else 
	access_str << ")";
    }
  string access_str_nospace = Utils::removeString(access_str.str(), " ");

  return access_str_nospace; 
}

bool AccessAnalysis::isPartofFuncCallExp(SgNode* node)
{
  SgNode* parent = node -> get_parent();
  if(parent != NULL)
    {
      if(isSgFunctionCallExp(parent) && !FunctionCallHandler::isFortranMathFunction(parent))
	{
	  return true; 
	}
      else 
	return isPartofFuncCallExp(parent);
    }
  return false;
}

bool findVarRefInIndexAndReplace(string varName, SgExpression *index, SgExpression* rhs, vector<SgExpression*> &subscripts, int idxOffset){
  SgExpression* copiedIndex= index;//deepCopy(index);
  copiedIndex->set_parent(index->get_parent());
  Rose_STL_Container<SgNode*> varList = NodeQuery::querySubTree(copiedIndex, V_SgVarRefExp);
  Rose_STL_Container<SgNode*>::iterator it = varList.begin();
  bool found=false;
  while(it!=  varList.end()){//scan over all varRefs and replace matched ones with RHS
    SgExpression* indexVar = isSgExpression(*it);
    ROSE_ASSERT(isSgVarRefExp(indexVar));
    SgInitializedName* iname = convertRefToInitializedName (indexVar);
    ROSE_ASSERT(iname);
    if(strcmp(varName.c_str(), iname->get_name().getString().c_str())==0){
      if(!found) found=true;
      SgExpression* copiedRHS = rhs;//copyExpression(rhs);
      if(indexVar!=copiedIndex)replaceExpression(indexVar, copiedRHS, true);
      else copiedIndex= copiedRHS;
    }
    it++;
  }
  if(found){
    subscripts[idxOffset]=copiedIndex; 
  }
  //else deepDelete(copiedIndex);
  return found;
}

void recursivePropagation(Rose_STL_Container<SgStatement*> stmtList, Rose_STL_Container<SgStatement*>::iterator it, SgStatement* enclosingStmt, SgExpression *index, vector<SgExpression*> &subscripts, int idxOffset){
  if(it==stmtList.end()) return;
  if(isSgStatement(*it) == enclosingStmt) return;
  recursivePropagation(stmtList, it+1, enclosingStmt, index, subscripts, idxOffset);
  SgExprStatement *stmt= isSgExprStatement(*it);
  if(stmt){
    SgExpression* expStmt =isSgExprStatement( stmt)->get_expression();
    ROSE_ASSERT(expStmt);
    if(isSgAssignOp(expStmt))
    {
       SgExpression * lhs = isSgAssignOp(expStmt)->get_lhs_operand();
       SgExpression * rhs = isSgAssignOp(expStmt)->get_rhs_operand();
       SgVarRefExp *varRef = isSgVarRefExp(lhs);
       if(varRef&&rhs){
         findVarRefInIndexAndReplace(varRef->get_symbol()->get_name().getString(), index, rhs, subscripts, idxOffset);
       }  
    }
  }else{
    SgVariableDeclaration* declStmt= isSgVariableDeclaration(*it);
    if(declStmt){
      SgInitializedName* initializedName = *(declStmt->get_variables().begin());
      if(initializedName){
        SgInitializer* initPtr = initializedName->get_initptr();
	if(initPtr){
          findVarRefInIndexAndReplace(initializedName->get_name().getString(), index, ((SgAssignInitializer*)initPtr)->get_operand_i(), subscripts, idxOffset);
        }
      }
    }
  }
  return;
}

void AccessAnalysis::arrayAccessAnalysis(SgNode* node, 
					 ArraysAccessInfo_t& arrList, 
					 FirstAccessMap_t& firstAccessMap)
{
#if 0
  Rose_STL_Container<SgNode*> prefs = NodeQuery::querySubTree(node, V_SgPointerDerefExp);
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
  //Node is a loop body or basic block 
  //Find all the array references in the block 
  //Group them into array-component pairs 
  //Each array-component pair has a list of access patterns
  //Each pattern has a number of attributes such as readcounts, write counts and is relative etc 

  Rose_STL_Container<SgNode*> refs = NodeQuery::querySubTree(node, V_SgPntrArrRefExp);
  Rose_STL_Container<SgNode*>::iterator ref = refs.begin();

  for(; ref != refs.end(); ref++)
  {
    //we group references into 2 buckets: flat and hierarchical
    //flat bucket contains array references represented in a flat fashion, e.g. A as in A[i], A and f as in A[f[i]] in C/C++ and A as in A(i,j,k) A(fi(i),fj(j),fk(k)) in Fortran
    //hierachical bucket contains array references represented in a hierarchical fashion, e.g. A[i][j][k] is represented as ((A[i])[j])[k]
    bool refOfInterest= !isSgPntrArrRefExp((*ref)->get_parent()); 
    if(!refOfInterest){ //parent is an arrayRef, still consider ref if it is not the left hand side of its parent
      refOfInterest = (isSgPntrArrRefExp((*ref)->get_parent())->get_lhs_operand() != (*ref));
    }
    if(refOfInterest){
      if(NodeQuery::querySubTree(isSgExpression(*ref), V_SgPntrArrRefExp).size()==1)//definitely considered flat
      {
        flatArrayRefList.push_back((*ref));
      }else{
        Rose_STL_Container<SgNode*> subrefs = NodeQuery::querySubTree(isSgExpression(*ref), V_SgPntrArrRefExp);
        Rose_STL_Container<SgNode*>::iterator subref = subrefs.begin();
	bool found=false;
        for(; subref != subrefs.end(); subref++){
          if(isSgPntrArrRefExp((*subref)->get_parent())){
            if(isSgPntrArrRefExp((*subref)->get_parent())->get_lhs_operand()->unparseToString().compare((*subref)->unparseToString())==0) found=true;
	  }
        }
        if(found)hArrayRefList.push_back((*ref));
        else flatArrayRefList.push_back((*ref));
      }
    }//no else
  }

  for(list<SgNode*>::iterator h_it= hArrayRefList.begin(); h_it != hArrayRefList.end(); h_it++)
  {
      SgExpression* exp = isSgExpression(*h_it);
      //ROSE_ASSERT(exp);
      SgExpression* arrNameExp = NULL;
      std::vector<SgExpression*> subscripts; 
      if(isArrayReference(exp, &arrNameExp, &subscripts)) //no else 
      {
          std::vector<SgExpression*> unifiedSubscripts; 
          SgInitializedName* array_name = convertRefToInitializedName (arrNameExp);
          ROSE_ASSERT(array_name);
          Rose_STL_Container<SgNode*> subrefs = NodeQuery::querySubTree(exp, V_SgPntrArrRefExp);
          Rose_STL_Container<SgNode*>::iterator subref = subrefs.begin();
          for(; subref != subrefs.end(); subref++){
            SgExpression* arrNameExp1 = NULL;
            std::vector<SgExpression*> subscripts1; 
            if(isArrayReference(isSgExpression(*subref), &arrNameExp1, &subscripts1)) //no else 
            {
              //SgInitializedName* array_name1 = convertRefToInitializedName (arrNameExp1);
              unifiedSubscripts.push_back(*subscripts1.begin());
            }
          }
	  //assumes the component is in the last subscript 
	  string component  = getComponent(unifiedSubscripts);

	  //create a array name - component pair 
	  NameComp_t namecomp (array_name, component);

	  //relies on the expression order that AST query returns
	  //first statements are handled first 
	  firstWrittenOrReadAnalysis(firstAccessMap, exp, namecomp);

	  //compute attributes related to this access pattern
	  AccessListPerArray_t accessList;
	  AccessPattern_t this_access;

	  string access_str = getAccessString(unifiedSubscripts);

	  //if the exp is part of a function call, we 
	  if(isPartofFuncCallExp(exp)){ 
	    this_access.readCount = 0 ; 
	    this_access.writeCount = 0 ; 
	  }
	  else {
	    this_access.readCount = isRead(exp) ? 1 : 0 ;
	    this_access.writeCount = (this_access.readCount == 0) ? 1 : 0 ;
	  }

	  ArraysAccessInfo_t::iterator arrayCompPair = arrList.find(namecomp);	  

	  //this is the first time the array-comp pair  appears in the node
	  if( arrayCompPair == arrList.end())
	    {
	      computeDependentLoopVar(this_access, unifiedSubscripts);
	      computeOffsets(this_access, unifiedSubscripts);

	      accessList.insert(make_pair(access_str, this_access));
	      
	      arrList.insert( make_pair(namecomp, accessList));
	    }
	  else//we found the entry for the component 
	    {
	      //get array's access list 
	      AccessListPerArray_t* accList=&((arrayCompPair)->second);  

	      //do we have the entry for the access?
	      AccessListPerArray_t::iterator acIt = accList->find(access_str);

	      //no, this access is new
	      if(acIt == accList->end())
		{
		  computeDependentLoopVar(this_access, unifiedSubscripts);
		  computeOffsets(this_access, unifiedSubscripts);

		  accList->insert(make_pair(access_str, this_access));
		}
	      else
		{ //we have seen this access before, update the read or write count
		  AccessPattern_t* existing_access = &(acIt->second);
		  if(this_access.readCount == 1)
		    existing_access->readCount++;
		  else if(this_access.writeCount == 1)
		    existing_access->writeCount++;
		}
	    } 
      }
  }

  for(list<SgNode*>::iterator f_it= flatArrayRefList.begin(); f_it != flatArrayRefList.end(); f_it++)
  {
      SgExpression* exp = isSgExpression(*f_it);
      ROSE_ASSERT(exp);

      SgExpression* arrNameExp = NULL;
      std::vector<SgExpression*> subscripts; 

      if(isArrayReference(exp, &arrNameExp, &subscripts)) //no else 
	{

          //expand subscripts
          int indexNo = 0 ;
          vector<SgExpression*>::iterator it;
          for(it  = subscripts.begin(); it != subscripts.end(); it++, indexNo++)
          {
            SgExpression* index = (*it);
            SgNode* parent= index->get_parent();
            SgStatement *enclosingStmt = getEnclosingStatement(index);
            if(enclosingStmt){
              SgBasicBlock *bb = isSgBasicBlock(enclosingStmt->get_parent());
              if(bb){
                SgStatementPtrList stmtList = bb->get_statements();
                recursivePropagation(stmtList, stmtList.begin(), enclosingStmt, index, subscripts, indexNo);
              }
            }
          }

	  SgInitializedName* array_name = convertRefToInitializedName (arrNameExp);	  
	  ROSE_ASSERT(array_name);

	  //assumes the component is in the last subscript 
	  string component  = getComponent(subscripts);

	  //create a array name - component pair 
	  NameComp_t namecomp (array_name, component);

	  //relies on the expression order that AST query returns
	  //first statements are handled first 
	  firstWrittenOrReadAnalysis(firstAccessMap, exp, namecomp);

	  //compute attributes related to this access pattern
	  AccessListPerArray_t accessList;
	  AccessPattern_t this_access;

	  string access_str = getAccessString(subscripts);

	  //if the exp is part of a function call, we 
	  if(isPartofFuncCallExp(exp)){ 
	    this_access.readCount = 0 ; 
	    this_access.writeCount = 0 ; 
	  }
	  else {
	    this_access.readCount = isRead(exp) ? 1 : 0 ;
	    this_access.writeCount = (this_access.readCount == 0) ? 1 : 0 ;
	  }

	  ArraysAccessInfo_t::iterator arrayCompPair = arrList.find(namecomp);	  

	  //this is the first time the array-comp pair  appears in the node
	  if( arrayCompPair == arrList.end())
	    {
	      computeDependentLoopVar(this_access, subscripts);
	      computeOffsets(this_access, subscripts);

	      accessList.insert(make_pair(access_str, this_access));
	      
	      arrList.insert( make_pair(namecomp, accessList));
	    }
	  else//we found the entry for the component 
	    {
	      //get array's access list 
	      AccessListPerArray_t* accList=&((arrayCompPair)->second);  

	      //do we have the entry for the access?
	      AccessListPerArray_t::iterator acIt = accList->find(access_str);

	      //no, this access is new
	      if(acIt == accList->end())
		{
		  computeDependentLoopVar(this_access, subscripts);
		  computeOffsets(this_access, subscripts);

		  accList->insert(make_pair(access_str, this_access));
		}
	      else
		{ //we have seen this access before, update the read or write count
		  AccessPattern_t* existing_access = &(acIt->second);
		  if(this_access.readCount == 1)
		    existing_access->readCount++;
		  else if(this_access.writeCount == 1)
		    existing_access->writeCount++;
		}
	    } 
	  
	} //end of isArrayRefence
    } //end of pntrefexp 
}


bool AccessAnalysis::isRead(SgExpression* exp)
{
  //not very robust but works for now. 
  SgStatement* stmt = getEnclosingStatement(exp);
 
  //ROSE_ASSERT(stmt);
  if(isSgExprStatement(stmt))
    {
      SgExpression* expStmt =isSgExprStatement( stmt)->get_expression();
      ROSE_ASSERT(expStmt);

      if(isSgAssignOp(expStmt))
	{
	  SgExpression * lhs = isSgAssignOp(expStmt)->get_lhs_operand();
	  
	  if(lhs == exp)
	    return false;

	  if(isSgDotExp(lhs)) //A.x , return false for x and A
	    {
	      SgExpression* field = isSgDotExp(lhs)-> get_rhs_operand();
	      SgExpression* struct_= isSgDotExp(lhs)-> get_lhs_operand();

	      if(field == exp || struct_ == exp){
		return false; //it is a write
	      }
	    }
	}
    }
  else if (isSgFortranDo(stmt))
    {// do n = nspecies : return n as readwrite var
      SgExpression * init = isSgFortranDo(stmt)->get_initialization();

      if(isSgAssignOp(init))
	{
	  SgExpression * lhs = isSgAssignOp(init)->get_lhs_operand();
	  SgExpression * rhs = isSgAssignOp(init)->get_rhs_operand();
	  
	  if(lhs == exp)
	    return false;
	  else 
	    return true;
	}
    }
  //assumes it is a read if we cannot figure out
  return true;
}


int expandClause(SgBinaryOp* op){ //this routine is useful when the lhs and/or rhs of a multiply op is not a varRef
 SgNode* parent = op->get_parent();
 SgExpression* lhs = op->get_lhs_operand(); 
 SgExpression* rhs = op->get_rhs_operand(); 
 if(isSgVarRefExp(lhs)){
   if(isSgAddOp(rhs)){
     SgMultiplyOp* newLhs= buildMultiplyOp(copyExpression(lhs), copyExpression(isSgBinaryOp(rhs)->get_lhs_operand()));
     SgMultiplyOp* newRhs= buildMultiplyOp(copyExpression(rhs), copyExpression(isSgBinaryOp(rhs)->get_rhs_operand()));
     SgAddOp* newOp = buildAddOp(newLhs, newRhs);
     newOp->set_parent(parent);
     delete op;
     op= newOp;
     return 1;
   }
   if(isSgAddOp(lhs)){
     SgMultiplyOp* newLhs= buildMultiplyOp(copyExpression(isSgBinaryOp(lhs)->get_lhs_operand()),copyExpression(rhs));
     SgMultiplyOp* newRhs= buildMultiplyOp(copyExpression(isSgBinaryOp(lhs)->get_rhs_operand()),copyExpression(rhs));
     SgAddOp* newOp = buildAddOp(newLhs, newRhs);
     newOp->set_parent(parent);
     delete op;
     op= newOp;
     return 1;
   }
   if(isSgSubtractOp(rhs)){
     SgMultiplyOp* newLhs= buildMultiplyOp(copyExpression(lhs), copyExpression(isSgBinaryOp(rhs)->get_lhs_operand()));
     SgMultiplyOp* newRhs= buildMultiplyOp(copyExpression(rhs), copyExpression(isSgBinaryOp(rhs)->get_rhs_operand()));
     SgSubtractOp* newOp = buildSubtractOp(newLhs, newRhs);
     newOp->set_parent(parent);
     delete op;
     op= newOp;
     return -1;
   }
   if(isSgSubtractOp(lhs)){
     SgMultiplyOp* newLhs= buildMultiplyOp(copyExpression(isSgBinaryOp(lhs)->get_lhs_operand()),copyExpression(rhs));
     SgMultiplyOp* newRhs= buildMultiplyOp(copyExpression(isSgBinaryOp(lhs)->get_rhs_operand()),copyExpression(rhs));
     SgSubtractOp* newOp = buildSubtractOp(newLhs, newRhs);
     newOp->set_parent(parent);
     delete op;
     op= newOp;
     return -1;
   }
 }
}

void AccessAnalysis::normalizeTerms(SgNode *node, bool sign){
 if(node==NULL) return;
 //we haven't processed composit operations such as dot and arrow
 Rose_STL_Container<SgNode*> opList= NodeQuery::querySubTree(node, V_SgArrowExp); 
 if(opList.size()>0) return;
 Rose_STL_Container<SgNode*> opList1= NodeQuery::querySubTree(node, V_SgDotExp); 
 if(opList1.size()>0) return;
 
 if(isSgVarRefExp(node)){  
   clause *c= new clause;
   c->sign= sign;
   c->coeff=1;
   c->lhs= isSgVarRefExp(node);
   c->rhs=NULL; 
   c->loopVar= false; 
   clauseList.push_back(c);
   return;
 }
 if(isSgIntVal(node)){
   clause *c= new clause;
   c->sign= sign;
   c->coeff=1;
   c->lhs= isSgIntVal(node);
   c->rhs=NULL;
   c->loopVar= false; 
   clauseList.push_back(c);
   return;
 }
 if(isSgMultiplyOp(node)){
   if(!isSgVarRefExp(isSgBinaryOp(node)->get_lhs_operand())&&!isSgIntVal(isSgBinaryOp(node)->get_lhs_operand()) || !isSgVarRefExp(isSgBinaryOp(node)->get_rhs_operand())&&!isSgIntVal(isSgBinaryOp(node)->get_rhs_operand()))
   {
     if(isSgAddOp(isSgBinaryOp(node)->get_lhs_operand()) || isSgSubtractOp(isSgBinaryOp(node)->get_lhs_operand()) || isSgAddOp(isSgBinaryOp(node)->get_rhs_operand()) || isSgSubtractOp(isSgBinaryOp(node)->get_rhs_operand())){
       int returnedSign= expandClause(isSgBinaryOp(node));
       normalizeTerms(isSgBinaryOp(node)->get_lhs_operand(), sign);
       normalizeTerms(isSgBinaryOp(node)->get_lhs_operand(), sign && returnedSign);
       return;
     }
     if(isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())){
       clause *c= (clause*)malloc(sizeof(clause));
       c->sign= sign;
       c->loopVar= false; 
       if(isSgIntVal(isSgBinaryOp(node)->get_rhs_operand())){
         c->coeff= isSgIntVal(isSgBinaryOp(node)->get_rhs_operand())->get_value();
         c->lhs= isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_lhs_operand();
         c->rhs= isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_rhs_operand(); 
         clauseList.push_back(c);
         return;
       }
       if(isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_lhs_operand())){
         c->coeff= isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_lhs_operand())->get_value();
         c->lhs= isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_rhs_operand();
         c->rhs= isSgBinaryOp(node)->get_rhs_operand(); 
         clauseList.push_back(c);
         return;
       }
       if(isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_rhs_operand())){
         c->coeff= isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_rhs_operand())->get_value();
         c->lhs= isSgMultiplyOp(isSgBinaryOp(node)->get_lhs_operand())->get_lhs_operand();
         c->rhs= isSgBinaryOp(node)->get_rhs_operand(); 
         clauseList.push_back(c);
         return;
       }
     }
     if(isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())){
       clause *c= (clause*)malloc(sizeof(clause));
       c->sign= sign;
       c->loopVar= false; 
       if(isSgIntVal(isSgBinaryOp(node)->get_lhs_operand())){
         c->coeff= isSgIntVal(isSgBinaryOp(node)->get_lhs_operand())->get_value();
         c->lhs= isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_lhs_operand();
         c->rhs= isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_rhs_operand();
         clauseList.push_back(c);
         return;
       }
       if(isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_lhs_operand())){
         c->coeff= isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_lhs_operand())->get_value();
         c->lhs= isSgBinaryOp(node)->get_lhs_operand();
         c->rhs= isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_rhs_operand();
         clauseList.push_back(c);
         return;
       }
       if(isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_rhs_operand())){
         c->coeff= isSgIntVal(isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_rhs_operand())->get_value();
         c->lhs= isSgBinaryOp(node)->get_lhs_operand();
         c->rhs= isSgMultiplyOp(isSgBinaryOp(node)->get_rhs_operand())->get_lhs_operand();
         clauseList.push_back(c);
         return;
       }
     }
     ROSE_ASSERT(false); //we haven't implemented the other cases 
   }
   SgExpression* lhs= isSgBinaryOp(node)->get_lhs_operand();
   SgExpression* rhs= isSgBinaryOp(node)->get_rhs_operand();
   clause *c= (clause*)malloc(sizeof(clause));
   c->sign= sign;
   c->coeff=1;
   if(isSgIntVal(lhs)) c->coeff=isSgIntVal(lhs)->get_value();
   else if(isSgIntVal(rhs)) c->coeff=isSgIntVal(rhs)->get_value();
   c->lhs= lhs;
   c->rhs= rhs; 
   c->loopVar= false; 
   clauseList.push_back(c);
   return;
 }
 if(isSgAddOp(node)) {
   SgExpression* lhs = isSgBinaryOp(node)->get_lhs_operand();
   SgExpression* rhs = isSgBinaryOp(node)->get_rhs_operand();
   if(lhs)normalizeTerms(isSgNode(lhs), sign);
   if(rhs)normalizeTerms(isSgNode(rhs), sign);
   return;
 }
 if(isSgSubtractOp(node)) {
   SgExpression* lhs = isSgBinaryOp(node)->get_lhs_operand();
   SgExpression* rhs = isSgBinaryOp(node)->get_rhs_operand();
   normalizeTerms(lhs, sign);
   normalizeTerms(rhs, sign==true?false:true);
   return;
 }
 if(isSgMinusOp(node)){
   SgExpression* operand = isSgUnaryOp(node)->get_operand();
   normalizeTerms(operand, sign==true?false:true);
   return;
 }
}

//void AccessAnalysis::sortVarList(){
     
//}

//Our analysis supports 3 common memory access representations (generated by ROSE frontend)
//Fortran representation: A(i,j,k) consists of 1 array name and 3 subscripts
//C/C++ representation: A[k][j][i] consists of 3 array names, each going with a subscript
//Flattened array representation A[ijk] consists of 1 array name and 1 subscript

void AccessAnalysis::computeOffsets(AccessPattern_t& this_access,
				    vector<SgExpression*> subscripts)
{
  //subscripts go from fastest varying dim to slowest
  //but loops are from outmost to inner most 	
  int indexNo = 0 ;

  vector<SgExpression*>::iterator it;
  //sequentially analyze all subscripts
  for(it = subscripts.begin(); it != subscripts.end(); it++, indexNo++)
    {
      //set the default values for this access
      this_access.offset.push_back("");

      SgExpression* index = (*it);
      ROSE_ASSERT(index);
      Rose_STL_Container<SgNode*> varList = NodeQuery::querySubTree(index, V_SgVarRefExp);
      
      if(varList.size() ==  0)//no varRef in the subscript
	{ 
	  //constant index, this code handles all 3 memory access schemes
	  if(isSgIntVal(index))
	    this_access.offset[indexNo] = index->unparseToString();
	}
      else if(varList.size() >= 1)
	{
          eraseClauseList();
          //expand the subscript expression to sum of products
          normalizeTerms(isSgNode(index), true);
	  for(vector<SgInitializedName*>::iterator it= loopIndices->begin(); it!= loopIndices->end(); it++){
            SgInitializedName* loopName = (*it);
	    if(loopName != NULL){
              for(list<clause*>::iterator it1= clauseList.begin(); it1!=clauseList.end(); it1++){
	        clause* c1 = (*it1); 
	        if(isSgVarRefExp(c1->lhs)){
		  SgSymbol *s= isSgVarRefExp(c1->lhs)->get_symbol();
		  if(strcmp(loopName->get_name().getString().c_str(), s->get_name().getString().c_str())==0)
		    c1->loopVar=true;
		}
                if(isSgVarRefExp(c1->rhs)){
                  SgSymbol *s= isSgVarRefExp(c1->rhs)->get_symbol();
                  if(strcmp(loopName->get_name().getString().c_str(), s->get_name().getString().c_str())==0)
                    c1->loopVar=true;
                }
	      }
            }
          }
#if 1
          std::map<string, int> countMap;

          bool loopVarMatched=false;
          //check if a varRef matches with one of the loop variables
	  for(vector<SgInitializedName*>::iterator it= loopIndices->begin(); it!= loopIndices->end(); it++){
	    string strideName= "";
            SgInitializedName* loopName = (*it);
	    if(loopName != NULL){
              for(list<clause*>::iterator it1= clauseList.begin(); it1!=clauseList.end(); it1++){
 	        //first, we check the flattened array style is used. If it is, we should find a stride for slow varying dimensions
 	        //stride can be empty for the most varying dimension
	        clause* c1 = (*it1); 
	        if(isSgVarRefExp(c1->lhs)){
		  SgSymbol *s= isSgVarRefExp(c1->lhs)->get_symbol();
		  if(strcmp(loopName->get_name().getString().c_str(), s->get_name().getString().c_str())==0)
	            if(isSgVarRefExp(c1->rhs)){
		      strideName = isSgVarRefExp(c1->rhs)->get_symbol()->get_name().getString();
		      c1->loopVar=true;
		      loopVarMatched=true;
                   }
	        }
	        if(isSgVarRefExp(c1->rhs)){
	  	  SgSymbol *s= isSgVarRefExp(c1->rhs)->get_symbol();
		  if (strcmp(loopName->get_name().getString().c_str(), s->get_name().getString().c_str())==0 )
	            if(isSgVarRefExp(c1->lhs)){
		      strideName = isSgVarRefExp(c1->lhs)->get_symbol()->get_name().getString();
		      c1->loopVar=true;
		      loopVarMatched=true;
		    }
                }
	        if(strideName!=""){ //we found a stride
		  //now we search all clauses to find offsets associated with this stride
		  //for example -k*kStride + 2*kStride where k is loop variable => kStride is stride for k and offset = -2
		 if(countMap.find(strideName)== countMap.end()){
		  int count=0;
                  for(list<clause*>::iterator it2= clauseList.begin(); it2!=clauseList.end(); it2++){
		    clause* c2= (*it2);
	            if(isSgVarRefExp(c2->lhs)){
		      if(c1!=c2 && strcmp(isSgVarRefExp(c2->lhs)->get_symbol()->get_name().getString().c_str(), strideName.c_str())==0){
		        c2->loopVar=true;
	                if(c2->rhs == NULL){
		          if(c2->sign)count++;
		  	  else count--;
                        }else if(isSgIntVal(c2->rhs)){
		   	  SgIntVal* value = isSgIntVal(c2->rhs);
			  if(c2->sign)count+= atoi(value->unparseToString().c_str());
			  else count-= atoi(value->unparseToString().c_str());
		        }
		        //else if (isSgVarRefExp(c2->lhs)) this_access.offset[indexNo] += ;
		      }
		    }
	            if(isSgVarRefExp(c2->rhs)){
		      if(c1!=c2 && strcmp(isSgVarRefExp(c2->rhs)->get_symbol()->get_name().getString().c_str(), strideName.c_str())==0){
		        c2->loopVar=true;
	                if(c2->lhs == NULL){
		          if(c2->sign)count++;
		          else count--;
                        }else{ 
			  if(isSgIntVal(c2->lhs)){
		          SgIntVal* value = isSgIntVal(c2->lhs);
      	                  if(c2->sign)count+= atoi(value->unparseToString().c_str());
			    else count-= atoi(value->unparseToString().c_str());
			  }
                        }
                      }
		    }
		  }
                  countMap[strideName]=count;
                 }
	        }//end found a stride
              }//end loop over clauses
	    }//end found loop name 
//process the infomation we have learned
            if((it+1)== loopIndices->end()) 
	    {
	      if(this_access.offset[indexNo] != "")this_access.offset[indexNo] = "," + this_access.offset[indexNo];
  	      bool modified=false;
              for(list<clause*>::reverse_iterator it1= clauseList.rbegin(); it1!=clauseList.rend(); it1++){
	        clause* c1 = (*it1); 
	        if(c1->loopVar==false)
	        {
	          if(c1->rhs){
                    ostringstream offset_str ;
		    offset_str << "*";  
                    offset_str << (c1->rhs)->unparseToString().c_str();
		    if(this_access.offset[indexNo] != "" && modified)
		      offset_str << "+";  
                    this_access.offset[indexNo] = offset_str.str() + this_access.offset[indexNo]; 
		    modified=true;
		  }
		  if(!c1->sign){
                    ostringstream offset_str; 
                    offset_str << "-" ; //minus sign
                    offset_str << (c1->lhs)->unparseToString().c_str(); 
		    if(this_access.offset[indexNo] != "" && modified && !c1->rhs)
		      offset_str << "+";  
                    this_access.offset[indexNo] = offset_str.str() + this_access.offset[indexNo]; 
		    modified=true;
		  }else{
                    ostringstream offset_str ; 
                    offset_str << (c1->lhs)->unparseToString().c_str(); 
		    if(this_access.offset[indexNo] != "" && modified && !c1->rhs)
		      offset_str << "+";  
                    this_access.offset[indexNo] = offset_str.str() + this_access.offset[indexNo]; 
		    modified=true;
		  }
	        }
	      }
	      if(!modified) {
	        this_access.offset[indexNo] = "0" + this_access.offset[indexNo];
	      }
            }else{
	      if(strideName != ""){
	        //if(varList.size()>1)
	          //this_access.offset[indexNo] += ", ";
                char countChar[10];
	        sprintf(countChar, "%d", countMap[strideName]);
	        if(this_access.offset[indexNo] != "")this_access.offset[indexNo] = "," + this_access.offset[indexNo];
	        this_access.offset[indexNo] = string(countChar)+this_access.offset[indexNo]; 
	      }else{
#if 0 
              if(loopVarMatched){
	        if(varList.size()>1)
	          this_access.offset[indexNo] += ", ";
	        this_access.offset[indexNo] += string(countChar); 
	      }
#endif
	      }
            }
          }
#endif
        }//varList.size >= 1
      } //end of for loop
}

SgInitializedName* AccessAnalysis::getNonNullLoopIndex(int position)
{
  if(position < 0)
    return NULL; 

  SgInitializedName* loop_iname = loopIndices->at(position);
  if(loop_iname == NULL){
    return getNonNullLoopIndex(position-1);
  }
  return loop_iname; 
}




void AccessAnalysis::computeDependentLoopVar(AccessPattern_t& this_access,
					    vector<SgExpression*> &subscripts)
				
{
  //subscripts go from fastest varying dim to slowest
  //but loops are from outmost to inner most 	
  int indexNo = 0 ;

  vector<SgExpression*>::iterator it;
  for(it  = subscripts.begin(); it != subscripts.end(); it++, indexNo++)
    {

      //set the default values for this access
      this_access.dependentLoopVar.push_back(""); 

      SgExpression* index = (*it);
      ROSE_ASSERT(index);
      eraseClauseList();
      normalizeTerms(isSgNode(index), true);
      Rose_STL_Container<SgNode*> varList = NodeQuery::querySubTree(/*index*/subscripts[indexNo], V_SgVarRefExp);

      if(varList.size() ==  0)
	{ //constant index, default values: loop dependentLoopVar is " " because it doesn't depend on the loop
	  //TODO: no VarRefExp may not mean constant index exp, it might be func(q(:,:,:,n))
	  //do nothing for now
	}
      else if(varList.size() >= 1)
	{
	  bool foundOneLoopVar=false;
          bool firstLoopVar=true;
	  std::vector<SgInitializedName*>::iterator indexItr;
	  for(indexItr= loopIndices->begin(); indexItr!= loopIndices->end(); indexItr++){
            SgInitializedName* loopName = (*indexItr);
	    if(loopName != NULL){
              string coeff="";
              for(list<clause*>::iterator it1= clauseList.begin(); it1!=clauseList.end(); it1++){
	        clause* c1 = (*it1); 
	        if(isSgVarRefExp(c1->lhs)){
		  SgSymbol *s= isSgVarRefExp(c1->lhs)->get_symbol();
		  if(strcmp(loopName->get_name().getString().c_str(), s->get_name().getString().c_str())==0){
                    foundOneLoopVar=true;
                    char buffer [10];
                    sprintf (buffer, "%d * ", c1->coeff);
                    if(c1->coeff!=1) coeff+= buffer;
                    break;
                  }
                }
	        if(isSgVarRefExp(c1->rhs)){
	  	  SgSymbol *s= isSgVarRefExp(c1->rhs)->get_symbol();
		  if (strcmp(loopName->get_name().getString().c_str(), s->get_name().getString().c_str())==0){
                    foundOneLoopVar=true;
                    char buffer [10];
                    sprintf (buffer, "%d * ", c1->coeff);
                    if(c1->coeff!=1) coeff+= buffer;
                    break;
                  }
                }
              }
	    if(foundOneLoopVar){//index variable matches with the loop control variable 
	      //get the loop no that the variable depends on
	      //outermost loop is at position zero 
	      int position = indexItr - loopIndices->begin();
	      ROSE_ASSERT(position >= 0);
	      SgInitializedName* loopName = loopIndices->at(position);
	      ROSE_ASSERT(loopName);
              if(!firstLoopVar) //varList.size() > 1 && indexItr!=(loopIndices->end()-1))
  	        this_access.dependentLoopVar[indexNo]+= ",";
	      this_access.dependentLoopVar[indexNo] += coeff;
	      this_access.dependentLoopVar[indexNo] += loopName->get_name().str();
	      foundOneLoopVar=false;
              firstLoopVar=false;
	    }//loop variable and index variable match
         }
       }
     }




#if 0
         Rose_STL_Container<SgNode*>::iterator varListIt= varList.begin();
	 bool foundOneLoopVar=false;
         while(varListIt!= varList.end()){ 
	  SgVarRefExp* indexVar = isSgVarRefExp(*varListIt);
	  ROSE_ASSERT(indexVar);
	  
	  SgInitializedName* iname = convertRefToInitializedName (indexVar);
	  ROSE_ASSERT(iname);
	  
	  std::vector<SgInitializedName*>::iterator indexItr;
	  //indexItr = std::find(this->loopIndices->begin(), this->loopIndices->end(), iname);
          for(indexItr= this->loopIndices->begin();indexItr!= this->loopIndices->end(); indexItr++){
	    if(*indexItr)
	    if(strcmp((*indexItr)->get_name().getString().c_str(),iname->get_name().getString().c_str())==0){
	       break;
            }
	  }
	  
	  if(indexItr == loopIndices->end())
	    {
	    }
	  else {//index variable matches with the loop control variable 
	    //get the loop no that the variable depends on
	    //outermost loop is at position zero 
	    int position = indexItr - loopIndices->begin();
	    ROSE_ASSERT(position >= 0);
	    SgInitializedName* loopName = loopIndices->at(position);
	    ROSE_ASSERT(loopName);
            if(foundOneLoopVar==true) //varList.size() > 1 && indexItr!=(loopIndices->end()-1))
  	      this_access.dependentLoopVar[indexNo]+= ",";
	    this_access.dependentLoopVar[indexNo] += loopName->get_name().str();
	    foundOneLoopVar=true;
	  }//loop variable and index variable match
          varListIt++;
         }
	}//varList.size >= 1
      //}
#endif

    } //end of for loop
}

string AccessAnalysis::getComponent(vector<SgExpression*> subscripts) 
{
  //this function assumes the component is in the last subs
  //if there are fewer than DEF_NUM_DIM loops, then return an empty component
 
  string component = "";
  
  //subscripts go from fastest varying dim to slowest
  //but loops are from outmost to inner most 	
  if(subscripts.size() > DEF_NUM_DIM) {

    vector<SgExpression*>::iterator it = subscripts.end() -1 ; 
    SgExpression* index = (*it);
    ROSE_ASSERT(index);
    
    //might be something like this (a,b,c,QN:)
    if(isSgSubscriptExpression(index))
      {
	SgSubscriptExpression* sub_expr = isSgSubscriptExpression(index);
	SgExpression* ub = sub_expr->get_upperBound() ;
	SgExpression* lb = sub_expr->get_lowerBound() ;
	SgExpression* st = sub_expr->get_stride() ;
      }
    //component check 
    string indexStr = index->unparseToString(); 
    
    component.resize(indexStr.size());
    component = indexStr;
  }

  return component; 
}

void AccessAnalysis::firstWrittenOrReadAnalysis(FirstAccessMap_t& firstAccessmap,
						SgExpression* exp, 
						NameComp_t namecomp)
					
{
  //Checking the first reference to this namecomp is a read or write 

  FirstAccessMap_t::iterator fa = firstAccessmap.find(namecomp);
  
  //no entry found, this is the first time we see this array-comp pair 
  if(fa == firstAccessmap.end())
    {
      FirstAccess_t first_entry;
      first_entry.stmt = getEnclosingStatement(exp);

      first_entry.firstRead = isRead(exp);

      firstAccessmap.insert(make_pair(namecomp, first_entry));
    }
  else // an entry found. if it is written, it may actually first read or vice versa
    {
      FirstAccess_t* first_entry= &(fa->second); 
      
      //first written, not read, this might change if the stmts are the same A[][] = A[][]
      if(!first_entry->firstRead && isRead(exp))
	{
	  SgStatement* cur_stmt = getEnclosingStatement(exp);
	  if(first_entry->stmt == cur_stmt)//they are the same statements
	    {
	      //need to change first written into first read
	      first_entry->firstRead = true;
	    }
	}
    }
}
  

// eof















