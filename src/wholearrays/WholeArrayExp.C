/*
  Didem Unat
 */

#include <rose.h>

#include "WholeArrayExp.h"
#include "../utils/ArrayUtils.h"
#include "../utils/Utils.h"
#include "stringify.h"

using namespace std;
using namespace SageInterface;
using namespace Utils; 
using namespace SageBuilder;
using namespace ArrayUtils;


void WholeArrayExp::convertWholeArraysToFortranDos(SgNode* node)
{
  //In order to treat whole array statements as same as fortran do loops, we convert all 
  //the whole array references into nested loop

  Rose_STL_Container<SgNode*> stmts = NodeQuery::querySubTree(node, V_SgStatement);
  Rose_STL_Container<SgNode*>::iterator st ;
  
  set<SgStatement*> selected ;
  
  for(st= stmts.begin(); st != stmts.end(); st++)
    {
      SgStatement* stmt = isSgStatement(*st);
      switch(stmt->variantT())
	{
	  //There are sliced array expressions in a function call, what do we do for those?
	case V_SgFunctionCallExp:
	  break;
	case V_SgAssignStatement:
	case V_SgExprStatement:
	    {
	      Rose_STL_Container<SgNode*> funcCalls = NodeQuery::querySubTree(stmt, V_SgFunctionCallExp);

	      //we don't handle if the whole array statement includes a function call A= funcCall(B)
	      if(funcCalls.size()>0)
		  break;

	      int loopDepth = 0;
	      if(WholeArrayExp::isWholeArrayStatement(stmt, &loopDepth))
		selected.insert(stmt);

	      break;
	    }
	default:
	  break;
	}//end of switch
    }//end of for loop

  //Now, we need to convert the whole array statements into fortran do statement
  for(set<SgStatement*>::iterator stmt= selected.begin(); stmt!= selected.end(); stmt++)
    {
      convertWholeArrayStatementToFortranDoStatement(*stmt); 
    }

  selected.clear();
}

void WholeArrayExp::convertWholeArrayStatementToFortranDoStatement(SgStatement* stmt)
{
  /*
    Convert the following statement into:
    A = B + C 

    do j0 = 0, ubound(A, 2) - lbound(A, 2), stride
    jA = ...
      do i0 = 0, ubound(A, 1) - lbound(A, 1), stride  
        iA = lbound(A, 1) + i0
        iB = lbound(B, 1) + i0 
        ic = lbound(C, 1) + i0 
        A(iA, jA, k) = B (iB, ...) + C (iC)
      end do
    end do
  
    TODO: how about ascending order ?
   */
  
  //step 1: find all the array references in the selected statement that needs to be changed 
  //step 2: find the lower, upper bounds and stride of each whole array expression
  //step 3: build the do loop by using the lower, upper and stride expressions (pick one of the array to do that)
  //step 4: create index variables for each array using their lowers and loop index vars
  //step 5: expand the array expressions using the index variables created in step 4

  //Step 1. find all the whole array references

  int depth = -1; 
  vector<SgExpression*> selectedRefs; 
  
  Rose_STL_Container<SgNode*> varList = NodeQuery::querySubTree(stmt, V_SgExpression);
  Rose_STL_Container<SgNode*>::iterator var;
  for(var = varList.begin(); var != varList.end(); var++)
    {
      SgExpression* varExp = isSgExpression(*var);
      int dim = 0;
      if(isWholeArrayExpression(varExp, &dim))
	{
	  //this is a whole array expression so add it to the selected ones
	  //note that we find ALL
	  selectedRefs.push_back(varExp);

	  if(depth == -1)
	      depth = dim ; 
	  // all the arrays (sliced arrays) should have the same rank
	  ROSE_ASSERT(depth == dim);
	}
    }
  //step 2: find the lower, upper bounds and stride of each array expression
  //inner most first
  vector<SgExpression*>::iterator sel; 
  LoopBounds_t* loopBounds = new LoopBounds_t[selectedRefs.size()]; 
  for (sel= selectedRefs.begin(); sel != selectedRefs.end(); sel++) 
    {
      SgExpression* selExp = (*sel);
      RangeInfo_t lInfo; 

      getLoopBounds(selExp, lInfo.lower, lInfo.upper, lInfo.stride);
      loopBounds->insert(make_pair(selExp, lInfo));
    }

  //step 3: build the do loop by using the lower, upper and stride expressions
  SgScopeStatement* scope = stmt->get_scope();
  ROSE_ASSERT(scope);
  SgBasicBlock* newLoopBlock = buildBasicBlock();

  insertStatementAfter(stmt, newLoopBlock, false); //false means don't move the preprocessing info, otherwise I get a warning
  //newLoopBlock-> set_file_info(scope->get_file_info());
  newLoopBlock-> set_file_info(stmt->get_file_info());

  //we are going to use the bounds of the first selected expression to create forloop bounds 
  LoopBounds_t::iterator l= loopBounds-> begin(); 
  if( l!= loopBounds->end())
    {
      RangeInfo_t* lInfo = &(l->second);

      buildFortranDoLoop(newLoopBlock, lInfo);
    }

  // //Query all the Do loops that we generated in the body
  vector<SgFortranDo*> loopList; 
  Rose_STL_Container<SgNode*> loops = NodeQuery::querySubTree(newLoopBlock, V_SgFortranDo);
  for(Rose_STL_Container<SgNode*>::reverse_iterator loop_i = loops.rbegin(); loop_i != loops.rend(); loop_i++){

      SgFortranDo* loop = isSgFortranDo(*loop_i);
      ROSE_ASSERT(loop);
      loopList.push_back(loop);
    }

  //step 4 
  //LoopBounds is a whole array exp map
  // exp -> lower, upper and stride 
  // get their loop index variable
  // create an index variable for each selected whole array expressions
  // add a declaration for each of these index variables
  // initialize these index variables as iA = lower(A) + i
#if 0
  SgFunctionDefinition* func = getEnclosingFunctionDefinition(stmt);
  SgBasicBlock* func_body = func->get_body();
  ROSE_ASSERT(func_body);

  for(LoopBounds_t::iterator lmap= loopBounds-> begin(); lmap != loopBounds-> end() ; lmap++ )
    {
      SgExpression* arrayExp = lmap-> first; 
      RangeInfo_t* lInfo = &(lmap->second);
      
      vector<SgExpression*>* lower = &(lInfo->lower);
      vector<SgExpression*>* upper = &(lInfo->upper);
      vector<SgExpression*>* stride = &(lInfo->stride);
      vector<SgVariableDeclaration*>* indexVar= &(lInfo->indexVar);

      int indexNo = 0;
      vector<SgExpression*>::iterator lb= lower-> begin(); 
      vector<SgExpression*>::iterator st= stride-> begin(); 
      vector<SgExpression*>::iterator up= upper -> begin(); 

      for(; lb != lower-> end(); lb++, st++, up++)
  	{
  	  if( *lb == NULL) { // explicit index 
  	    indexVar-> push_back(NULL);
  	  }
  	  else {

  	    SgVariableDeclaration* var_i = buildIndexVariable(func_body);
  	    indexVar-> push_back(var_i);

  	    SgFortranDo* fortranDo = loopList[indexNo];

  	    SgInitializedName* loopvar = getLoopIndexVariable(isSgNode(fortranDo));

  	    //SgVarRefExp* indexVar = buildVarRefExp(loopvar);
  	    //ROSE_ASSERT(indexVar);
  	    //SgExpression* boundExp = buildAddOp(*lb, indexVar);  //lower(A) + i  
  	    indexVar-> push_back(buildVarRefExp(loopvar)); //Tan added
  	    //SgStatement* init = buildAssignStatement(buildVarRefExp(var_i), boundExp); // iA = lower(A) + i

  	    //SgBasicBlock* doBody = fortranDo-> get_body();
  	    //ROSE_ASSERT(doBody);

  	    //prependStatement(init, doBody);

   	    indexNo++;

  	    //we do not handle strides other than 1, report if that's the case
  	    if((*st) -> unparseToString() != "1"){
  	      cerr << " Whole Array Conversion handles only arrays with stride 1 " << endl; 
  	      cerr<<  " Sorry, we cannot handle this expression -> " << arrayExp-> unparseToString() << " Maybe one day, we will ..."<<endl;
  	    }
	    //upper and stride are not used, so deep delete them to avoid dangling pointers
	    deepDelete(*st);
	    deepDelete(*up);
  	  }
  	}
    }
#endif

  //step 5: expand whole array expressions to explicit array expression using the indices we created in step 4
  std::map<SgExpression*, SgExpression*> replaceList; 

  for(LoopBounds_t::iterator l= loopBounds-> begin(); l != loopBounds->end(); l++)
    {
      SgExpression* arrayExp = l-> first; 
      RangeInfo_t* lInfo = &(l->second);
      vector<SgExpression*>* lower = &(lInfo->lower);
      vector<SgExpression*>* upper = &(lInfo->upper);
      vector<SgExpression*>* stride = &(lInfo->stride);
      vector<SgExpression*>::iterator lb= lower-> begin(); 
      vector<SgExpression*>::iterator st= stride-> begin(); 
      vector<SgExpression*>::iterator up= upper -> begin(); 


      //vector<SgVariableDeclaration*>* indexVar = &(lInfo-> indexVar);
     
      if(isSgPntrArrRefExp(arrayExp))// this is the partial form A(:, :, k)
  	{
  	  //need to replace (:,j,k) with (i,j,k), we need to keep i and j as it is
  	  SgPntrArrRefExp* ptrExp = isSgPntrArrRefExp(arrayExp);
  	  SgExpression* subExp = ptrExp -> get_rhs_operand();

  	  SgExpressionPtrList& list = isSgExprListExp(subExp)->get_expressions();	
  	  //vector<SgVariableDeclaration*>:: iterator ind = indexVar-> begin(); 

  	  //ROSE_ASSERT(indexVar-> size() == list.size());

          int indexNo = 0;
  	  for( SgExpressionPtrList::iterator it = list.begin(); it != list.end(); it++,lb++, st++, up++/*, ind++*/) 
  	    {
  	      SgExpression* exp = (*it);
  	      //SgVariableDeclaration* indVar = *ind ; 

  	      if(exp->variantT() == V_SgSubscriptExpression && *lb != NULL/*&& indVar != NULL*/){
                SgFortranDo* fortranDo = loopList[indexNo];
                indexNo++;
                SgInitializedName* loopvar = getLoopIndexVariable(isSgNode(fortranDo));
                SgVarRefExp* indexVar = buildVarRefExp(loopvar);
                ROSE_ASSERT(indexVar);
  		replaceList.insert(make_pair(exp, buildAddOp(*lb, indexVar) /*buildVarRefExp(indVar)*/));
  	      }
  	    }
  	}// end of partial whole array refs
      else //truely whole array ref, need to create all the indices 
  	{
  	  SgInitializedName* iname = convertRefToInitializedName(arrayExp); 
  	  vector<SgExpression*> indList; 

  	  //vector<SgVariableDeclaration*>:: iterator ind; 
  	  //for (ind = indexVar-> begin();  ind != indexVar-> end(); ind++)
  	    //{
          int indexNo = 0;
          for(; lb != lower-> end(); lb++, st++, up++)
            {
  	      //SgVariableDeclaration* indVar = *ind;
  	      if(*lb != NULL /*indVar != NULL*/){
                SgFortranDo* fortranDo = loopList[indexNo];
                indexNo++;
                SgInitializedName* loopvar = getLoopIndexVariable(isSgNode(fortranDo));
                SgVarRefExp* indexVar = buildVarRefExp(loopvar);
                ROSE_ASSERT(indexVar);
  		indList.push_back(buildAddOp(*lb, indexVar));
	       }
  		//indList.push_back(buildVarRefExp(indVar));
  	    }
  	  SgExprListExp * expr_list = SageBuilder::buildExprListExp(indList);
  	  ROSE_ASSERT(expr_list);
  	  SgPntrArrRefExp* ptrExp = buildPntrArrRefExp(buildVarRefExp(iname), expr_list);	    
  	  ROSE_ASSERT(ptrExp);

  	  replaceList.insert(make_pair(arrayExp, isSgExpression(ptrExp)));

  	} // end of truely whole array refs

    }//end of creating replace list 

  // //FINAL step 
  // //replace the old expression with the new array pointer expression
  std::map<SgExpression*, SgExpression*>::iterator e = replaceList.begin();
  for(; e != replaceList.end(); e++)
    {
      SgExpression* old = e->first; 
      SgExpression* newe = e->second; 

      replaceExpression(old, newe, false); //true means keep the old exp, I should set it to false but it fails if I do
    }

  //append the new body to the loop body
  //This is the inner most 
  SgFortranDo* fDo = *(loopList.begin());
  ROSE_ASSERT(fDo);
  SgBasicBlock* loopBody = fDo -> get_body();
  ROSE_ASSERT(loopBody);

  SgStatement* newLoopBody = deepCopy(stmt);
  appendStatement(newLoopBody, loopBody);

  //remove the old statement
  removeStatement(stmt);

  //clean up 
  for(LoopBounds_t::iterator l= loopBounds-> begin(); l != loopBounds->end(); l++)
    {
      RangeInfo_t* lInfo = &(l->second);
      lInfo->indexVar.clear();
      lInfo->lower.clear();
      lInfo->upper.clear();
      lInfo->stride.clear();
    }
  delete [] loopBounds;
}

void WholeArrayExp::buildFortranDoLoop(SgBasicBlock* body, RangeInfo_t* lInfo)
{ 
  SgFunctionDefinition* func = getEnclosingFunctionDefinition(body);
  ROSE_ASSERT(func);
  SgBasicBlock* func_body = func->get_body();
  ROSE_ASSERT(func_body);
  
  vector<SgExpression*>* lower = &(lInfo-> lower);
  vector<SgExpression*>* upper = &(lInfo-> upper);
  vector<SgExpression*>* stride = &(lInfo-> stride);
  
  //Outmost loop is created first, so use reverse iterator
  vector<SgExpression*>::reverse_iterator lb = lower-> rbegin();
  vector<SgExpression*>::reverse_iterator ub = upper-> rbegin();
  vector<SgExpression*>::reverse_iterator st = stride-> rbegin();
  
  for (; lb != lower->rend(); lb++, ub++, st++)
    {
      if(*lb == NULL ) // this is an explicit loop, skip it 
	continue;
      else 
	{
	  //building loop from lower, upper and stride
	  // do var_i = 0, upper-lower, stride
	  
	  SgVariableDeclaration* var_i = buildIndexVariable(func_body);
	  ROSE_ASSERT(var_i);
	  
	  SgExpression* init = buildAssignOp(buildVarRefExp(var_i), buildIntVal(0));
	  ROSE_ASSERT(init);
	  
	  SgExpression* bound= buildSubtractOp(copyExpression(*ub), copyExpression(*lb));
	  ROSE_ASSERT(bound);
	  
	  SgExpression* increment = copyExpression(*st);
	  ROSE_ASSERT(increment);
	  
	  SgBasicBlock* bb = buildBasicBlock();
	  ROSE_ASSERT(bb != NULL);
	  
	  SgFortranDo* fortranDo = new SgFortranDo(init, bound, increment, bb);
	  //add END DO
	  fortranDo->set_has_end_statement(true);

	  bb->set_parent(fortranDo);
	  //fortranDo-> set_file_info(func_body->get_file_info());
	  fortranDo-> set_file_info(body->get_file_info());
	  
	  appendStatement(fortranDo, body);
	  body = fortranDo -> get_body();
	  body-> set_file_info(fortranDo->get_file_info());
	}
    }
}

SgVariableDeclaration* WholeArrayExp::buildIndexVariable(SgBasicBlock* body)
{

  SgFunctionDefinition* func = getEnclosingFunctionDefinition(body);
  ROSE_ASSERT(func);
  SgBasicBlock* func_body = func->get_body();
  ROSE_ASSERT(func_body);

  string var_name = Utils::generateVarName(func_body);
  
  SgVariableDeclaration* var_i =  buildVariableDeclaration(var_name,
							   buildIntType(), NULL , func_body);
  ROSE_ASSERT(var_i);
  SgStatement* first_stmt = NULL; //getFirstStatement(func_body);      
  
  //ROSE's routine getFirstStatement has a bug, so try getting the first statement 	  
  //we also do not want the use or implicit statements 
  if(first_stmt == NULL || isSgUseStatement(first_stmt) || isSgImplicitStatement(first_stmt)) {
    
    SgStatementPtrList& srcStmts = func_body-> get_statements();
    for (SgStatementPtrList::iterator i = srcStmts.begin(); i != srcStmts.end(); i++)
      {
	SgStatement* s = isSgStatement(*i);
	
	if(isSgUseStatement(s))
	  continue;
	if(isSgImplicitStatement(s))
	  continue;
	else {
	  first_stmt = s; 
	  break;
	}
      }
  }
  ROSE_ASSERT(first_stmt) ; 	      
  insertStatementBefore(first_stmt, var_i, false); 
  
  return var_i;
}

void WholeArrayExp::getLoopBounds(SgExpression* array, 
			       vector<SgExpression*>& lower, 
			       vector<SgExpression*>& upper, 
			       vector<SgExpression*>& stride)
{
  /* 
     Didem: Dec 11, 2013: Swapped the order of the steps, now we do step 2 and 3 first and then merged Step 4 and 1.
     The reason why I changed the order is that I had a problem with Rose functionCallExpressions. When I build the exp and delete
     later, and rebuild it, I get function symbol error. I traced the error back to nondefiningfunctiondeclaration_t in Rose but couldn't
     completely understood what was happening. The order shouldn't make any difference in theory but in practice it does. 

     New version (Dec 11, 2013) tries to get the array range from alloc, array declaration or from the slicing. If these do not work, then get the information 
     by using the default method by calling lbound and ubound (that's last resort)
     
     Step 2= if the array is dynamically allocated, try to get its bounds from the alloc statement

     Step 3= if we can find the dimension information at the array declaration, then use this 

     Step 4= if the array is sliced, we shouldnot use the dimension information at the declaration or at the 
     alloc statement. Parse the subscripts and get the slicing low and highs. Also some indices might be 
     explicit. We don't create a for-loop for these dimensions. (e.g  A(:,:,K) thrid dimension is explicit, set
     its lows and highs to NULL

     Step 1= set the default bounds by calling lbound (ArrayExp, dim) and ubound(ArrayExp, dim)
     this version is not preferred because the index expressions are very complex 
     We only do this if Step 2 and 3 fails and step 4 has a null subscript expression.

     Step 1 and 4 are essential. 2 and 3 are there to simplify the index expressions 
  */

  SgGlobal* global = getGlobalScope(array);
  ROSE_ASSERT(global);

  //init the list with NULL first 
  int DIM = get_array_rank(array);
  for(int index=0; index < DIM ; index++)
    {
      lower.push_back(NULL);
      upper.push_back(NULL);
      stride.push_back(NULL);
    }

  //step 2: look for an alloc statement
  SgInitializedName* iname = convertRefToInitializedName(array); 
  ROSE_ASSERT(iname);
  
  SgDeclarationStatement* var_decl = iname ->get_declaration();
  ROSE_ASSERT(var_decl);
  
  //if it is dynamically allocated, we do not know the size of the array at the declaration
  if(var_decl -> get_declarationModifier().get_typeModifier().isAllocatable())
    {
      SgExpression* exp = getAllocExpression(array);

      if( exp != NULL){ //NULL means couldn't find the alloc statement for this var, it might be in another function or lib

	vector<SgExpression*> subscripts; 
	SgExpression* arrNameExp = NULL;
	
	if(isArrayReference(exp, &arrNameExp, &subscripts))
	  {
	    int index =0;
	    vector<SgExpression*>::iterator sub;
	    
	    for(sub = subscripts.begin(); sub != subscripts.end() ; sub++, index++)
	      parseSubscript(*sub, lower, upper, stride, index);
	  }
      }//exp is NULL, this is the case when the allocation happends in another func or file etc 
    }
  //step 3: Look at the declaration 
  else // if not, we can get the dimension from the declaration
   {
      SgType* type = iname ->get_type();
      ROSE_ASSERT(type);

      if ( !isSgArrayType(type) ) // check if the array is part of a named type struct, class or union
        type = getTypeInNamedType(iname, array); //including class, union or struct

      SgExprListExp* expList = isSgArrayType(type)-> get_dim_info();
      SgExpressionPtrList& list = expList->get_expressions();

      int index=0;
      for( SgExpressionPtrList::iterator it = list.begin(); it != list.end(); it++, index++) 
	parseSubscript(isSgExpression(*it), lower, upper, stride, index);
      
   } // end of declaration 

  /*
    Didem (Dec 11/2013): merged step 1 and 4. Now we use the default case when we cannot get the ranges 
    from other methods such as alloc or declaration 
   */

  //Step 4:
  //Sliced Array references: only some of the dimensions are whole array, then we need to set the explicit ones to NULL 
  //e.g A(:, a:b, k)

  //step 1: set the default loop bounds
  /*
    lower(dim1) = lbound (array, dim1) ...
    upper(dim1) = ubound (array, dim1) ...
   */

  Rose_STL_Container<SgNode*> vList = NodeQuery::querySubTree(array, V_SgPntrArrRefExp);
  if (vList.size() > 0)
    {
      Rose_STL_Container<SgNode*>::iterator v = vList.begin();

      SgPntrArrRefExp* ptrExp = isSgPntrArrRefExp(*v);
      SgExpression* subExp = ptrExp->get_rhs_operand();
     
      SgExpressionPtrList& list = isSgExprListExp(subExp)->get_expressions();	
      int index=0;

      for( SgExpressionPtrList::iterator it = list.begin(); it != list.end(); it++, index++) 
  	{
	  SgExpression* exp = (*it);
	  ROSE_ASSERT(exp);

	  // this index is already explicit, no need to replace it  
	  if(exp->variantT() != V_SgSubscriptExpression) //A(k)
  	    {
	      //if we set a bound already thru alloc or declare, we should delete those 
	      if(lower[index] != NULL)
		deepDelete(lower[index]);
	      if(upper[index] != NULL)
		deepDelete(upper[index]);
	      if(stride[index] != NULL)
		deepDelete(stride[index]);

	      //need to set it to NULL, so that we don't create a loop for this dimension
	      lower[index]=(NULL);
	      upper[index]=(NULL);
	      stride[index]=(NULL);
  	    }
	  else // index is implicit: two options 1) sliced A(a:b) or  2) Null subscript A(:) 
	    {
	      SgSubscriptExpression* sub_expr = isSgSubscriptExpression(exp);
	      SgExpression* ub = sub_expr->get_upperBound() ;
	      SgExpression* lb = sub_expr->get_lowerBound() ;
	      SgExpression* st = sub_expr->get_stride() ;

	      //UPPER BOUND 
	      //if it is null then we can create the default function call exp lbound ubound 
	      if (!isSgNullExpression(ub)){
		if(upper[index] != NULL)
		  deepDelete(upper[index]);		
		upper[index] =(copyExpression(ub));
	      }
	      else //means (:) 
		{
		  if(upper[index] == NULL) { 
		    //fall back to default case
		    //this means we couldn't get the range from alloc or declaration 
		    //so need to use the default function call exp for ubound
		    SgExprListExp* paramListU= buildExprListExp(deepCopy(array), buildIntVal(index+1));
		    ROSE_ASSERT(paramListU);
		    
		    SgExpression* funcExp  = buildFunctionCallExp("ubound", buildIntType(), paramListU, global);
		    ROSE_ASSERT(funcExp);
		    upper[index] = funcExp; 
		  }
		  //else do nothing because we got the info from step 2 or 3
		}

	      //LOWER BOUND
	      if(!isSgNullExpression(lb)){ 
		if(lower[index] != NULL)
		  deepDelete(lower[index]);
		lower[index] =(copyExpression(lb));
	      }
	      else //means (:)
		{
		  if(lower[index] == NULL){ 
		    //fall back to default case 
		    //this means we couldn't get the range from alloc or declaration
		    //need to use the default function call exp for lbound
		    SgExprListExp* paramListL= buildExprListExp(deepCopy(array), buildIntVal(index+1));
		    ROSE_ASSERT(paramListL);
		    
		    SgExpression* funcExp  = buildFunctionCallExp("lbound", buildIntType(), paramListL, global);
		    ROSE_ASSERT(funcExp);
		    lower[index] = funcExp; 
		  }
		  //else do nothing because we got the info from step 2 or 3
		}

	      if(!isSgNullExpression(st)){
		if(stride[index] != NULL)
		  deepDelete(stride[index]);
		stride[index] = (copyExpression(st));
	      }
	      else 
		{
		  if(stride[index] == NULL){ //fall back to default case
		    SgIntVal* val = buildIntVal(1);
		    stride[index] = val; 
		  }
		  //else do nothing 
		}
	    }
  	}
    }// end of sliced array 
}

SgExpression* WholeArrayExp::getAllocExpression(SgExpression* array)
{
  SgInitializedName* iname = convertRefToInitializedName(array); 
  ROSE_ASSERT(iname);

  SgFunctionDefinition* func = getEnclosingFunctionDefinition(array);
  ROSE_ASSERT(func);
  SgBasicBlock* func_body = func->get_body();
  
  Rose_STL_Container<SgNode*> allocList = NodeQuery::querySubTree(func_body, V_SgAllocateStatement);
  Rose_STL_Container<SgNode*>::iterator as;
  
  for(as = allocList.begin(); as != allocList.end(); as++)
    {
      SgAllocateStatement* allocStmt = isSgAllocateStatement(*as); 
      
      SgExprListExp* expList = allocStmt ->get_expr_list();
      SgExpressionPtrList& list = expList->get_expressions();
      
      for( SgExpressionPtrList::iterator it = list.begin(); it != list.end(); it++) 
	{
	  SgExpression* exp = (*it);
	  if(exp->variantT() == V_SgPntrArrRefExp)
	    {
	      vector<SgExpression*> subscripts; 
	      SgExpression* arrNameExp = NULL;
	      if(isArrayReference(exp, &arrNameExp, &subscripts ))
		{
		  SgInitializedName* array_name = convertRefToInitializedName (arrNameExp);
		  ROSE_ASSERT(array_name);
		  
		  if(array_name == iname)
		    {
		      return exp; 
		    }
		}
	    }
	}
    }
  return NULL;
}

void WholeArrayExp::parseSubscript(SgExpression* subExp,     
				 vector<SgExpression*>& lower, 
				 vector<SgExpression*>& upper, 
				 vector<SgExpression*>& stride, int index)
{

  if(lower[index] != NULL)
    deepDelete(lower[index]);
  if(upper[index] != NULL)
    deepDelete(upper[index]);
  if(stride[index] != NULL)
    deepDelete(stride[index]);

  if(subExp->variantT() == V_SgSubscriptExpression)
    {
      SgSubscriptExpression* sub_expr = isSgSubscriptExpression(subExp);
      SgExpression* ub = sub_expr->get_upperBound() ;
      SgExpression* lb = sub_expr->get_lowerBound() ;
      SgExpression* st = sub_expr->get_stride() ;
      
      lower[index]=(copyExpression(lb));
      upper[index]=(copyExpression(ub));
      stride[index]= (copyExpression(st));
    }
  else //A(N) means 1:N
    {
      lower[index]=buildIntVal(1);
      upper[index]=copyExpression(subExp);
      stride[index]=buildIntVal(1);
    }
}

bool WholeArrayExp::isWholeArrayStatement(SgStatement* stmt, int* depth)
{ 
//if there is at least one whole array reference, return true
// 1. either the statement has a (:) subscript expression in it 
// 2. all the arrays in the expression (or statement) are whole arrays references 

  Rose_STL_Container<SgNode*> varExp = NodeQuery::querySubTree(stmt, V_SgExpression);
  Rose_STL_Container<SgNode*>::iterator it;
 
  for(it = varExp.begin(); it != varExp.end(); it++)
    {
      SgExpression* var = isSgExpression(*it); 

      if(isWholeArrayExpression(var, depth)){
	return true;
      }
    }
  return false; 
}


bool WholeArrayExp::isWholeArrayExpression(SgExpression* exp, int* depth)
			
{
  //Didem: I have seen three types of whole array statements so far
  //1. This one appears as PntrArrRefExp with null lower and upper bound (SgSubscriptExpression)  A(:,:)= B(:,:,m)
  //So we are going to search for all expressions and decide if they fall into case 1 or 2
  //2. A = B + C pure whole statements, they are simply represented as a varrefexp
  //and act very similar to scalars.
  //3. This is a special case of 2. the array might be a member of a struct, or class, so getting its type is a bit tricky

  //Case 1: A(:,:,:) = B(:,:,:,m) ...
  //Case 2: A = B + C  
  //Case 3: mystruct % x = A

  //case 1 : 
  
  if(isSgPntrArrRefExp(exp))
    {
      SgPntrArrRefExp* ptrExp = isSgPntrArrRefExp(exp);
      SgExpression* varExp  = ptrExp->get_lhs_operand(); 
      SgExpression* subExp = ptrExp->get_rhs_operand();
  
      if(isSgVarRefExp(varExp) && isSgExprListExp(subExp))
	{
	  //this should be the array name 
	  SgInitializedName* name = convertRefToInitializedName(varExp);
	  ROSE_ASSERT(name); 
	  
	  SgType* type = name->get_type();
	  
	  if(isSgArrayType(type) || getTypeInNamedType(name, exp))
	    {
	      // SubscriptExpression means this (:), we need to have at least one of them
	      Rose_STL_Container<SgNode*> subList = NodeQuery::querySubTree(subExp, V_SgSubscriptExpression);

	      if(subList.size() > 0 ){

		if (depth != NULL)
		  *depth = subList.size();

		return true;
	      }
	    } //if it is array type or name array type 
	}
    }// end of isSgPntrArrRefExp
  //case 2 and 3
  else if(isSgVarRefExp(exp))
    {
      //check if the variable reference is part of an index expression 
      SgNode* node = exp-> get_parent();

      if(node!= NULL && isSgPntrArrRefExp(node)){
	return false; 
      }
      SgVarRefExp* varExp = isSgVarRefExp(exp);
      ROSE_ASSERT(varExp);
  
      SgInitializedName* name = convertRefToInitializedName(varExp);
      ROSE_ASSERT(name); 

      SgType* type = name->get_type();
      
      if ( !isSgArrayType(type) ) // check if the array is part of a named type struct, class or union
        type = getTypeInNamedType(name, exp); //including class, union or struct

      if(type != NULL &&  isSgArrayType(type)) 
	{
	  SgExprListExp* expList = isSgArrayType(type)->get_dim_info();
	  SgExpressionPtrList& list = expList->get_expressions();	

	  if (depth != NULL)
	    *depth = list.size() ; 

	  return true; 
	}
    }
  return false; 
}

// eof


