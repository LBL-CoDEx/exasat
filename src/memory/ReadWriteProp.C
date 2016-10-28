/*
  Didem Unat
 */

#include <rose.h>
#include "sage3basic.h"
#include <iostream>
#include <list>
#include <string>
//user defined classes
#include "ReadWriteProp.h"
#include "../utils/Utils.h"
#include "../utils/ExaSatOptions.h"
#include "../utils/ArrayUtils.h"
#include "../utils/OutputFormatting.h"
#include "./FuncSideEffect.h"

#include <algorithm> //for set_difference

//ROSE classes
#include "ArrayAnnot.h"
#include "ArrayInterface.h"
#include "LoopTransformInterface.h"
#include "AstInterface_ROSE.h"
#include "DepInfoAnal.h"

using namespace std;
using namespace SageInterface;
using namespace ArrayUtils;
using namespace OutputFormat;

void ReadWriteProp::getReadOnlyVars(set<NameComp_t>& readOV)
{
  for(std::set<NameComp_t>::iterator it = readOnlyVars.begin(); 
      it!= readOnlyVars.end(); it++)
    {
      SgInitializedName* init_name = (*it).first;
      string component = (*it).second;
      SgType* type = init_name ->get_type();

      //According to Fortran, all arrays have 1 dimension
      if(isArrayOrPointerType(type))// && getDimension(init_name, type) >= 3)
	{
	  readOV.insert(NameComp_t(init_name,component));
	}
    }
}

void ReadWriteProp::getWriteOnlyVars(set<NameComp_t>& writeOV)
{
  for(std::set<NameComp_t>::iterator it = writeOnlyVars.begin(); 
      it!= writeOnlyVars.end(); it++)
    {
      SgInitializedName* init_name = (*it).first;
      string component = ""; 
      component = (*it).second;
      SgType* type = init_name ->get_type();

      //According to Fortran, all arrays have 1 dimension
      if(isArrayOrPointerType(type))// && getDimension(init_name, type) >= 3)
	{
	  writeOV.insert(NameComp_t(init_name,component));
	}
    }
}

void ReadWriteProp::getBothReadWrittenVars(set<NameComp_t>& writereadOV)
{
  for(std::set<NameComp_t>::iterator it = bothReadWrittenVars.begin(); 
      it!= bothReadWrittenVars.end(); it++)
    {
      SgInitializedName* init_name = (*it).first;
      string component = ""; 
      component = (*it).second;
      SgType* type = init_name ->get_type();

      //According to Fortran, all arrays have 1 dimension
      if(isArrayOrPointerType(type))// && getDimension(init_name, type) >= 3)
	{
	  writereadOV.insert(NameComp_t(init_name,component));
	}
    }
}



bool ReadWriteProp::collectReadWriteRefs(SgStatement* stmt, 
					 std::vector<SgNode*>& readRefs, 
					 std::vector<SgNode*>& writeRefs)
{   
  // The type cannot be SgExpression since variable declarations have SgInitializedName as the reference, not SgVarRefExp.                                                                  
  ROSE_ASSERT(stmt !=NULL);

  // convert a request for a defining function declaration to its function body                                                                                                               
  SgFunctionDeclaration* funcDecl = isSgFunctionDeclaration(stmt);
  if (funcDecl != NULL)
    {
      funcDecl= isSgFunctionDeclaration(funcDecl->get_definingDeclaration ());
      if (funcDecl == NULL)
	{
	  cerr<<"In collectReadWriteRefs(): cannot proceed without a function body!"<<endl;
	}
      stmt = funcDecl->get_definition()->get_body();
    }

  // get function level information                                                                                                                                                           
  SgFunctionDefinition* funcDef = SageInterface::getEnclosingFunctionDefinition(stmt);
  ROSE_ASSERT(funcDef != NULL);
  SgBasicBlock* funcBody = funcDef->get_body();
  ROSE_ASSERT(funcBody!= NULL);


  // why? Didem: this was in ROSE and don't understand why they need this.
  // prepare Loop transformation environment                                                                                                                                                  
  AstInterfaceImpl faImpl(funcBody);
  AstInterface fa(&faImpl);

  ArrayAnnotation* annot = ArrayAnnotation::get_inst();
  ArrayInterface array_interface (*annot);

  //The following line caused termination in ROSE, I commented out.
  //it performs alias analysis 
  //array_interface.initialize(fa, AstNodePtrImpl(funcDef));

  array_interface.observe(fa);

  FunctionSideEffectAnalysis* func_sideEffectInfo = new FunctionSideEffectAnalysis(); 
  LoopTransformInterface::set_arrayInfo(&array_interface);
  LoopTransformInterface::set_astInterface(fa);
  LoopTransformInterface::set_sideEffectInfo(func_sideEffectInfo);

  // variables to store results     
  DoublyLinkedListWrap<AstNodePtr> rRef1, wRef1;
  CollectDoublyLinkedList<AstNodePtr> crRef1(rRef1),cwRef1(wRef1);
  AstNodePtr s1 = AstNodePtrImpl(stmt);

  // Actual side effect analysis                                                                            
  if (!AnalyzeStmtRefs(fa, s1, cwRef1, crRef1))
    {
      //cerr<<"error in side effect analysis!"<<endl;
      //  return false;
    }
  delete func_sideEffectInfo; 

  // transfer results into STL containers.                                                                                               
  for (DoublyLinkedEntryWrap<AstNodePtr>* p = rRef1.First(); p != 0; )
    {
      DoublyLinkedEntryWrap<AstNodePtr>* p1 = p;
      p = rRef1.Next(p);
      AstNodePtr cur = p1->GetEntry();
      SgNode* sgRef = AstNodePtrImpl(cur).get_ptr();
      ROSE_ASSERT(sgRef != NULL);
      readRefs.push_back(sgRef);
    }

  for (DoublyLinkedEntryWrap<AstNodePtr>* p = wRef1.First(); p != 0; )
    {
      DoublyLinkedEntryWrap<AstNodePtr>* p1 = p;
      p = wRef1.Next(p);
      AstNodePtr cur = p1->GetEntry();
      SgNode* sgRef = AstNodePtrImpl(cur).get_ptr();
      ROSE_ASSERT(sgRef != NULL);
      writeRefs.push_back(sgRef);
    }
  return true;
}

string getComponentAndDim(SgNode* current, int& dim)
{
  string component="";
  SgPntrArrRefExp* exp = isSgPntrArrRefExp(current);

  if(exp != NULL)
    {
      //get which component we are working on 
      SgExpression* array_name = NULL;
      std::vector<SgExpression*> subscripts; 
      
      if(ArrayUtils::isArrayReference(isSgExpression(exp), &array_name, &subscripts))
	{
	  dim = subscripts.size(); 

	  if(dim > 3 ){
	    //last one is the component, not always the case, parameterize this	   
	    SgExpression* componentExp = subscripts[dim-1];
	    //we keep the string, because it can be integer or varref
	    component = componentExp->unparseToString();
	  }
	}    
    }
  else //might be whole array reference expression (F90 only)
    {
      /*
      if(ArrayUtils::isWholeArrayReference(isSgExpression(current), NULL, &dim))
	{
	  if(dim > 3)
	    {
	      //TODO: Need to handle the components, do I? 
	    }
	}
      */
    }
  return component; 
}

bool ReadWriteProp::collectReadWriteVariables(SgStatement* stmt, 
					      set<NameComp_t>& readVars, 
					      set<NameComp_t>& writeVars)
{
  ROSE_ASSERT(stmt != NULL);
  vector <SgNode* > readRefs, writeRefs;

  //ROSE version breaks, I have my own version of collect read write refs 
  //if (!SageInterface::collectReadWriteRefs(stmt, readRefs, writeRefs))
  if (!ReadWriteProp::collectReadWriteRefs(stmt, readRefs, writeRefs))
    {
      return false;
    }
  // process read references                                                                                                                                                                  
  vector<SgNode*>::iterator iter = readRefs.begin();
  for (; iter!=readRefs.end();iter++)
    {
      SgNode* current = *iter;
      int dim=0; 
      if(!isSgFunctionRefExp(current)){
        SgInitializedName* name;
        if(isSgPointerDerefExp(current))
        {
	  SgVarRefExp* pointerVarRef;
          Rose_STL_Container<SgNode*> varList = NodeQuery::querySubTree(current, V_SgVarRefExp);
	  for(Rose_STL_Container<SgNode*>::iterator it= varList.begin(); it!= varList.end(); it++){
	    SgType* type = isSgVarRefExp(*it)->get_type();
	    if(isSgPointerType(type)){
	      pointerVarRef= isSgVarRefExp(*it);
	      break;
	    }
          }
          name = convertRefToInitializedName(pointerVarRef);
        }else name = convertRefToInitializedName(current);

        string component = getComponentAndDim(current, dim);

        // We use std::set to ensure uniqueness now                                                                    
        if(dim >= 1)
	  readVars.insert(NameComp_t(name,component));
      }//Tan: convertRefToInitializedName does not support function pointer  
      //
    }
  // process write references                                                                                                                                                                 
  vector<SgNode*>::iterator iterw = writeRefs.begin();
  for (; iterw!=writeRefs.end();iterw++)
    {
      SgNode* current = *iterw;
      int dim=0;
      if(!isSgFunctionRefExp(current)){
        string component = getComponentAndDim(current, dim);
        SgInitializedName* name;
        if(isSgPointerDerefExp(current))
        {
          SgVarRefExp* pointerVarRef;
          Rose_STL_Container<SgNode*> varList = NodeQuery::querySubTree(current, V_SgVarRefExp);
          for(Rose_STL_Container<SgNode*>::iterator it= varList.begin(); it!= varList.end(); it++){
            SgType* type = isSgVarRefExp(*it)->get_type();
            if(isSgPointerType(type)){
              pointerVarRef= isSgVarRefExp(*it);
              break;
            }
          }
          name = convertRefToInitializedName(pointerVarRef);
        }else name = convertRefToInitializedName(current);
        // We use std::set to ensure uniqueness now                    
        if(dim >= 1)
	  writeVars.insert(NameComp_t(name,component));
      }
    }//Tan: convertRefToInitializedName does not support function pointer  
  return true;
}

void printSet (const std::set<NameComp_t> vars, string category, string tag) 
{
  
  cout << category ; 
  
  for(std::set<NameComp_t>::const_iterator it = vars.begin(); it!= vars.end(); it++)
    {
      SgInitializedName* init_name = (*it).first;
      string component = ""; 
      component = (*it).second;
      SgType* type = init_name ->get_type();

      //According to Fortran, all arrays have 1 dimension
      if(isArrayOrPointerType(type))// && getDimension(init_name, type) >= 3)
	{
	  string name_str = init_name -> get_name().str();

	  if(!component.empty())
	    cout << name_str << "(" << component << ")"<< ", "; 
	  else 
	    cout << name_str << ", ";
	  
	}
    }
  cout << endl;
}

void ReadWriteProp::run(SgNode* node)
{
  //the following call will generate an error for Fortran routines 
  //SageInterface::collectReadWriteVariables(stmt, readVars, writeVars);
  
  SgStatement* stmt = isSgStatement(node);
  ROSE_ASSERT(stmt);
  bool result = ReadWriteProp::collectReadWriteVariables(isSgStatement(node), readVars, writeVars);
  ROSE_ASSERT(result);

  //get W U R 
  set_union(writeVars.begin(), writeVars.end(),
	    readVars.begin(), readVars.end(),
	    std::inserter( allVars, allVars.begin()));


  //get read only vars R-W
  set_difference(readVars.begin(), readVars.end(),
		 writeVars.begin(), writeVars.end(),
		 std::inserter(readOnlyVars, readOnlyVars.begin()));


  //get write only vars W-R
  set_difference(writeVars.begin(), writeVars.end(),
		 readVars.begin(), readVars.end(),
		 std::inserter( writeOnlyVars, writeOnlyVars.begin()));

  //get both written and read vars W intersection R
  set_intersection(writeVars.begin(), writeVars.end(),
		   readVars.begin(), readVars.end(),
		   std::inserter( bothReadWrittenVars, bothReadWrittenVars.begin()));

  //printSet(allVars, "     All Vars              : ");
  if(!ExaSatOptions::xmlOutput)
    {
    printSet(readOnlyVars, "     Read Only Vars        : ", "readonly");
    
    printSet(writeOnlyVars, "     Write Only Vars       : ", "writeonly");
    
    printSet(bothReadWrittenVars,"     Both Written/Read Vars: ", "readwrite" );
  }
}

//eof
