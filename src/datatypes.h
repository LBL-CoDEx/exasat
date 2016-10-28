#ifndef DATATYPES_H

#define DATATYPES_H

#include <map>
#include <vector>
#include <string>
#include "./handler/NodeHandler.h"


#define MAX_DIM 7 //max rank in Fortran arrays
using namespace std;
using namespace SageBuilder;
using namespace SageInterface;

typedef struct{
  long int lineno;
  SgInitializedName* indexVar;
  SgExpression* lower; 
  SgExpression* upper;
  SgExpression* stride;
}LoopAttributes_t;

typedef struct{
  long int lineno;
  SgStatement* conditional;
}IfAttributes_t;

typedef struct{
  SgBasicBlock* bb;
}SpecialBBAttributes_t;

typedef struct{
  int addOps; 
  int mulOps; 
  int divOps; 
  int specialOps; 
}Flops_t;

typedef struct{
  //we keep the statement for cases when an var appears both sides of an assignment
  //in those cases, write appears first in AST and then read
  //if the stmt is encountered before, then read-first overwrites the write-first 
  SgStatement* stmt;
  bool firstRead;
}FirstAccess_t;

typedef struct{
  //-1 if it doesn't depend on a loop
  //innermost loop is zero
  std::vector<string> dependentLoopVar; 
  std::vector<string> offset;
  int readCount;
  int writeCount; 
}AccessPattern_t;


typedef struct{
  std::vector <SgVariableDeclaration* > indexVar;
  std::vector <SgExpression* > lower;
  std::vector <SgExpression* > upper; 
  std::vector <SgExpression* > stride;
}RangeInfo_t;

typedef struct{
  long int lineno;
  string func_name;
  string func_rename; 
  string module_name;
  string flops;
  string bytes;
  std::map<SgInitializedName*, string> arg_map;
}FuncCallInfo_t;


//access string is used as a hash key
typedef std::map<string, AccessPattern_t> AccessListPerArray_t; 

//SgInitializedName is the array name 
//string is the component name, conponent can be a constant or a variable name 
typedef std::pair< SgInitializedName*, string> NameComp_t;

//Name component pair maps to the access list 
typedef std::map < NameComp_t, AccessListPerArray_t > ArraysAccessInfo_t ; 
    
//Name component pair maps to first access info 
//There is only one info per name-comp pair, not per access pattern
typedef std::map <NameComp_t, FirstAccess_t> FirstAccessMap_t;


typedef std::map <SgInitializedName*, AccessPattern_t> ScalarRWMap_t;

typedef std::map< NodeHandler*, vector<NodeHandler*> > ParentToChildrenMap_t; 

typedef std::map< SgExpression*, RangeInfo_t > LoopBounds_t;
#endif
