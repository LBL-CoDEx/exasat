#ifndef _DAG_H
#define _DAG_H

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>
#include "ASTtools.hh"
#include "CallGraph.h"
#include "staticSingleAssignment.h"
#include "uniqueNameTraversal.h"
#include "/opt/rose/SSA_ext/src/arraySSA/SSA_extension.h"
#include "/opt/rose/SSA_ext/src/arraySSA/reachingDef_ext.h"

#undef EXASAT_DEBUG
#undef ARRAY_EARLY_TERMINATION

using namespace SageInterface;
#define foreach BOOST_FOREACH
using namespace std;
using namespace boost;
using namespace ssa_private;
using namespace std;

enum operatorType{add_op, subtract_op, multiply_op, divide_op, assign_op, arrayRef_op, null_op};

struct DagNode
{
  struct{
    operatorType _type;
    SgNode* _astNode;
  } _operator;
  vector<DagNode*> _dependency; 
  vector<DagNode*> _children;
  vector<string> _operands;
  SgType* expType;
  int _id;
  bool _isHead;//we need one head node for each basic block
};

//Dag of each basic block
class localDagTree
{
  private:
    DagNode* _head;
    std::list<DagNode*> _leaves;
    SgBasicBlock* _bb; 
    int _currentID;
    int _count;
    SSA_extension *ssa;
    typedef std::vector<SgInitializedName*> VarName;
    typedef boost::shared_ptr<ReachingDef> ReachingDefPtr;
    typedef std::map<VarName, ReachingDefPtr> NodeReachingDefTable;
  public:
    localDagTree(){
      _head=NULL;
      _count=0;
      _currentID=-1;
    }
    localDagTree(SgBasicBlock* b, int dagIdBase){
      _head=NULL;
      _count=0;
      _currentID=dagIdBase;
      _bb=b; 
    }
    DagNode* getHead(){
      return _head;
    }
    std::list<DagNode*> getLeaves(){
      return _leaves;
    }
    int size(){
      return _count;
    }
    bool mayAliasFortran(SgNode*, SgNode*);
#ifdef EXASAT_DEBUG
    void printDefsAtNode(SgNode*);
    void printUses(SgNode*);
#endif
    void setBasicBlock(SgBasicBlock* b){_bb = b;}
    void build();
    void setSSA(SSA_extension* /*StaticSingleAssignment*/ inSSA){ssa = inSSA;}
    void findDependencyAndConnect(DagNode*, vector<DagNode*> &);
    void findDependencyAndConnect(DagNode*, SgNode*, vector<DagNode*> &);
};


class SSA_DAG
{
 private:
  typedef vector<SgSourceFile*> fileVec;
  typedef CFGNode cfgNode;
  typedef CFGEdge cfgEdge;
  typedef vector<SgFunctionDefinition*> funcDefVec;
  SgProject* project;
  SSA_extension* /*StaticSingleAssignment */ssa;
  map<SgBasicBlock*, DagNode*> localdagMap; //each entry of this map is the head node of a DAG associated with a basic block
  void recurDAGPrint(DagNode*, ofstream &, vector<DagNode*>&);
  int dagIdBase;

 public:
  SSA_DAG():dagIdBase(-1){};
  void setBase(int base){dagIdBase =base;}
  void setProject(SgProject* p){project=p;}
  void generateSSA();
  void generateDAGs(map<string,SgFunctionDefinition*> reachableFuncDefList, map<string,SgBasicBlock*> basicBlockList);
  void ssaToDOT();
  void dagsToDOT(map<string,SgFunctionDefinition*>, map<string,SgBasicBlock*>);
};

#endif //_DAG_H
