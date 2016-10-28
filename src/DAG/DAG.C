#include <rose.h>
#include "DAG.h"

using namespace std;
using namespace SageInterface;


#ifdef EXASAT_DEBUG
void localDagTree::printDefsAtNode(SgNode* current){
  stringstream defUse;
  defUse << "Op:  "<< current->unparseToString().c_str();
  defUse << "   Def  ";
  typedef std::vector<SgInitializedName*> VarName;
  typedef boost::shared_ptr<ReachingDef> ReachingDefPtr;
  typedef std::map<VarName, ReachingDefPtr> NodeReachingDefTable;
  foreach(const NodeReachingDefTable::value_type& varDefPair, ssa->getDefsAtNode(current))
  {
    defUse << " [" << ssa->varnameToString(varDefPair.first) << "]: ";
    defUse << varDefPair.second->getRenamingNumber();
  }
  printf("%s", defUse.str().c_str());
}

void localDagTree::printUses(SgNode* current){
  stringstream defUse;
  defUse << "   Use  ";
  typedef std::vector<SgInitializedName*> VarName;
  typedef boost::shared_ptr<ReachingDef> ReachingDefPtr;
  typedef std::map<VarName, ReachingDefPtr> NodeReachingDefTable;
  foreach(const NodeReachingDefTable::value_type& varDefPair, ssa->getUsesAtNode(current))
  {
    defUse << " [" << ssa->varnameToString(varDefPair.first) << "]: ";
    defUse << varDefPair.second->getRenamingNumber();
  }
  printf("%s\n", defUse.str().c_str());
}
#endif

void localDagTree::findDependencyAndConnect(DagNode* dNode, vector<DagNode*> &accepted){
  bool found=false;
  foreach(const NodeReachingDefTable::value_type& varDefPair, ssa->getUsesAtNode(dNode->_operator._astNode))
  {
    vector<DagNode*>::iterator acceptedNode= accepted.begin();
    while(acceptedNode!=accepted.end()){
      foreach(const NodeReachingDefTable::value_type& varDefPair1, ssa->getDefsAtNode((*acceptedNode)->_operator._astNode))
      {
        if(ssa->varnameToString(varDefPair.first)==ssa->varnameToString(varDefPair1.first) && varDefPair.second->getRenamingNumber() == varDefPair1.second->getRenamingNumber())                     {
          if(count(dNode->_dependency.begin(), dNode->_dependency.end(), *acceptedNode)==0)
	  {
            dNode->_dependency.push_back(*acceptedNode);
            (*acceptedNode)->_children.push_back(dNode);
            (*acceptedNode)->_operands.push_back(ssa->varnameToString(varDefPair.first));
          }
          found=true;
	  #ifdef EXASAT_DEBUG
          printf("%s depends on %s\n", isSgBinaryOp(dNode->_operator._astNode)->unparseToString().c_str(), isSgBinaryOp((*acceptedNode)->_operator._astNode)->unparseToString().c_str());
          #endif
          //break;
        }
      }
      acceptedNode++;
    }
  }
  if(!found){
    dNode->_dependency.push_back(_head);
    _head->_children.push_back(dNode);
    _head->_operands.push_back("");
  }
}

void localDagTree::findDependencyAndConnect(DagNode* dNode, SgNode* branch, vector<DagNode*> &accepted){
  foreach(const NodeReachingDefTable::value_type& varDefPair, ssa->getUsesAtNode(branch))
  {
    bool found=false;
    vector<DagNode*>::iterator acceptedNode= accepted.begin();
    while(acceptedNode!=accepted.end()){
      foreach(const NodeReachingDefTable::value_type& varDefPair1, ssa->getDefsAtNode((*acceptedNode)->_operator._astNode))
      {
        if(ssa->varnameToString(varDefPair.first)==ssa->varnameToString(varDefPair1.first) && varDefPair.second->getRenamingNumber() == varDefPair1.second->getRenamingNumber())                     {
          if(count(dNode->_dependency.begin(), dNode->_dependency.end(), *acceptedNode)==0)
          {
            dNode->_dependency.push_back(*acceptedNode);
            (*acceptedNode)->_children.push_back(dNode);
            (*acceptedNode)->_operands.push_back(ssa->varnameToString(varDefPair.first));
          }
          found=true;
          #ifdef EXASAT_DEBUG
          printf("%s depends on %s\n", isSgBinaryOp(dNode->_operator._astNode)->unparseToString().c_str(), isSgBinaryOp((*acceptedNode)->_operator._astNode)->unparseToString().c_str());
          #endif
        }
      }
      acceptedNode++;
    }
  }
}

bool localDagTree::mayAliasFortran(SgNode* node1, SgNode* node2){
  return true;
}

bool isNodeOfInterest(SgNode* node){
  //return (isSgAssignOp(node) || isSgAddOp(node) || isSgSubtractOp(node) || isSgMultiplyOp(node) || isSgDivideOp(node) || isSgPntrArrRefExp(node)) || isSgGreaterThanOp(node) || isSgLessThanOp(node) || isSgEqualityOp(node);
  return true;
}

//we should generate a dag for each basic block
//If a basic block contains loops, treat each loop as a statement
void localDagTree::build(){
  typedef std::vector<SgInitializedName*> VarName;
  typedef boost::shared_ptr<ReachingDef> ReachingDefPtr;
  typedef std::map<VarName, ReachingDefPtr> NodeReachingDefTable;
  map<SgNode*, DagNode*> dagMap;
  ROSE_ASSERT(_bb);
  _head = new DagNode; 
  _count ++;
  _head->_id=++_currentID;
  _head->_isHead=true;
  _head->_operator._type=null_op;
  _head->_operator._astNode= _bb;
  _head->expType= NULL;
  dagMap[_bb] = _head;
  vector<SgBinaryOp*> ops= querySubTree<SgBinaryOp>(_bb, V_SgBinaryOp); 
  vector<DagNode*> accepted;
  vector<vector<SgBinaryOp*>::iterator> visited;
  stack<vector<SgBinaryOp*>::iterator> traverse;
  vector<SgBinaryOp*>::iterator current= ops.begin();
  traverse.push(current);
  while (!traverse.empty()){
    current = traverse.top();
	  #ifdef EXASAT_DEBUG
	  printDefsAtNode(*current);
          printUses(*current);
	  #endif
    SgNode* lhs= isSgBinaryOp(*current)->get_lhs_operand ();
    SgNode* rhs= isSgBinaryOp(*current)->get_rhs_operand ();
    if (count(visited.begin(), visited.end(), current) == 0)
    {//haven't visited before
      bool backtrack=false; 
      #ifdef EXASAT_DEBUG
       printf("BINARY OP: Type %s %s LHS %s, RHS %s\n", (*current)->get_type()->unparseToString().c_str(), (*current)->unparseToString().c_str(), (*current)->get_lhs_operand ()->unparseToString().c_str(), (*current)->get_rhs_operand ()->unparseToString().c_str());
      #endif
      if (isNodeOfInterest(*current))
      {
        if(querySubTree<SgBinaryOp>(rhs, V_SgBinaryOp).size()==0 && querySubTree<SgBinaryOp>(lhs, V_SgBinaryOp).size()==0){//terminated node
          DagNode* dNode= new DagNode; 
	  _count++;
          dNode->_isHead=false;
          if(isSgAddOp(*current))dNode->_operator._type=add_op;
          else if(isSgSubtractOp(*current))dNode->_operator._type=subtract_op;
               else if(isSgMultiplyOp(*current))dNode->_operator._type=multiply_op;
          	    else if(isSgDivideOp(*current))dNode->_operator._type=divide_op;
	                 else if(isSgPntrArrRefExp(*current))dNode->_operator._type=arrayRef_op;
          dNode->_operator._astNode= *current;
          dNode->_id=++_currentID;
          dNode->expType= (*current)->get_type();
          accepted.push_back(dNode);
          dagMap[*current] = dNode;
          traverse.pop();
          backtrack=true;
          if(ssa->getUsesAtNode(*current).size()==0)//the execution of this node does not depend on any other nodes
          {
            dNode->_dependency.push_back(_head); 
	    _head->_children.push_back(dNode);
	    _head->_operands.push_back("");
	  }else{//look for dependent nodes
	    findDependencyAndConnect(dNode, accepted);
          }
        }else{
	  if(querySubTree<SgBinaryOp>(lhs, V_SgBinaryOp).size()>0){
            vector<SgBinaryOp*>::iterator it= current;
            it++;
            while(it!=ops.end() && *it != *(querySubTree<SgBinaryOp>(lhs, V_SgBinaryOp).begin())){
              it++;
 	    }
  	    if(it!= ops.end())traverse.push(it);
          }
	  if(querySubTree<SgBinaryOp>(rhs, V_SgBinaryOp).size()>0){
            vector<SgBinaryOp*>::iterator it= current;
            it++;
            while(it!=ops.end() && *it != *(querySubTree<SgBinaryOp>(rhs, V_SgBinaryOp).begin())){
              it++;
 	    }
  	    if(it!= ops.end())traverse.push(it);
          }
        }
      }
      visited.push_back(current);
      if((/*!backtrack ||*/ traverse.size()==0) && (current+1)!=ops.end()){//next node will be children of the current node (if any) otherwise next statement
	traverse.push(++current);
      }//no else for back tracking
    }else{//have visited before
      #ifdef EXASAT_DEBUG
      printf("VISITED BINARY OP: %s LHS %s, RHS %s\n", (*current)->unparseToString().c_str(), (*current)->get_lhs_operand ()->unparseToString().c_str(), (*current)->get_rhs_operand ()->unparseToString().c_str());
      #endif
      DagNode* dNode= new DagNode;
      _count++;
      dNode->_isHead=false;
      //dNode->_operator._type=assign_op;
      dNode->_operator._astNode= *current;
      dNode->_id=++_currentID;
      dNode->expType= (*current)->get_type();
      accepted.push_back(dNode);
      dagMap[*current] = dNode;
      if(querySubTree<SgBinaryOp>(lhs, V_SgBinaryOp).size()>0){
        dNode->_dependency.push_back(dagMap[lhs]); 
	dagMap[lhs]->_children.push_back(dNode);
	dagMap[lhs]->_operands.push_back("_subExp");
	#ifdef EXASAT_DEBUG
        printf("%s depends on %s\n", isSgBinaryOp(dNode->_operator._astNode)->unparseToString().c_str(), lhs->unparseToString().c_str());
        #endif
      }else{//the lhs may still contain a varRef
        findDependencyAndConnect(dNode, lhs, accepted);
      }
      if(querySubTree<SgBinaryOp>(rhs, V_SgBinaryOp).size()>0){
        dNode->_dependency.push_back(dagMap[rhs]); 
	dagMap[rhs]->_children.push_back(dNode);
	dagMap[rhs]->_operands.push_back("_subExp");
	#ifdef EXASAT_DEBUG
        printf("%s depends on %s\n", isSgBinaryOp(dNode->_operator._astNode)->unparseToString().c_str(), rhs->unparseToString().c_str());
        #endif
      }else{//rhs may still contain a varRef
        findDependencyAndConnect(dNode, rhs, accepted);
      }
      //at this time, we have done with children of the current node. So we pop it and look for the next non-children node
      traverse.pop();
      if(traverse.size()==0){
        current++;
        while(current!=ops.end() && count(visited.begin(), visited.end(), current) > 0){
          current++;
	}
	if(current!= ops.end())traverse.push(current);
      }
    }
  }//end while
}


bool basicBlockOfInterest(SgBasicBlock* bb){
#if 0
  //here we are interested in basic blocks that don't contain any loops
  vector<SgNode*> stmtList = NodeQuery::querySubTree(bb, V_SgStatement);
  for (vector<SgNode*>::iterator istmt= stmtList.begin(); istmt!= stmtList.end(); istmt++){
    if(isSgForStatement(*istmt)) return false;
  }
#endif
  //we are interested in all basic blocks, including those that contain loops
  return true;
}

void SSA_DAG::generateDAGs(map<string,SgFunctionDefinition*> reachableFuncDefList, map<string,SgBasicBlock*> basicBlockList){ 
  for(map<string,SgBasicBlock*>::iterator bbItr= basicBlockList.begin(); bbItr != basicBlockList.end(); bbItr++)
  {
    SgBasicBlock* bb= (*bbItr).second;
    ROSE_ASSERT(bb);
    if(basicBlockOfInterest(bb)){//we just care about basic blocks that don't contain loops
      //if(isSgForStatement(bb->get_parent())){
        localDagTree *ldag= new localDagTree(bb, dagIdBase);
        ldag->setSSA(ssa);
        ldag->build();
        localdagMap[bb]= ldag->getHead();
        dagIdBase += ldag->size();
      //}
    }
  }
#if 0
  for(int f=0; f<project->numberOfFiles(); f++){
    Rose_STL_Container<SgNode*> bbList = NodeQuery::querySubTree(project->get_fileList()[f], V_SgBasicBlock);
    for(Rose_STL_Container<SgNode*>::iterator it= bbList.begin(); it!= bbList.end(); it++){
      SgBasicBlock* bb= isSgBasicBlock(*it);
      ROSE_ASSERT(bb);
      if(basicBlockOfInterest(bb)){//we just care about basic blocks that don't contain loops
        if(isSgForStatement(bb->get_parent())){
          localDagTree *ldag= new localDagTree(bb, dagIdBase);
          ldag->setSSA(ssa);
          ldag->build();
          localdagMap[bb]= ldag->getHead();
          dagIdBase += ldag->size();
        }
      }
    }
  }
#endif
}


void SSA_DAG::recurDAGPrint(DagNode* dNode, ofstream & outFile, vector<DagNode*> &visited){
if(count(visited.begin(), visited.end(), dNode)==0){
  visited.push_back(dNode);
  //Print this node
  outFile << dNode->_id;
  outFile <<  " [label=\"";
  if(isSgBasicBlock(dNode->_operator._astNode))
    outFile <<  "BasicBlock";
  else outFile << dNode->_operator._astNode->unparseToString() << "  Type:" << dNode->expType->unparseToString();
  outFile <<  "\"]\n";
  //Now print the out edges
  vector<DagNode*> children= dNode->_children;
  for(vector<DagNode*>::iterator childrenItr= children.begin(); childrenItr != children.end(); childrenItr++)
  {
     outFile << dNode->_id;
     outFile << " -> ";
     outFile << (*childrenItr)->_id;
     outFile <<  " [label=\"";
     outFile << dNode->_operands[childrenItr-children.begin()];
     outFile <<  "\"]\n";
     recurDAGPrint(*childrenItr, outFile, visited);
  }
}
}

void SSA_DAG::dagsToDOT(map<string,SgFunctionDefinition*> reachableFuncDefList, map<string,SgBasicBlock*> basicBlockList){
  //fileVec files = SageInterface::querySubTree<SgSourceFile > (project, V_SgSourceFile);
  ofstream outFile("DAG_OUTPUT.dot");
  outFile << "digraph DAGGraph {\n";
  vector<DagNode*> visited;
  for(map<string,SgFunctionDefinition*>::iterator funcListItr= reachableFuncDefList.begin(); funcListItr != reachableFuncDefList.end(); funcListItr++)
  {
    SgFunctionDefinition* funcDef= (*funcListItr).second;
    ROSE_ASSERT(funcDef);
    vector<SgNode*> bbList= SageInterface::querySubTree<SgNode>(funcDef, V_SgBasicBlock);
    for(vector<SgNode*>::iterator ibb= bbList.begin(); ibb!=bbList.end(); ibb++){
      SgBasicBlock* bb= isSgBasicBlock(*ibb);
      if(bb)
        if(localdagMap[bb]){
          recurDAGPrint(localdagMap[bb], outFile, visited);
        }
    }
  }
  for(map<string,SgBasicBlock*>::iterator bbItr= basicBlockList.begin(); bbItr != basicBlockList.end(); bbItr++)
  {
    SgBasicBlock* bb= (*bbItr).second;
    ROSE_ASSERT(bb);
    if(localdagMap[bb]){
      recurDAGPrint(localdagMap[bb], outFile, visited);
    }
  }
  outFile << "}\n";
}


void SSA_DAG::generateSSA(){
  ssa= new SSA_extension(project);
  ssa->run(false, true);
}

void SSA_DAG::ssaToDOT(){
    fileVec files = SageInterface::querySubTree<SgSourceFile > (project, V_SgSourceFile);

    foreach(fileVec::value_type& file, files)
    {
        ofstream outFile((StringUtility::stripPathFromFileName(file->getFileName())
                + "_SSA_DOT").c_str());
        outFile << "digraph SSAGraph {\n";
        funcDefVec funcs = SageInterface::querySubTree<SgFunctionDefinition > (file, V_SgFunctionDefinition);
        foreach(funcDefVec::value_type& func, funcs)
        {
          vector<cfgNode> visited;
          stack<cfgNode> traverse;
          cfgNode current = cfgNode(func->cfgForBeginning());
          traverse.push(current);
          while (!traverse.empty())
          {
            current = traverse.top();
            if (count(visited.begin(), visited.end(), current) == 0)
            {
                string id = current.id();
                string nodeColor = "black";

                bool uniqueName = current.getNode()->attributeExists(UniqueNameTraversal::varKeyTag);

                if (isSgStatement(current.getNode()))
                    nodeColor = "blue";
                else if (isSgExpression(current.getNode()))
                    nodeColor = "green";
                else if (isSgInitializedName(current.getNode()))
                    nodeColor = "red";

                string name = "";
                if (uniqueName)
                {
                    VarUniqueName *attr = StaticSingleAssignment::getUniqueName(current.getNode());
                    ROSE_ASSERT(attr);
                    name = attr->getNameString();
                }
                stringstream defUse;
    	        typedef std::vector<SgInitializedName*> VarName;
    	        typedef boost::shared_ptr<ReachingDef> ReachingDefPtr;
    	        typedef std::map<VarName, ReachingDefPtr> NodeReachingDefTable;
                foreach(const NodeReachingDefTable::value_type& varDefPair, ssa->getOutgoingDefsAtNode(current.getNode()))
                {
                    defUse << "Def [" << ssa->varnameToString(varDefPair.first) << "]: ";
                    defUse << varDefPair.second->getRenamingNumber() << " - "
                            << (varDefPair.second->isPhiFunction() ? "Phi" : "Concrete") << "\\n";
                }
    	        //typedef std::map<VarNameArray, ReachingDefPtr_ext> NodeReachingDefTableArray;
                foreach(const SSA_extension::NodeReachingDefTable_array::value_type& varDefPair, ssa->getOutgoingDefsAtNode_array(current.getNode()))
                {
                    defUse << "ARRAY Def [" << varDefPair.first->arrayName[0]->unparseToString().c_str() << " ( ";
		    for(int d=0; d<varDefPair.first->rank[0]; d++){
                      defUse << varDefPair.first->maskVector[0]->def.lo.x[d].s->unparseToString().c_str() << " : " << varDefPair.first->maskVector[0]->def.hi.x[d].s->unparseToString().c_str(); 
		      if(d+1 < varDefPair.first->rank[0]) defUse << " , ";
		    }
                    defUse << " ) " <<endl;
//                            << (varDefPair.second->isPhiArrFunction() ? "PhiArr" : "Concrete") << "\\n";
                }
                string defUseStr = defUse.str().substr(0, defUse.str().size() - 2);
                string label = escapeString(current.getNode()->class_name());
                if (isSgFunctionDefinition(current.getNode()))
                    label += ":" + escapeString(isSgFunctionDefinition(current.getNode())->get_declaration()->get_name());

                //Print this node
                outFile << id << " [label=\"<" << label << ">:" << current.getNode()
                        //Now we add the unique name information
                        << ((name != "") ? "\\n" : "") << name
                        << ((defUseStr != "") ? "\\n" : "") << defUseStr
                        << "\", color=\"" << nodeColor << "\", style=\""
                        << (current.isInteresting() ? "solid" : "dotted") << "\"];\n";

                //Now print the out edges
                vector<cfgEdge> outEdges = current.outEdges();

                foreach(vector<cfgEdge>::value_type& edge, outEdges)
                {
                    outFile << edge.source().id() << " -> " << edge.target().id()
                            << " [label=\"" << escapeString(edge.toString())
                            << "\"];\n";
                }
            }//end if
            visited.push_back(current);

            vector<cfgEdge> outEdges = current.outEdges();
            foreach(vector<cfgEdge>::value_type& edge, outEdges)
            {
                //If we haven't seen the target of this node yet, process the node
                if (count(visited.begin(), visited.end(), edge.target()) == 0)
                {
                    traverse.push(edge.target());
                    //break;
                }
            }

            //If there are no new out edges to explore
            if (traverse.top() == current)
            {
                vector<cfgEdge> inEdges = current.inEdges();
                foreach(vector<cfgEdge>::value_type& edge, inEdges)
                {
                    //If we haven't seen the target of this node yet, process the node
                    if (count(visited.begin(), visited.end(), edge.target()) == 0)
                    {
                        traverse.push(edge.source());
                        //break;
                    }
                }
            }

            //No out or in edges left to print, pop this node
            if (traverse.top() == current)
            {
                traverse.pop();
            }
          }//end while
        }//end for each
        outFile << "}";
  }
}


