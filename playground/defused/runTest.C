S/****************************************** 
 * Category: DFA
 * Application to test the Def-Use Analysis
 * modified by Didem to be used in
 * register-global memory optimizations
 *****************************************/
#include "rose.h"
#include "DefUseAnalysis.h"
#include <string>
#include <iostream>

#include "MintArrayInterface.h"
using namespace std;

void testOneFunction( std::string funcParamName, 
		      vector<string> argvList,
		      bool debug, int nrOfNodes, 
		      multimap <string, int> results,
		      multimap <string, int> useresults) {
  if (debug)
    cout << " \n\n------------------------------------------\nrunning ... " << argvList[1] << endl;
  // Build the AST used by ROSE
  SgProject* project = frontend(argvList);
  // Call the Def-Use Analysis
  DefUseAnalysis* defuse = new DefUseAnalysis(project);
  int val = defuse->run(false);//debug);
  
  if (debug)
    std::cout << "Analysis run is : " << (val ?  "failure" : "success" ) << " " << val << std::endl;
  if (val==1) exit(1);
  
  /*
  if (debug==false)
    defuse->dfaToDOT();
  */

  //std::list<SgNode*> vars = NodeQuery::querySubTree(project, V_SgFunctionDefinition); 
  //std::list<SgNode*>::const_iterator i = vars.begin();
  NodeQuerySynthesizedAttributeType vars = NodeQuery::querySubTree(project, V_SgFunctionDefinition); 
  NodeQuerySynthesizedAttributeType::const_iterator i = vars.begin();
  for (; i!=vars.end();++i) {
    SgFunctionDefinition* func = isSgFunctionDefinition(*i);
    std::string name = func->class_name();
    string funcName = func->get_declaration()->get_qualified_name().str();

    /*
    int maxNodes = defuse->getDefSize();
    int nodeNr = defuse->getIntForSgNode(func);
    if (nodeNr == -1)
      continue;
    */
    //cout << " checking function : " << funcName << endl;
    if (funcName!=funcParamName)
      continue;

    
    defuse->start_traversal_of_one_function(func);

    SgBasicBlock* func_body = func->get_body();

    NodeQuerySynthesizedAttributeType vars = NodeQuery::querySubTree(func_body, V_SgVarRefExp);

    NodeQuerySynthesizedAttributeType::const_iterator i = vars.begin();

    for (; i!=vars.end();++i) {

      SgVarRefExp* varRef = isSgVarRefExp(*i);

      SgInitializedName* initName = isSgInitializedName(varRef->get_symbol()->get_declaration());

      std::string name = initName -> get_qualified_name().str();

      
      if(name == "flux" || name == "cons" || name == "Eprev" )
	{// || name == n){
      //Find reaching defition of initname at the control flow node 
      
      vector<SgNode*> vec = defuse-> getDefFor(varRef, initName);
      vector<SgNode*> vec2 = defuse-> getUseFor(varRef, initName);

      ROSE_ASSERT( vec.size() > 0);

      std::cout << "---------------" << endl;

      std::cout << vec.size() << " definition entries for " << varRef-> unparseToString() <<
	" @ line " << varRef->get_file_info()->get_line() << ":" << varRef->get_file_info()->get_col() << endl ;


      for(size_t j= 0 ; j < vec.size(); j++)
	{

	  //cout << vec[j] -> class_name() << " " << vec[j]->unparseToString() << endl;
	  SgStatement* def_stmt = SageInterface::getEnclosingStatement(vec[j]);
	


	  cout<< def_stmt -> unparseToString() << 
	    " @ line " << def_stmt->get_file_info()->get_line() << ":" << def_stmt->get_file_info()->get_col() << endl ;
	}

      std::cout << vec2.size() << " use entries for " << varRef-> unparseToString() <<
	" @ line " << varRef->get_file_info()->get_line() << ":" << varRef->get_file_info()->get_col() << endl ;

      SgScopeStatement* maxScope=NULL ; 
      for(size_t j= 0 ; j < vec2.size(); j++)
	{

	  //cout << vec[j] -> class_name() << " " << vec[j]->unparseToString() << endl;
	  SgStatement* def_stmt = SageInterface::getEnclosingStatement(vec2[j]);

	  SgScopeStatement* temp= def_stmt->get_scope();

	  if(maxScope == NULL)
	    maxScope = temp ; 

	  else if(maxScope > temp ) 
	    maxScope = temp;
	  
	  cout<< def_stmt -> unparseToString() << 
	    " @ line " << def_stmt->get_file_info()->get_line() << ":" << def_stmt->get_file_info()->get_col() << endl ;
	}
      }

    }


    /*

    Rose_STL_Container<SgNode*> sList = NodeQuery::querySubTree(func_body, V_SgVarRefExp);
    
    Rose_STL_Container<SgNode*>::iterator s;
    
    for(s = sList.begin(); s != sList.end(); s++)
      {
	ROSE_ASSERT(*s);
	
	SgVarRefExp* sExp = isSgVarRefExp(*s);
	
	cout << "Mint: " << sExp->unparseToString() << endl;

	cout << AstToString(isSgAstNodePtr(*s)) << endl;

	std::vector <std::pair < SgInitializedName*, SgNode*> > map = defuse->getDefMultiMapFor(*s);
	
	if (map.size()>0) {
	  
	  std::vector < std::pair <SgInitializedName*, SgNode*> >::const_iterator j = map.begin();
	  
	  for (;j!=map.end();++j) {
	    
	    SgInitializedName* in_node = j->first;
	    
	    SgNode* node = j->second;
	    string name= in_node->get_qualified_name().str();
	    if (debug)
	      cout << " DEF checking :  " << name << endl;
	  }

	  std::vector <std::pair < SgInitializedName*, SgNode*> > map2 = defuse->getUseMultiMapFor(*s);

	  if (map2.size()>0) {
	   
	    std::vector < std::pair <SgInitializedName*, SgNode*> >::const_iterator jj = map2.begin();
	   
	    for (;jj!=map2.end();++jj) {
	      
	      SgInitializedName* in_node = jj->first;
	      SgNode* node = jj->second;
	      string name= in_node->get_qualified_name().str();
	      if (debug)
		cout << " USED checking :  " << name << endl;
	    }	  
	  }
	}
	}*/
  }
}

    /*
    if (debug)
      cout << "\n------------------------\nchecking for " << name << " -- " << funcName << " -- " << nodeNr << endl;
    if (maxNodes!=nrOfNodes) {
      if (debug)
	cout << " Error: Test should have " << nrOfNodes << " nodes. found: " << maxNodes << endl;
      abort();
    }
    if (debug)
      cout << " Test has nodes:  " << nrOfNodes <<  endl;
    if (debug)
      cout <<"\nChecking all definitions ... " << endl;
    // check nodes in multimap
  
    std::vector <std::pair < SgInitializedName*, SgNode*> > map = defuse->getDefMultiMapFor(func);
    if (map.size()>0) {
      std::vector < std::pair <SgInitializedName*, SgNode*> >::const_iterator j = map.begin();
      unsigned int hit=0;
      for (;j!=map.end();++j) {
	SgInitializedName* in_node = j->first;
	SgNode* node = j->second;
	string name= in_node->get_qualified_name().str();
	if (debug)
	  cout << " ... checking :  " << name << endl;
	multimap <string, int>::const_iterator k =results.begin();
	for (;k!=results.end();++k) {
	  string resName = k->first;
	  int resNr = k->second;
	  int tableNr = defuse->getIntForSgNode(node);
	  if (name==resName)
	    if (debug)
	      cout << " ... defNr: " << resNr << "  inTable: " << tableNr <<  endl; 
	  if (name==resName && tableNr==resNr) {
	    hit++;
	    if (debug)
	      cout << " Hit " << hit << "/" << results.size() << " - (" << name << "," << resNr << ")" << endl;
	  }
	}

      }
      if (hit!=results.size()) {
	cout << " Error: No hit! ... DFA values of node " << nrOfNodes << " are not correct! " << endl;
	exit(1);
      }
    } else {
      if (results.size()!=0) {
	cout << " Error: Test node " << defuse->getIntForSgNode(func) << " should have a multimap. " << endl;
	exit(1);
	}

  }
  
    if (debug)
      cout <<"\nChecking all uses ... " << endl;
    // check nodes in multimap
    map = defuse->getUseMultiMapFor(func);
    if (map.size()>0) {
      std::vector <std::pair <SgInitializedName*, SgNode*> >::const_iterator j = map.begin();
      size_t hit=0;
      for (;j!=map.end();++j) {
	SgInitializedName* in_node = j->first;
	SgNode* node = j->second;
	string name= in_node->get_qualified_name().str();
	if (debug)
	  cout << " ... checking :  " << name << endl;
	multimap <string, int>::const_iterator k =useresults.begin();
	for (;k!=useresults.end();++k) {
	  string resName = k->first;
	  int resNr = k->second;
	  int tableNr = defuse->getIntForSgNode(node);
	  if (name==resName)
	    if (debug)
	      cout << " ... defNr: " << resNr << "  inTable: " << tableNr <<  endl; 
	  if (name==resName && tableNr==resNr) {
	    hit++;
	    if (debug)
	      cout << " Hit " << hit << "/" << useresults.size() << " - (" << name << "," << resNr << ")" << endl;
	  }
	}

      }
      if (hit!=useresults.size()) {
	cout << " Error: No hit! ... DFA values of node " << nrOfNodes << " are not correct! " << endl;
	exit(1);
      }
      } // if

      }
  if (debug)
    std::cout << "Analysis test is success." << std::endl;
}
    */

void runCurrentFile(vector<string> &argvList, bool debug, bool debug_map) {
  // Build the AST used by ROSE
  if (debug)
    std::cout << ">>>> Starting ROSE frontend ... " << endl;
  SgProject* project = frontend(argvList);
  if (debug)
    std::cout << ">>>> generate PDF " << endl;
  generatePDF ( *project );
  if (debug)
    std::cout << ">>>> start def-use analysis ... " << endl;

  // Call the Def-Use Analysis
  DFAnalysis* defuse = new DefUseAnalysis(project);
  int val = defuse->run(debug);
  /*
  if (debug)
    std::cout << "Analysis is : " << (val ?  "failure" : "success" ) << " " << val <<std::endl;
  if (val==1) exit(1);
  if (debug==false)
    defuse->dfaToDOT();
  */
  //example usage
  // testing
  NodeQuerySynthesizedAttributeType vars = NodeQuery::querySubTree(project, V_SgInitializedName); 
  NodeQuerySynthesizedAttributeType::const_iterator i = vars.begin();
  for (; i!=vars.end();++i) {
    SgInitializedName* initName = isSgInitializedName(*i);
    std::string name = initName->get_qualified_name().str();
  
    vector<SgNode* > vec = defuse->getDefFor(initName, initName);
    if (vec.size()>0)
      if (debug)
	std::cout << "  DEF>> Vector entries for " << name <<  " ( " << 
	  initName << " ) : " << vec.size() << std::endl;

  }
  
  // testing
  vars = NodeQuery::querySubTree(project, V_SgReturnStmt); 
  i = vars.begin();
  typedef std::vector <std::pair <SgInitializedName*, SgNode* > > maptype; 
  for (; i!=vars.end();++i) {
    SgReturnStmt* ret = isSgReturnStmt(*i);
    ROSE_ASSERT(ret);
    
    maptype map = defuse->getUseMultiMapFor(ret);    
    maptype::const_iterator j;
    j = map.begin();
    for (; j!=map.end();++j) {
      SgInitializedName* initName = isSgInitializedName(j->first);
      ROSE_ASSERT(initName);
      std::string name = initName->get_qualified_name().str();
      
      vector<SgNode* > vec = defuse->getUseFor(ret, initName);
      if (vec.size()>0)
	std::cout << "  USE>> Vector entries for " << name <<  " ( " << 
	  ret << " ) : " << vec.size() << std::endl;

    }
    }
  
  // print resulting table
  if (debug_map) {
      cout << "\nDEFMAP" << endl;
    defuse->printDefMap();
      cout << "\nUSEMAP" << endl;
    defuse->printUseMap();
  }
  delete project;
  delete defuse;
}

void usage() {
  cout << " Usage: " << endl;
  cout << "   runTest all [startNr]   --- to test all testcases" << endl;
  cout << "   runTest [file]          --- to test one file" << endl;
  exit(1);
}

int strToInt(string str) {
  int b;
  std::stringstream ss(str);
  ss >> b;
  if (!ss) return 0;
  return b;
}

int main( int argc, char * argv[] )
{
  if (argc==1) {
    usage();
  }
  string pass = argv[1];
  if (pass=="--help") {
    usage();
  }
  bool allTests=false;
  if (pass=="all") {
    allTests=true;
  }

  vector<string> argvList(argv, argv + argc);

  if (allTests==false) {
    multimap <string, int> results;
    multimap <string, int> useresults;

    // change here if you want debug
    bool debug_map = false;
    bool debug = false;
    testOneFunction("::func", argvList, debug, 31, results, useresults);
    //runCurrentFile(argvList, debug, debug_map);
  }  
  
  if (allTests==true) {
    bool debug =false;
    string startNr = "";
    if (argc>2)
      startNr = argv[2];
    bool testAll = false;
    int startNrInt = 0;
    if (startNr=="")
      testAll=true;
    else 
      startNrInt = strToInt(startNr);
    //cout << "start: " << startNrInt << "  all: " << testAll << endl;

    argvList.resize(2);
    argvList[0]=argv[0];
    multimap <string, int> results;
    multimap <string, int> useresults;

    char* srcdirVar = getenv("SRCDIR");
    ROSE_ASSERT (srcdirVar);
    std::string srcdir = srcdirVar;
    srcdir += "/";

    if (startNrInt<=1 || testAll) {
      // ------------------------------ TESTCASE 1 -----------------------------------------
      argvList[1]=srcdir+"tests/test1.C";
      results.clear();      useresults.clear();
      results.insert(pair<string,int>("a", 8));
      results.insert(pair<string,int>("a", 17));    
      results.insert(pair<string,int>("c", 18));
      testOneFunction("::main",argvList, debug, 21, results,useresults);
    }

    if (startNrInt<=2 || testAll) {
      // ------------------------------ TESTCASE 2 -----------------------------------------
      argvList[1]=srcdir+"tests/test2.C";
      results.clear();      useresults.clear();
      results.insert(pair<string,int>("a", 8));
      results.insert(pair<string,int>("a", 17));    
      results.insert(pair<string,int>("c", 18));
      testOneFunction("::main",argvList, debug, 19, results,useresults);
    }

    if (startNrInt<=3 || testAll) {
      // ------------------------------ TESTCASE 3 -----------------------------------------
      argvList[1]=srcdir+"tests/test3.C";
      results.clear();       useresults.clear();
      results.insert(pair<string,int>("d", 15));
      testOneFunction("::main", argvList, debug, 17, results,useresults);
    }

    if (startNrInt<=4 || testAll) {
      // ------------------------------ TESTCASE 4 -----------------------------------------
      argvList[1]=srcdir+"tests/test4.C";
      results.clear();       useresults.clear();
      results.insert(pair<string,int>("::globalvar", 12));
      testOneFunction("::main", argvList, debug, 14, results, useresults);
    }

    if (startNrInt<=5 || testAll) {
      // ------------------------------ TESTCASE 5 -----------------------------------------
      argvList[1]=srcdir+"tests/test5.C";
      results.clear();       useresults.clear();
      results.insert(pair<string,int>("x", 10));
      results.insert(pair<string,int>("y", 16));    
      useresults.insert(pair<string,int>("x", 13));    
      testOneFunction("::main", argvList, debug, 18, results, useresults);
    }

    if (startNrInt<=6 || testAll) {
      // ------------------------------ TESTCASE 6 -----------------------------------------
      argvList[1]=srcdir+"tests/test6.C";
      results.clear();  useresults.clear();
      results.insert(pair<string,int>("x", 14));
      results.insert(pair<string,int>("x", 25));
      results.insert(pair<string,int>("z", 26));    
      results.insert(pair<string,int>("z", 10));    
      useresults.insert(pair<string,int>("x", 16));    
      useresults.insert(pair<string,int>("z", 22));    
      testOneFunction("::main", argvList, debug, 26, results,useresults);
    }

    if (startNrInt<=7 || testAll) {
      // ------------------------------ TESTCASE 7 -----------------------------------------
      argvList[1]=srcdir+"tests/test7.C";
      results.clear();  useresults.clear();
      results.insert(pair<string,int>("i", 6));
      results.insert(pair<string,int>("i", 15));
      useresults.insert(pair<string,int>("i", 8));    
      useresults.insert(pair<string,int>("i", 14));    
      testOneFunction("::main", argvList, debug, 16, results,useresults);
    }

    if (startNrInt<=8 || testAll) {
      // ------------------------------ TESTCASE 8 -----------------------------------------
      argvList[1]=srcdir+"tests/test8.C";
      results.clear();  useresults.clear();
      results.insert(pair<string,int>("x", 24));
      results.insert(pair<string,int>("i", 12));
      results.insert(pair<string,int>("i", 32));
      useresults.insert(pair<string,int>("i", 15));    
      useresults.insert(pair<string,int>("p", 26));    
      useresults.insert(pair<string,int>("x", 28));    
      testOneFunction("::main", argvList, debug, 31, results,useresults);
    }

    if (startNrInt<=9 || testAll) {
      // ------------------------------ TESTCASE 9 -----------------------------------------
      argvList[1]=srcdir+"tests/test9.C";
      results.clear();  useresults.clear();
      results.insert(pair<string,int>("array", 4));
      results.insert(pair<string,int>("i", 12));
      results.insert(pair<string,int>("i", 31));
      useresults.insert(pair<string,int>("i", 15));    
      // useresults.insert(pair<string,int>("array", 20)); // tps: fixed this    
      testOneFunction("::main", argvList, debug, 30, results,useresults);
    }

    if (startNrInt<=10 || testAll) {
      // ------------------------------ TESTCASE 10 -----------------------------------------
      argvList[1]=srcdir+"tests/test10.C";
      results.clear();  useresults.clear();
      results.insert(pair<string,int>("z", 13));
      results.insert(pair<string,int>("z", 22));
      results.insert(pair<string,int>("res", 23));
      testOneFunction("::main", argvList, debug, 25, results,useresults);
      results.clear();  useresults.clear();
      results.insert(pair<string,int>("x", 3));
      results.insert(pair<string,int>("y", 4));
      results.insert(pair<string,int>("z", 5));
      testOneFunction("::f2", argvList, debug, 25, results,useresults);
    }

  }
  
  argvList.clear();
  return 0;
}
