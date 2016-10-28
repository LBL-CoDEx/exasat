
/*
  Author: Didem Unat 
  Contact: dunat@lbl.gov

  Date Created       : Jun 2012
  Major Modifications: Oct 2012
  Major Modifications: Nov 2013

  * input  : Fortran source code
  * output : XML code description  
  * description: This program is the compiler front-end of the ExaSAT framework
  
  The output includes module, function and loop level information. The program collects ghost cell information, 
  number of flops/point, memory reads and writes. The collected information is useful for optimizing code, modeling performance,  
  selecting arch specs and performance prediction. 
*/  

#include "rose.h"
#include "wholeAST.h"
#include "./functions/FunctionHandler.h"
#include "./modules/ModuleHandler.h"
#include "./utils/ExaSatOptions.h"
#include "./utils/OutputFormatting.h"
#include "./handler/ScopeStmtHandler.h"
#include "./DAG/DAG.h"
#include <boost/python.hpp>
//#define CHILL_ENABLED 1
#ifdef CHILL_ENABLED
extern FILE * yyin;
#include <FlexLexer.h>
#include "loop.hh"
#include <omega.h>
using namespace omega;
//#include "parser.tab.cc"
#include "parser.h"
//namespace {
  enum COMPILER_IR_TYPE {CIT_NULL, CIT_ROSE};
extern  char *source_filename;
extern  COMPILER_IR_TYPE cit_name;
extern  char* procedure_name;

extern  int loop_num_start;
extern  int loop_num_end;
extern  Loop *myloop;
//}

#endif

using namespace std;
using namespace SageInterface;
using namespace OutputFormat;
using namespace boost::python;


string filenameFromPath(string path){
  char* token_prev= new char[path.length()];
  memcpy(token_prev, path.c_str(), path.length());
  char* token = strtok(token_prev, "/");
  while(token){
    token_prev= token;
    token = strtok(NULL, "/");
  }
  string retstr(token_prev);
  return retstr;
}


//reachableFunctions takes a list of all available function definitions in the project and a function definition as the starting point
//It then uses depth first search to return a list of all functions that can be potentially reachable from the initial function
void reachableFunctions(map<string,SgFunctionDefinition*> &funcDefList, SgFunctionDefinition* initFuncDef, map<string, SgFunctionDefinition*> allFuncDefList){
  if(initFuncDef){
     SgBasicBlock* body = initFuncDef->get_body();
     Rose_STL_Container<SgNode*> functionCalls = NodeQuery::querySubTree(body, V_SgFunctionCallExp);
     for(Rose_STL_Container<SgNode*>::iterator it= functionCalls.begin(); it!= functionCalls.end(); it++){
        SgFunctionCallExp* fc= isSgFunctionCallExp(*it);
        if(fc)
        {
           string functionName;
           SgFunctionRefExp* func_exp = isSgFunctionRefExp(fc->get_function());
           if(func_exp){
              SgFunctionSymbol* functionSymbol = NULL;
              functionSymbol = func_exp->get_symbol();
              if(functionSymbol){
                 functionName = functionSymbol->get_name().getString();
              }
           }

           SgFunctionDeclaration* funcDecl= isSgFunctionDeclaration(fc->getAssociatedFunctionDeclaration()->get_definingDeclaration());
           if(!funcDecl)funcDecl= isSgFunctionDeclaration(fc->getAssociatedFunctionDeclaration());
           if(funcDecl)
           {
              string key=functionName;
              SgInitializedNamePtrList argList = funcDecl->get_args();
              Rose_STL_Container<SgInitializedName*>::iterator argIt = argList.begin();
              while(argIt!=argList.end()){
                SgInitializedName* argInitName= isSgInitializedName(*argIt);
                if(argInitName !=NULL){
                  //SgName mangledName= argInitName->get_type()->get_mangled();
                  key+= argInitName->get_type()->unparseToString();//mangledName.getString();
                }
                argIt++;
              }

              SgFunctionDefinition* fd= funcDecl->get_definition();
              if(fd)
              {
                 if(funcDefList.find(key)== funcDefList.end()){
                    funcDefList[key]=fd;
                    reachableFunctions(funcDefList,fd,allFuncDefList);
                 }
              }
              else{//the definition of this function locates in a different file
                 SgFunctionDefinition* aFuncDef= allFuncDefList[key];
                 if(aFuncDef){
                    if(funcDefList.find(key)== funcDefList.end()){
                       funcDefList[key]=aFuncDef;
                       reachableFunctions(funcDefList,aFuncDef,allFuncDefList);
                    }
                 }
              }
           }
        }
     }
  }
}

class collectFunctions : public SgSimpleProcessing
   {
     public:
       void collectFunction(){
         isNoFilter = NULL;
       }
       void setNoFilter(char* val){
	 isNoFilter= val;
       }
       void setLists(map<string, SgFunctionDefinition*> *allList, map<string, SgFunctionDefinition*> *reachableList){
         allFuncDefList= allList;
	 reachableFuncDefList= reachableList;
       }
       virtual void visit (SgNode* node )
       {
         SgFunctionDefinition* funcDef= isSgFunctionDefinition(node);
         if(funcDef){
           SgFunctionDeclaration* funcDecl = funcDef -> get_declaration();
           if(funcDecl){
             SgName funcName = funcDecl->get_name();
             string funcNameStr= funcName.getString();
             string key= funcName.getString();
             //combine function name with type names of all arguments to deal with overloading functions
             SgInitializedNamePtrList argList = funcDecl->get_args();
             Rose_STL_Container<SgInitializedName*>::iterator argIt = argList.begin();
             while(argIt!=argList.end()){
               SgInitializedName* argInitName= isSgInitializedName(*argIt);
               if(argInitName !=NULL){
                 //SgName mangledName= argInitName->get_type()->get_mangled();
                 key+= argInitName->get_type()->unparseToString(); //mangledName.getString();
               }
               argIt++;
             }
             (*allFuncDefList)[key] = funcDef;
             if(isNoFilter){
	       (*reachableFuncDefList)[key]= funcDef;
             }
           }
         }
       }
     private:
       char* isNoFilter;
       map<string, SgFunctionDefinition*> *allFuncDefList;
       map<string, SgFunctionDefinition*> *reachableFuncDefList;
   };


#if 1
//Tan: I copied the loop unrolling implementation from ROSE source
//will implement a new one that works well with the SSA IR generator in the future
bool myloopUnrolling(SgForStatement* target_loop, size_t unrolling_factor)
{
  //Handle 0 and 1, which means no unrolling at all
  if (unrolling_factor <= 1)
    return true;
  // normalize the target loop first
  if (!forLoopNormalization(target_loop));
  {// the return value is not reliable
    //    cerr<<"Error in SageInterface::loopUnrolling(): target loop cannot be normalized."<<endl;
    //    dumpInfo(target_loop);
    //    return false;
  }
  // grab the target loop's essential header information
  SgInitializedName* ivar = NULL;
  SgExpression* lb = NULL;
  SgExpression* ub = NULL;
  SgExpression* step = NULL;
  SgStatement* orig_body = NULL;
  if (!isCanonicalForLoop(target_loop, &ivar, &lb, &ub, &step, &orig_body))
  {
    cerr<<"Error in SageInterface::loopUnrolling(): target loop is not canonical."<<endl;
    dumpInfo(target_loop);
    return false;
  }
  ROSE_ASSERT(ivar&& lb && ub && step);
  ROSE_ASSERT(isSgBasicBlock(orig_body));

   // generate the fringe loop
   bool needFringe = true;
   SgForStatement* fringe_loop = deepCopy<SgForStatement>(target_loop);
   insertStatementAfter(target_loop,fringe_loop);
   removeStatement(fringe_loop->get_for_init_stmt());
   fringe_loop->set_for_init_stmt(NULL);

  // _lu_iter_count = (ub-lb+1)%step ==0?(ub-lb+1)/step: (ub-lb+1)/step+1;
  SgExpression* raw_range_exp =buildSubtractOp(buildAddOp(copyExpression(ub),buildIntVal(1)),
            copyExpression(lb));
  raw_range_exp->set_need_paren(true);
  SgExpression* range_d_step_exp = buildDivideOp(raw_range_exp,copyExpression(step));//(ub-lb+1)/step
  SgExpression* condition_1 = buildEqualityOp(buildModOp(copyExpression(raw_range_exp),copyExpression(step)),buildIntVal(0)); //(ub-lb+1)%step ==0

  SgExpression* iter_count_exp = buildConditionalExp(condition_1,range_d_step_exp, buildAddOp(copyExpression(range_d_step_exp),buildIntVal(1)));
  // fringe = iteration_count%unroll_factor==0 ? 0:unroll_factor*step
  SgExpression* condition_2 = buildEqualityOp(buildModOp(iter_count_exp, buildIntVal(unrolling_factor)), buildIntVal(0));
  SgExpression* initor = buildConditionalExp(condition_2, buildIntVal(0), buildMultiplyOp(buildIntVal(unrolling_factor),copyExpression(step)));

   SgScopeStatement* scope = target_loop->get_scope();
   ROSE_ASSERT(scope != NULL);
   string fringe_name = "_lu_fringe_"+ StringUtility::numberToString(++gensym_counter);
   SgVariableDeclaration* fringe_decl = buildVariableDeclaration(fringe_name, buildIntType(),buildAssignInitializer(initor), scope);
   insertStatementBefore(target_loop, fringe_decl);
   attachComment(fringe_decl, "iter_count = (ub-lb+1)%step ==0?(ub-lb+1)/step: (ub-lb+1)/step+1;");
   attachComment(fringe_decl, "fringe = iter_count%unroll_factor==0 ? 0:unroll_factor*step");

  // compile-time evaluate to see if initor is a constant of value 0
  // if so, the iteration count can be divided even by the unrolling factor
  // and no fringe loop is needed
  // WE have to fold on its parent node to get a possible constant since
  // constant folding only folds children nodes, not the current node to a constant
   ConstantFolding::constantFoldingOptimization(fringe_decl,false);
   SgInitializedName * ivarname = fringe_decl->get_variables().front();
   ROSE_ASSERT(ivarname != NULL);
   // points to a new address if constant folding happens
   SgAssignInitializer * init1 = isSgAssignInitializer(ivarname->get_initializer());
   if (init1)
    if (isSgIntVal(init1->get_operand_i()))
     if (isSgIntVal(init1->get_operand_i())->get_value() == 0)
       needFringe = false;

  // rewrite loop header ub --> ub -fringe; step --> step *unrolling_factor
   SgBinaryOp* ub_bin_op = isSgBinaryOp(ub->get_parent());
   ROSE_ASSERT(ub_bin_op);
   if (needFringe)
     ub_bin_op->set_rhs_operand(buildSubtractOp(copyExpression(ub),buildVarRefExp(fringe_name,scope)));
   else
   {
     ub_bin_op->set_rhs_operand(copyExpression(ub));
     removeStatement(fringe_decl);
   }

   SgBinaryOp* step_bin_op = isSgBinaryOp(step->get_parent());
   ROSE_ASSERT(step_bin_op != NULL);
   step_bin_op->set_rhs_operand(buildMultiplyOp(copyExpression(step),buildIntVal(unrolling_factor)));

   bool isPlus = false;
   if (isSgPlusAssignOp(step_bin_op))
     isPlus = true;
    else if (isSgMinusAssignOp(step_bin_op))
      isPlus = false;
    else
    {
      cerr<<"Error in SageInterface::loopUnrolling(): illegal incremental exp of a canonical loop"<<endl;
      dumpInfo(step_bin_op);
      ROSE_ASSERT(false);
    }

   // copy loop body factor -1 times, and replace reference to ivar  with ivar +/- step*[1 to factor-1]
   for (size_t i =1; i<unrolling_factor; i++)
   {
     SgBasicBlock* body = isSgBasicBlock(deepCopy(fringe_loop->get_loop_body())); // normalized loop has a BB body
     ROSE_ASSERT(body);
     // replace reference to ivar with ivar +/- step*i
     SgExpression* new_exp = NULL;
     std::vector<SgVarRefExp*> refs = querySubTree<SgVarRefExp> (body, V_SgVarRefExp);
     for (std::vector<SgVarRefExp*>::iterator iter = refs.begin(); iter !=refs.end(); iter++)
     {
       SgVarRefExp* refexp = *iter;
       if (refexp->get_symbol()==ivar->get_symbol_from_symbol_table())
       {
         //build replacement  expression if it is NULL
         if (new_exp == NULL)
         {
           if (isPlus) //ivar +/- step * i
           new_exp = buildAddOp(buildVarRefExp(ivar,scope),buildMultiplyOp(copyExpression(step),buildIntVal(i)));
           else
           new_exp = buildSubtractOp(buildVarRefExp(ivar,scope),buildMultiplyOp(copyExpression(step),buildIntVal(i)));

         }
         // replace it with the right one
         replaceExpression(refexp, new_exp);
       }
     }
     // copy body to loop body, this should be a better choice
     // to avoid redefinition of variables after unrolling (new scope is introduced to avoid this)
     appendStatement(body,isSgBasicBlock(orig_body));
    // moveStatementsBetweenBlocks(body,isSgBasicBlock(orig_body));
   }

   // remove the fringe loop if not needed finally
   // it is used to buffering the original loop body before in either cases
   if (!needFringe)
     removeStatement(fringe_loop);

   // constant folding for the transformed AST
   ConstantFolding::constantFoldingOptimization(scope,false);
   //ConstantFolding::constantFoldingOptimization(getProject(),false);
  return true;
}
#endif

int main(int argc, char * argv[])
{

  if(argc <= 1){
    cerr << "ERROR :no input file "<<endl; 
    return -1;
  }
  //Getting program options
  vector<string> argvList(argv, argv + argc);
  ExaSatOptions::GetInstance()->SetExaSatOptions(argvList);

  //Frontend parsing (entirely done in ROSE)
  SgProject *project = frontend (argc, argv);
  if(!SageInterface::is_Fortran_language() && !SageInterface::is_Cxx_language() && !SageInterface::is_C99_language() && !SageInterface::is_C_language()){
    cerr << "ERROR: Input is not a Fortran/C/Cxx program!" << endl;
    return -1;
  }

  char* isNoFilter = getenv("NOFILTER");
  char* initFunctions = getenv("FILTER");
/* True table
 * FILTER             0                     1
 * NOFILTER
 *        0          main/program tree      subtrees
 *        1          Everything             Error
 */


  if(isNoFilter && initFunctions){ //error
    printf("NOFILTER and FILTER can't be set at the same time\n");
    exit(0);
  }

  char* excludeHeader = getenv("NOHEADER");
  //collect all function definitions scattered across source files
  map<string, SgFunctionDefinition*> allFuncDefList;
  map<string, SgFunctionDefinition*> reachableFuncDefList;
  map<string, SgBasicBlock*> basicBlockList;
  if(excludeHeader){//look at specified files only
    collectFunctions c;
    c.setNoFilter(isNoFilter);
    c.setLists(&allFuncDefList, &reachableFuncDefList);
    c.traverseInputFiles (project, preorder);
  }else{//also look at header files 
    for(int f=0; f<project->numberOfFiles(); f++){
      Rose_STL_Container<SgNode*> funcDefList = NodeQuery::querySubTree(project->get_fileList()[f], V_SgFunctionDefinition);
      for(Rose_STL_Container<SgNode*>::iterator it= funcDefList.begin(); it!= funcDefList.end(); it++){
        SgFunctionDefinition* funcDef= isSgFunctionDefinition(*it);
        if(funcDef){
          SgFunctionDeclaration* funcDecl= funcDef->get_declaration();
          if(funcDecl){
            SgName funcName = funcDecl->get_name();
            string funcNameStr= funcName.getString();
            string key= funcName.getString();
            //combine function name with type names of all arguments to deal with overloading functions
            SgInitializedNamePtrList argList = funcDecl->get_args();
            Rose_STL_Container<SgInitializedName*>::iterator argIt = argList.begin();
            while(argIt!=argList.end()){
              SgInitializedName* argInitName= isSgInitializedName(*argIt);
              if(argInitName !=NULL){
                //SgName mangledName= argInitName->get_type()->get_mangled();
                key+= argInitName->get_type()->unparseToString(); //mangledName.getString();
              }
              argIt++;
            }
            allFuncDefList[key] = funcDef;
            if(isNoFilter)reachableFuncDefList[key]= funcDef;
          }
        }
      }
    }
  }

  //if the FILTER variable is set, we look for calling trees started from functions that have that name 
  if(initFunctions)  //subtrees rooted from specified functions
  {
    char* initFunction= strtok(initFunctions, ", ");
    while(initFunction){
      for(int f=0; f<project->numberOfFiles(); f++){
        string token= initFunction;
        size_t pos= token.find(":"); 
        if(pos == token.npos){
          Rose_STL_Container<SgNode*> funcDefList = NodeQuery::querySubTree(project->get_fileList()[f], V_SgFunctionDefinition);
          for(Rose_STL_Container<SgNode*>::iterator it= funcDefList.begin(); it!= funcDefList.end(); it++){
            SgFunctionDefinition* funcDef= isSgFunctionDefinition(*it);
            if(funcDef){
              SgFunctionDeclaration* funcDecl= funcDef->get_declaration();
              if(funcDecl){
                SgName funcName = funcDecl->get_name();
                string funcNameStr= funcName.getString();
                string key=funcNameStr;
	        if(initFunction == funcNameStr){
                  SgInitializedNamePtrList argList = funcDecl->get_args();
                  Rose_STL_Container<SgInitializedName*>::iterator argIt = argList.begin();
                  while(argIt!=argList.end()){
                    SgInitializedName* argInitName= isSgInitializedName(*argIt);
                    if(argInitName !=NULL){
          	      key+= argInitName->get_type()->unparseToString(); //mangledName.getString();
                    }
                    argIt++;
                  }
                  reachableFuncDefList[key]= funcDef;
                  //filter out reachable functions
                  reachableFunctions(reachableFuncDefList, funcDef, allFuncDefList);
	        }
              }
            }
          }
        }else{
          string filename= token.substr(0, pos);
          string lineno= token.substr(pos+1, token.length()); 
          string fname= filenameFromPath(project->get_fileList()[f]->getFileName());
          int dotPos= filename.find(".");
          string namePortion= filename.substr(0, dotPos);
          int dotPos1= filename.find(".");
          string namePortion1= fname.substr(0, dotPos1);
          if(namePortion.compare(namePortion1)==0){//file names match
            //if(filename.substr(dotPos+1, dotPos+1+filename.length()-namePortion.length()).compare(fname.substr(dotPos1+1, dotPos1+1+filename.length()-namePortion1.length()))==0){
              Rose_STL_Container<SgNode*> bbList = NodeQuery::querySubTree(project->get_fileList()[f], V_SgBasicBlock);
              for(Rose_STL_Container<SgNode*>::iterator it= bbList.begin(); it!= bbList.end(); it++){
                SgBasicBlock* bb= isSgBasicBlock(*it);
                if(bb){
                  if(bb->get_file_info()->get_line() == atoi(lineno.c_str())){
                    string key= filename+lineno;
                    basicBlockList[key] = bb; 
                  }
                }
              }
            //}
          }
        }
      }
      initFunction= strtok(NULL, ", ");
    }
  }else if(!isNoFilter){//find main function in C and program in Fortran
    bool Cmain=false;
    bool FortranProgram=false;
    SgFunctionDeclaration* main = findMain(project);
    if(main){
      if(strcmp(strtok((char*)main->unparseToString().c_str(), " "), "PROGRAM")!=0){//we have a C driver
        SgFunctionDefinition* mainDef= main->get_definition();
        ROSE_ASSERT(mainDef);
        //put main in the reachable function list
        string functionName= "main";
        string key=functionName;
        SgInitializedNamePtrList argList = main->get_args();
        Rose_STL_Container<SgInitializedName*>::iterator argIt = argList.begin();
        while(argIt!=argList.end()){
          SgInitializedName* argInitName= isSgInitializedName(*argIt);
          if(argInitName !=NULL){
            //SgName mangledName= argInitName->get_type()->get_mangled();
            key+= argInitName->get_type()->unparseToString(); //mangledName.getString();
          }
          argIt++;
        }
        reachableFuncDefList[key]= main->get_definition();
        //filter out reachable functions
        reachableFunctions(reachableFuncDefList, mainDef, allFuncDefList);
        Cmain=true;
      }
    }
    if(!Cmain){
      //we may have a fortran driver
      //search for the program function
      for(int f=0; f<project->numberOfFiles(); f++){
        Rose_STL_Container<SgNode*> funcDefList = NodeQuery::querySubTree(project->get_fileList()[f], V_SgFunctionDefinition);
        for(Rose_STL_Container<SgNode*>::iterator it= funcDefList.begin(); it!= funcDefList.end(); it++){
          SgFunctionDefinition* funcDef= isSgFunctionDefinition(*it);
          if(funcDef){
 	    string funcCharStr(funcDef->unparseToString());
	    char program[]= "PROGRAM";
	    if(strncmp((char*)funcCharStr.c_str(), program,7)==0){//this is the main function of a fortran program
	      FortranProgram=true;
              SgFunctionDeclaration* funcDecl= funcDef->get_declaration();
              if(funcDecl){
                SgName funcName = funcDecl->get_name();
                string funcNameStr= funcName.getString();
                string key=funcNameStr;
                SgInitializedNamePtrList argList = funcDecl->get_args();
                Rose_STL_Container<SgInitializedName*>::iterator argIt = argList.begin();
                while(argIt!=argList.end()){
                  SgInitializedName* argInitName= isSgInitializedName(*argIt);
                  if(argInitName !=NULL){
          	    key+= argInitName->get_type()->unparseToString(); //mangledName.getString();
                  }
                  argIt++;
                }
                reachableFuncDefList[key]= funcDef;
                //filter out reachable functions
                reachableFunctions(reachableFuncDefList, funcDef, allFuncDefList);
              }
            }
          }
        } 
      }
    }
#if 0
  OutputFormat::outputBegin("program");
  for(map<string,SgFunctionDefinition*>::iterator funcListItr= reachableFuncDefList.begin(); funcListItr != reachableFuncDefList.end(); funcListItr++)
  {
    SgFunctionDefinition* func_def= (*funcListItr).second;
    ROSE_ASSERT(func_def);
    SgFunctionDeclaration* func_dec = func_def -> get_declaration();
    ROSE_ASSERT(func_dec);
    SgName func_name = func_dec->get_name();
    outputBegin("function", func_name.str());
    FunctionHandler* funcHandler = new FunctionHandler();
    funcHandler -> setFortranLanguage(false);
    funcHandler -> process(isSgNode(func_def));
    delete funcHandler;
    outputEnd("function");
  }
  OutputFormat::outputEnd("program");

    //currently we only handle single input file so get the first file 
    SgSourceFile* file = isSgSourceFile(project->get_fileList()[0]);

    OutputFormat::outputBegin("program");  

    //find fortran modules if there is any 
    Rose_STL_Container<SgNode*> modList = NodeQuery::querySubTree(file, V_SgModuleStatement);

    if(modList.size() > 0){
   
      Rose_STL_Container<SgNode*>::iterator m = modList.begin();    
      //process each module in the input file 
      for( ; m != modList.end(); m++)
      {
  	SgModuleStatement* ms = isSgModuleStatement(*m);
  	ROSE_ASSERT(ms);

  	SgName modName = ms->get_name();
  	outputBegin("module", modName.str());
	
  	ModuleHandler* module = new ModuleHandler();
  	module-> processModule(isSgNode(ms));
  	delete module;

  	outputEnd("module");
      }
    }
    else // remove this "else" when getEnclosingModuleStatement in ROSE is implemented
    {
      //find all the fortran subroutines and process them one by one
      Rose_STL_Container<SgNode*> funcList = NodeQuery::querySubTree(file, V_SgFunctionDefinition);
      Rose_STL_Container<SgNode*>::iterator funcListItr = funcList.begin();
      
      for(; funcListItr != funcList.end(); funcListItr++)
  	{
  	  SgFunctionDefinition* func_def= isSgFunctionDefinition(*funcListItr);
  	  ROSE_ASSERT(func_def);

  	  SgFunctionDeclaration* func_dec = func_def -> get_declaration();
  	  ROSE_ASSERT(func_dec);
      
  	  SgName func_name = func_dec->get_name();
  	  outputBegin("function", func_name.str()); 

  	  ModuleHandler::getDependentModuleList(isSgNode(func_def)); //of this function	  
  	  //This function (getEnclosingModuleStatement) is not implemented in Rose yet 
  	  //so instead I check the size of the modules (this means we either 
  	  //have only module or function not a combination of them 
  	  //if (getEnclosingModuleStatement(isSgNode(func)) == NULL)
  	  FunctionHandler* funcHandler = new FunctionHandler();
	  funcHandler -> process(isSgNode(func_def));
	  delete funcHandler; 
  	  //else we already processed this function in processModule routine
  	  outputEnd("function");
  	}
    }
    OutputFormat::outputEnd("program");
#endif
  }else{//include everything
    //and we donot need to do anything here, interesting huh?
  }
  //analyze and generate XML description of reachable functions
  OutputFormat::outputBegin("program");
  for(map<string,SgFunctionDefinition*>::iterator funcListItr= reachableFuncDefList.begin(); funcListItr != reachableFuncDefList.end(); funcListItr++)
  {
    SgFunctionDefinition* func_def= (*funcListItr).second;
    ROSE_ASSERT(func_def);
    SgFunctionDeclaration* func_dec = func_def -> get_declaration();
    ROSE_ASSERT(func_dec);

    FunctionHandler* funcHandler = new FunctionHandler();
    funcHandler -> process(isSgNode(func_def));
    delete funcHandler;
  }
  OutputFormat::outputEnd("program");


  if(!ExaSatOptions::xmlOutput)
    cout << endl <<"  ExaSat: exiting" << endl;
  

  //Generate AST 
  if(ExaSatOptions::generateAST){
    vector<string>  argvList (argv, argv+ argc);
    CustomMemoryPoolDOTGeneration::s_Filter_Flags* filter_flags = new CustomMemoryPoolDOTGeneration::s_Filter_Flags(argvList);
    string filename = SageInterface::generateProjectName(project);
    generateWholeGraphOfAST(filename+"_WholeAST", filter_flags);
  }

#ifdef CHILL_ENABLED
  char* chillScript = getenv("CHILLSCRIPT");
  std::ifstream script;
  if (chillScript) {
    script.open(chillScript);
    if (!script.is_open()) {
      printf("can't open script file \"%s\"\n", chillScript);
      exit(-1);
    }else{
      printf("Successfully opened script file \"%s\"\n", chillScript);
    }
    lexer.switch_streams(&script, &std::cout);
    is_interactive = false;
    ir_code = NULL;
    initializeOmega();
    yydebug = 0;//disable the debug mode
    if (yyparse(project) == 0) {
     if (ir_code != NULL && myloop != NULL && myloop->isInitialized()) {
      if (loop_num_start == loop_num_end) {
        ir_code->ReplaceCode(ir_controls[loops[loop_num_start]], myloop->getCode());
        ir_controls[loops[loop_num_start]] = NULL;
      }
      else {
        std::vector<IR_Control *> parm;
        for (int i = loops[loop_num_start]; i <= loops[loop_num_end]; i++)
          parm.push_back(ir_controls[i]);
        IR_Block *block = ir_code->MergeNeighboringControlStructures(parm);
        ir_code->ReplaceCode(block, myloop->getCode());
        for (int i = loops[loop_num_start]; i <= loops[loop_num_end]; i++) {
          delete ir_controls[i];
          ir_controls[i] = NULL;
        }
      } 
    }
    delete myloop;
    for (int i = 0; i < ir_controls.size(); i++)
      delete ir_controls[i];
    //delete ir_code;
    delete []source_filename;
   }
  }
#endif



  char* ntimes = getenv("UNROLL");
  if(ntimes)
  {
    for(int f=0; f<project->numberOfFiles(); f++){
      Rose_STL_Container<SgNode*> bbList = NodeQuery::querySubTree(project->get_fileList()[f], V_SgBasicBlock);
      for(Rose_STL_Container<SgNode*>::iterator it= bbList.begin(); it!= bbList.end(); it++){
        SgBasicBlock* bb= isSgBasicBlock(*it);
        ROSE_ASSERT(bb);
        if(NodeQuery::querySubTree(bb, V_SgForStatement).size()==0){//we just care about basic blocks that don't contain loops
          if(isSgForStatement(bb->get_parent())){
            bool result = myloopUnrolling(isSgForStatement(bb->get_parent()), atoi(ntimes));
            ROSE_ASSERT(result != false);
          }
        }
      }
    }
  }

  char* isGenerateDAG = getenv("DAG");
  if (isGenerateDAG)
  {
    Py_Initialize();
    SSA_DAG dagGen;
    dagGen.setProject(project);
    dagGen.generateSSA();
    dagGen.ssaToDOT();
    dagGen.generateDAGs(reachableFuncDefList, basicBlockList);
    dagGen.dagsToDOT(reachableFuncDefList, basicBlockList);
  }
  project->unparse();

 return 0;
}


