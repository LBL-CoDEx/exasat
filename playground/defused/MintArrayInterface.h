

#include "rose.h"
//#include "ASTtools.hh"
#include <vector>


using namespace std;
using namespace SageInterface;

class MintArrayInterface
{

 public:

  static bool isArrayReference(SgExpression* ref, 
			       SgExpression** arrayName/*=NULL*/, 
			       std::vector<SgExpression*>* subscripts/*=NULL*/);

  static void linearizeArrays(SgFunctionDeclaration* node);
  static void linearizeArrays(SgSourceFile* file);

 private:


};
