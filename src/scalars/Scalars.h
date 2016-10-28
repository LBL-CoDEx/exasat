/*
  by Didem Unat 
  */

#ifndef _SCALARS_H
#define _SCALARS_H

//#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>
#include "../datatypes.h"

using namespace std;
using namespace SageInterface;

class Scalars
{
public:

  static void getActiveScalarVars(SgNode* node, 
				  std::set<SgVarRefExp*>& scalars);

  static void  getScalarReadWriteCounts(ScalarRWMap_t& scalarRW, 
					FirstAccessMap_t& firstAccess,
					std::set<SgVarRefExp*>& diffScalars);

  static void getScalarInames(std::set<SgInitializedName*>& scalarIname, 
			      std::set<SgVarRefExp*>& diffScalars);
    
};

#endif //_SCALARS_H
