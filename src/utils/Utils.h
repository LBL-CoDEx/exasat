/*
  by Didem Unat 
  */

#ifndef _NOTAUTILS_H
#define _NOTAUTILS_H

//#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;
using namespace SageInterface;

namespace Utils
{
  string generateVarName(const SgStatement* stmt);
  void stringSplit(string str, string delim, vector<string>& results);
  string removeString( string stringIn, string rem );

  string convertIntVecToString(vector<int> int_in);
  string convertStringVecToString(vector<string> str_in);
};

#endif //_NOTAUTILS_H
