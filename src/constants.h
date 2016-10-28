
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <set>
using namespace std;

const string  math_func_table[]= {"sqrt", "dsqrt", "csqrt", 
				  "exp", "dexp", "cexp", 
				  "log","alog", "dlog", "clog",
				  "log10", "alog10", "dlog10", 
				  "sin", "dsin", "csin",
				  "cos", "dcos", "ccos",
				  "tan", "dtan", "asin", "dasin", "acos", "dacos",
				  "atan", "dtan", "atan2", "datan2", 
				  "sinh", "dsinh", "cosh", "dcosh", "tanh", "dtanh", "mod"};

const string  math_low_lat_func_table[]= {"abs", "max", "min", "floor", "ceiling"};

const string fortran_intrinsic_func_table[] ={"abort", 
					      "abs",
					      "access",
					      "achar",
					      "acos",
					      "acosh",
					      "adjustl","aimag", "aint", "alarm", "all", "allocated", "and", 
					      "anint", "any", "asin", "asinh", "associated", "atan", "atan2", "atanh",
					      "atomic_define", "atomic_ref", //23
					      "present",
					      "merge"
					      
};

const std::set<string >math_func_list (math_func_table, math_func_table + 36);
const std::set<string >math_func_low_lat_list (math_low_lat_func_table, math_low_lat_func_table + 5);
const std::set<string >fortran_intrinsic_func_list (fortran_intrinsic_func_table, fortran_intrinsic_func_table + 25);

const int DEFAULT_BOX_SIZE= 64; 

const double Mbytes = 9.53674e-07;

const bool CACHE_BYPASS = true ; //false; 

#endif
