/*
  Didem Unat
  
 */

#include <rose.h>
//#include <iostream>
//#include <string>

#include "Utils.h"

#include "NameGenerator.hh"     //unique variable name generator

using namespace std;
using namespace SageInterface;
using namespace SageBuilder;


static NameGenerator g_var_names ("i_", 0, "_");

//! Hash a string into an unsigned long integer.
static
unsigned long
hashStringToULong (const string& s)
{
  unsigned long sum = 0;
  for (size_t i = 0; i < s.length (); ++i)
    sum += (unsigned long)s[i];
  return sum;
}
 
string Utils::generateVarName(const SgStatement* stmt)
{
  // Generate a prefix.
  stringstream s;
  s << g_var_names.next ();
  const Sg_File_Info* info = stmt->get_startOfConstruct ();
  ROSE_ASSERT (info);

  s << hashStringToULong (info->get_raw_filename ()) ;

  return s.str ();
}

void Utils::stringSplit(string str, string delim, vector<string>& results)
{

  size_t cutAt;
  while( (cutAt = str.find_first_of(delim)) != str.npos )
    {
      if(cutAt > 0)
        {
          results.push_back(str.substr(0,cutAt));
        }
      str = str.substr(cutAt+1);
    }
  if(str.length() > 0)
    results.push_back(str);
}

string Utils::convertIntVecToString(vector<int> int_in)
{
  ostringstream string_out; 
  string_out << "("; 

  for(vector<int>::iterator it = int_in.begin() ; it != int_in.end(); it++)
    {      
      int value = *it ; 
      string_out << value ; 

      //for proper formatting 
      if(it != int_in.end() - 1)
	string_out << ",";
    }

  string_out << ")";
  return string_out.str();
}

string Utils::convertStringVecToString(vector<string> str_in)
{
  ostringstream string_out; 
  string_out << "("; 

  for(vector<string>::iterator it = str_in.begin() ; it != str_in.end(); it++)
    {      
      string value = *it ; 
      string_out << value ; 

      //for proper formatting 
      if(it != str_in.end() - 1)
	string_out << ",";
    }

  string_out << ")";
  return string_out.str();
}

string Utils::removeString( string stringIn, string rem )
{
  string::size_type pos = 0;
  bool spacesLeft = true;

  while( spacesLeft )
    {
      pos = stringIn.find(rem);
      if( pos != string::npos )
        stringIn.erase( pos, 1 );
      else
        spacesLeft = false;
    }

  return stringIn;
}


// eof


