/*
 * File automatically generated by
 * gengen 1.4.2 by Lorenzo Bettini
 * http://www.gnu.org/software/gengen
 */

#ifndef ENUM_DECL_GEN_CLASS_H
#define ENUM_DECL_GEN_CLASS_H

#include <string>
#include <iostream>

using std::string;
using std::ostream;

class enum_decl_gen_class
{
 protected:
  string enum_values;
  string var_arg;

 public:
  enum_decl_gen_class()
  {
  }

  enum_decl_gen_class(const string &_enum_values, const string &_var_arg) :
    enum_values (_enum_values), var_arg (_var_arg)
  {
  }

  static void
  generate_string(const string &s, ostream &stream, unsigned int indent)
  {
    if (!indent || s.find('\n') == string::npos)
      {
        stream << s;
        return;
      }

    string::size_type pos;
    string::size_type start = 0;
    string ind (indent, ' ');
    while ( (pos=s.find('\n', start)) != string::npos)
      {
        stream << s.substr (start, (pos+1)-start);
        start = pos+1;
        if (start+1 <= s.size ())
          stream << ind;
      }
    if (start+1 <= s.size ())
      stream << s.substr (start);
  }

  void set_enum_values(const string &_enum_values)
  {
    enum_values = _enum_values;
  }

  void set_var_arg(const string &_var_arg)
  {
    var_arg = _var_arg;
  }

  void generate_enum_decl(ostream &stream, unsigned int indent = 0);

};

#endif // ENUM_DECL_GEN_CLASS_H
