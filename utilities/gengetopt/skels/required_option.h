/*
 * File automatically generated by
 * gengen 1.2 by Lorenzo Bettini
 * http://www.gnu.org/software/gengen
 */

#ifndef REQUIRED_OPTION_GEN_CLASS_H
#define REQUIRED_OPTION_GEN_CLASS_H

#include <string>
#include <iostream>

using std::string;
using std::ostream;

class required_option_gen_class
{
 protected:
  bool checkrange;
  string mode_condition;
  string option_descr;
  string option_var_name;
  string package_var_name;

 public:
  required_option_gen_class() :
    checkrange (false)
  {
  }

  required_option_gen_class(bool _checkrange, const string &_mode_condition, const string &_option_descr, const string &_option_var_name, const string &_package_var_name) :
    checkrange (_checkrange), mode_condition (_mode_condition), option_descr (_option_descr), option_var_name (_option_var_name), package_var_name (_package_var_name)
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

  void set_checkrange(bool _checkrange)
  {
    checkrange = _checkrange;
  }

  void set_mode_condition(const string &_mode_condition)
  {
    mode_condition = _mode_condition;
  }

  void set_option_descr(const string &_option_descr)
  {
    option_descr = _option_descr;
  }

  void set_option_var_name(const string &_option_var_name)
  {
    option_var_name = _option_var_name;
  }

  void set_package_var_name(const string &_package_var_name)
  {
    package_var_name = _package_var_name;
  }

  void generate_required_option(ostream &stream, unsigned int indent = 0);

};

#endif // REQUIRED_OPTION_GEN_CLASS_H
