#pragma once
#include "hostutils_export.h"

#include <vector>
#include <string>

namespace Gadgetron {

  typedef enum 
    {
      COMMAND_LINE_STRING,
      COMMAND_LINE_INT,
      COMMAND_LINE_FLOAT,
      COMMAND_LINE_NO_ARG
    } CommandLineParameterType; 

  class EXPORTHOSTUTILS CommandLineParameter
  {
  public:
    CommandLineParameter(char com_switch, CommandLineParameterType type, unsigned int nr_values, const char* desc, bool required);
    ~CommandLineParameter();

    bool is_switch_equal_to(char com_switch);

    char** set_value(char** argv);

    int get_number_of_values();
    char get_switch();

    const char* get_string_value(unsigned int i = 0);
    int get_int_value(unsigned int i = 0);
    float get_float_value(unsigned int i = 0);

    bool get_is_set();
    bool get_is_required();
    std::string get_desc();

  private:
    CommandLineParameterType  m_type;
    char                      m_switch;
    unsigned int              m_nr_values;
    std::string               m_desc;
    bool                      m_is_set;
    bool                      m_is_required;
    int                      *m_int_value;
    std::string              *m_string_value;
    float                    *m_float_value;
  };

  class EXPORTHOSTUTILS ParameterParser
  {
  public:
    ParameterParser(int list_size = 10, int list_increment = 10);
    ~ParameterParser();

    int add_parameter(char com_switch,CommandLineParameterType type,  unsigned int nr_values, const char* desc, bool required);
    int add_parameter(char com_switch,CommandLineParameterType type,  unsigned int nr_values, const char* desc, bool required, const char* def);

    int parse_parameter_list(int argc, char** argv);

    int get_number_of_parameters();

    void print_usage();
    void print_parameter_list();

    bool all_required_parameters_set();

    CommandLineParameter* get_parameter(char com_switch);

  private:
    CommandLineParameter** m_parameter_list;
    int m_number_of_parameters;
    int m_list_size;
    int m_list_increment;
    int m_max_desc_length;
    int m_max_number_values;
    std::string m_command_name;
    void expand_list();
    void delete_list();
  };
}
