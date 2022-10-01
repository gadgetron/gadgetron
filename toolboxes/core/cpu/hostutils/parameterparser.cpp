#include "parameterparser.h"

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "log.h"
namespace Gadgetron {

  CommandLineParameter::CommandLineParameter(char com_switch,CommandLineParameterType type, unsigned int nr_values, const char* desc, bool required)
  {
    m_type = type;
    m_switch = com_switch;
    m_nr_values = nr_values;
    m_desc = std::string(desc);
    m_is_required = required;
    m_is_set = false;

    if (m_nr_values > 0){
      m_int_value = new int[m_nr_values];
      m_float_value = new float[m_nr_values];
      m_string_value = new std::string[m_nr_values];
    }
    else{
      m_int_value = 0;
      m_float_value = 0;
      m_string_value = 0;
    }
  }

  CommandLineParameter::~CommandLineParameter()
  {
    if (m_int_value != 0) delete [] m_int_value;
    if (m_float_value != 0) delete [] m_float_value;
    if (m_string_value != 0) delete [] m_string_value;
  } 

  const char* CommandLineParameter::get_string_value(unsigned int i)
  {
    if (i < m_nr_values){
      return m_string_value[i].c_str();
    }
    else{
      return 0;
    }
  }

  int CommandLineParameter::get_int_value(unsigned int i)
  {
    if (i < m_nr_values){
      return m_int_value[i];
    }
    else{
      return 0;
    }
  }

  float CommandLineParameter::get_float_value(unsigned int i)
  {
    if (i < m_nr_values){
      return m_float_value[i];
    }
    else{
      return 0.0f;
    }
  }

  bool CommandLineParameter::get_is_set()
  {
    return m_is_set;
  }

  bool CommandLineParameter::get_is_required()
  {
    return m_is_required;
  }

  bool CommandLineParameter::is_switch_equal_to(char com_switch)
  {
    return (m_switch == com_switch);
  }

  char** CommandLineParameter::set_value(char** argv)
  {
    int args = 0; 
    for (unsigned int i = 0; i < m_nr_values;i++){
      m_string_value[i] = std::string(argv[i]);
      if (m_type == COMMAND_LINE_FLOAT || m_type == COMMAND_LINE_INT){
        std::stringstream ss (std::stringstream::in | std::stringstream::out);
        ss << m_string_value[i];
        ss >> m_float_value[i];
        m_int_value[i] = static_cast<int>(m_float_value[i]);
      }
      else{
        m_int_value[i] = 1;
        m_float_value[i] = 1.0f;
      }
      args++;
    }
    m_is_set = true;

    return (argv+args);
  }

  int CommandLineParameter::get_number_of_values()
  {
    return m_nr_values;
  }

  char CommandLineParameter::get_switch()
  {
    return m_switch;
  }

  std::string CommandLineParameter::get_desc()
  {
    return m_desc;
  }

  ParameterParser::ParameterParser(int list_size, int list_increment)
  {
    m_list_size = list_size;
    m_list_increment = list_increment;
    m_parameter_list = new CommandLineParameter*[m_list_size];
    m_number_of_parameters = 0;
    m_max_desc_length = 0;
    m_max_number_values = 0;
  }

  ParameterParser::~ParameterParser()
  {
    delete_list();
  }

  void ParameterParser::expand_list()
  {
    int new_list_size = m_list_size + m_list_increment;
    CommandLineParameter **new_list = new CommandLineParameter*[new_list_size];

    for (int i = 0; i < m_number_of_parameters; i++){
      new_list[i] = m_parameter_list[i];
    }

    delete [] m_parameter_list;
    m_parameter_list = new_list;
  }

  void ParameterParser::delete_list()
  {
    for (int i = 0; i < m_number_of_parameters; i++){
      delete m_parameter_list[i];
    }
    delete [] m_parameter_list;
  }

  int ParameterParser::add_parameter(char com_switch,CommandLineParameterType type, unsigned int nr_values, 
                                     const char* desc, bool required, const char* def)
  {
    char** argv = new char*[nr_values];
    std::string *arg_list = new std::string[nr_values];

    add_parameter(com_switch, type, nr_values, desc, required);

    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << def;

    unsigned int args = 0; 
    while (args < nr_values){
      ss >> arg_list[args];
      argv[args] = (char*)arg_list[args].c_str();
      args++;
    }

    m_parameter_list[m_number_of_parameters-1]->set_value(argv);
 
    delete [] argv;
    delete [] arg_list;

    return 0;
  }

  int ParameterParser::add_parameter(char com_switch,CommandLineParameterType type, unsigned int nr_values, const char* desc, bool required)
  {
    CommandLineParameter *p = new CommandLineParameter(com_switch, type, nr_values, desc, required);
    for (int i = 0; i < m_number_of_parameters; i++){
      if (m_parameter_list[i]->is_switch_equal_to(com_switch)){
        GINFO(std::string("ParameterParser: Attempt to parameter twice " + std::to_string(com_switch) + "\n").c_str());
        delete p;
        return -1;
      }
    }
    if (m_number_of_parameters >= m_list_size) expand_list();
    m_parameter_list[m_number_of_parameters++] = p;
    if ((int)p->get_desc().length() > m_max_desc_length){
      m_max_desc_length = p->get_desc().length();
    }
    if ((int)p->get_number_of_values() > m_max_number_values) {
      m_max_number_values = p->get_number_of_values();
    }
    return 0;
  }

  int ParameterParser::parse_parameter_list(int argc, char** argv)
  {
    int a = 0;
    int ret = 0;
    m_command_name = std::string(argv[a++]);
    bool argument_found;
    while (a < argc){
      if (argv[a][0] != '-'){
        GINFO("ParameterParser: malformed argument list\n");
        print_usage();
      }

      argument_found = false;
      for (int i = 0; i < m_number_of_parameters; i++){
        if (m_parameter_list[i]->is_switch_equal_to(argv[a][1])){
          if (m_parameter_list[i]->get_number_of_values() <= argc-a-1){
            m_parameter_list[i]->set_value((argv+a+1));
            a += m_parameter_list[i]->get_number_of_values()+1;
            argument_found = true;
            break;
          }
          else{
            GINFO(std::string("\nParameterParser: malformed argument list: -" + std::to_string(argv[a][1]) + "\n").c_str());
            //print_usage();
            argument_found = true;
            a++;
            while (a < argc && argv[a][0] != '-') a++;
            ret = -1;
            break;
          }
        }
      }

      if (!argument_found){
        GINFO(std::string("\n" + std::string("ParameterParser: unknown argument: -") + std::to_string( argv[a][1]) + "\n").c_str());
        //print_usage();
        ret = -1;
        a++;
        while (a < argc && argv[a][0] != '-') a++;
      }
    }
    return ret;
  }
 
  void ParameterParser::print_usage()
  {
    int space_fill = 0;
    std::stringstream usage_stream;

    GINFO("----------------------------------------------------\n");
    usage_stream << "Usage: " << m_command_name << " -[";
    for (int i = 0; i < m_number_of_parameters; i++){
      usage_stream << m_parameter_list[i]->get_switch();
    }
    usage_stream << "]" << "\n";
 
    for (int i = 0; i < m_number_of_parameters; i++){
        usage_stream << " -" << m_parameter_list[i]->get_switch() << " ";
        if (m_max_number_values > 1){
          if (m_parameter_list[i]->get_number_of_values() > 1){
            usage_stream << m_parameter_list[i]->get_number_of_values() << "x "; 
          }
          else{
            usage_stream << " "; 
          }
        }
        if (m_parameter_list[i]->get_number_of_values() > 0){
          usage_stream << "[" << m_parameter_list[i]->get_desc() << "]";
          space_fill = (m_max_desc_length - m_parameter_list[i]->get_desc().length())+2;
        }
        else{
          space_fill = m_max_desc_length+2+2;
        }
        usage_stream << "\n";
      }
    usage_stream << "----------------------------------------------------\n "; 

    GINFO(usage_stream.str().c_str());
  }

  void ParameterParser::print_parameter_list()
  {
    std::stringstream stream;
    stream << "----------------------------------------------------\n";
    for (int i = 0; i < m_number_of_parameters; i++){
      stream << " ";
      stream << "(-" << m_parameter_list[i]->get_switch() << ") ";
      stream << std::setw(m_max_desc_length+2) << std::setiosflags(std::ios::left);
      stream << m_parameter_list[i]->get_desc() << ": ";
      if (m_parameter_list[i]->get_is_set()){
        for (int j = 0; j < m_parameter_list[i]->get_number_of_values(); j++){
          stream << m_parameter_list[i]->get_string_value(j) << " ";
          if (!m_parameter_list[i]->get_is_required()){
            stream << "(optional)";
          }
        }
      }
      else{
        if (m_parameter_list[i]->get_is_required()){
          stream << "(missing)";
        }
        else{
          stream << "N/A (optional)";
        }
      }
      stream << "\n";
    }
    stream << "----------------------------------------------------\n";

    GINFO(stream.str().c_str());
  }
 
  bool ParameterParser::all_required_parameters_set()
  {
    for (int i = 0; i < m_number_of_parameters; i++){
      if (!m_parameter_list[i]->get_is_set() && m_parameter_list[i]->get_is_required())
        return false;
    }
    return true;
  }

  CommandLineParameter* ParameterParser::get_parameter(char com_switch)
  {
    for (int i = 0; i < m_number_of_parameters; i++){
      if (m_parameter_list[i]->is_switch_equal_to(com_switch)){
        return m_parameter_list[i];
      }
    }
    std::stringstream ss;
    ss << "Parameter " << com_switch << " is undefined";
    throw std::runtime_error(ss.str());
  }
}
