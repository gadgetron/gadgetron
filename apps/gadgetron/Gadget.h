#ifndef GADGET_H
#define GADGET_H
#pragma once

#include <ace/OS_NS_stdlib.h>
#include <ace/Task.h>
#include <ace/Stream.h>
#include <ace/Module.h>
#include <ace/OS_Memory.h>
#include <ace/Svc_Handler.h>
#include <ace/SOCK_Stream.h>

#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "gadgetbase_export.h"
#include "GadgetContainerMessage.h"
#include "GadgetronExport.h"
#include "gadgetron_config.h"
#include "log.h"
#include <initializer_list>

#include <stdexcept>

#define GADGET_FAIL -1
#define GADGET_OK    0

namespace Gadgetron{

  class GadgetPropertyBase
  {
  public:
    GadgetPropertyBase(const char* name, const char* type_string, const char* description)
    : name_(name)
    , type_str_(type_string)
    , description_(description)
    , str_value_("")
    , is_reference_(false)
    , reference_gadget_("")
    , reference_property_("")
    {

    }

    virtual const char* name()
    {
      return name_.c_str();
    }

    virtual const char* string_value()
    {
      return str_value_.c_str();
    }

    virtual void string_value(const char* value)
    {
      str_value_ = value;
      size_t at_pos = str_value_.find('@');
      if (at_pos != std::string::npos) {
        //There was an add sign, which means look for that parameter on another gadget
        std::string reference_property_ = str_value_.substr(0,at_pos);
        std::string reference_gadget_   = str_value_.substr(at_pos+1);
        is_reference_ = true;
      }
    }

    virtual const char* type_string()
    {
      return type_str_.c_str();
    }

    virtual const char* description()
    {
      return description_.c_str();
    }

    virtual const char* limits_description()
    {
      return "";
    }

  protected:
    std::string name_;
    std::string type_str_;
    std::string description_;
    std::string str_value_;
    bool is_reference_;
    std::string reference_gadget_;
    std::string reference_property_;
  };

  //Forward declarations
  class GadgetStreamInterface;

  class EXPORTGADGETBASE Gadget : public ACE_Task<ACE_MT_SYNCH>
  {

  public:
    typedef ACE_Task<ACE_MT_SYNCH> inherited;

    enum
    {
      GADGET_MESSAGE_CONFIG = (ACE_Message_Block::USER_FLAGS << 1)
    };

    Gadget()
    : inherited()
    , desired_threads_(1)
    , pass_on_undesired_data_(false)
    , controller_(0)
    , parameter_mutex_("GadgetParameterMutex")
    {

      gadgetron_version_ = std::string(GADGETRON_VERSION_STRING) + std::string(" (") +
      std::string(GADGETRON_GIT_SHA1_HASH) + std::string(")");

    }

    virtual ~Gadget()
    {
      if (this->module()) {
        GDEBUG("Shutting down Gadget (%s)\n", this->module()->name());
      }
    }


    virtual int init(void)
    {
      return 0;
    }

    virtual int open(void* = 0)
    {
      return this->activate( THR_NEW_LWP | THR_JOINABLE, this->desired_threads() );
    }

    int put(ACE_Message_Block *m, ACE_Time_Value* timeout = 0)
    {
      return this->putq(m, timeout);
    }

    virtual unsigned int desired_threads()
    {
      return desired_threads_;
    }

    virtual void desired_threads(unsigned int t)
    {
      desired_threads_ = t;
    }

    virtual bool pass_on_undesired_data()
    {
      return pass_on_undesired_data_;
    }

    virtual void pass_on_undesired_data(bool d)
    {
      pass_on_undesired_data_ = d;
    }

    virtual void set_controller(GadgetStreamInterface* controller) {
      controller_ = controller;
    }

    virtual GadgetStreamInterface* get_controller()
    {
      return controller_;
    }

    virtual int close(unsigned long flags)
    {
      GDEBUG("Gadget (%s) Close Called with flags = %d\n", this->module()->name(), flags);
      int rval = 0;
      if (flags == 1) {
        ACE_Message_Block *hangup = new ACE_Message_Block();
        hangup->msg_type( ACE_Message_Block::MB_HANGUP );
        if (this->putq(hangup) == -1) {
          hangup->release();
          GDEBUG("Gadget (%s) failed to put hang up message on queue\n", this->module()->name());
          return GADGET_FAIL;
        }
        GDEBUG("Gadget (%s) waiting for thread to finish\n", this->module()->name());
        rval = this->wait();
        GDEBUG("Gadget (%s) thread finished\n", this->module()->name());
        controller_ = 0;
      }
      return rval;
    }

    virtual int svc(void)
    {
      for (ACE_Message_Block *m = 0; ;) {

        //GDEBUG("Waiting for message in Gadget (%s)\n", this->module()->name());
        if (this->getq(m) == -1) {
          GDEBUG("Gadget (%s) failed to get message from queue\n", this->module()->name());
          return GADGET_FAIL;
        }
        //GDEBUG("Message Received in Gadget (%s)\n", this->module()->name());

        //If this is a hangup message, we are done, put the message back on the queue before breaking
        if (m->msg_type() == ACE_Message_Block::MB_HANGUP) {
          //GDEBUG("Gadget (%s) Hangup message encountered\n", this->module()->name());
          if (this->putq(m) == -1) {
            GDEBUG("Gadget (%s) failed to put hang up message on queue (for other threads)\n", this->module()->name());
            return GADGET_FAIL;
          }
          //GDEBUG("Gadget (%s) breaking loop\n", this->module()->name());
          break;
        }


        //Is this config info, if so call appropriate process function
        if (m->flags() & GADGET_MESSAGE_CONFIG) {

          int success;
          try{ success = this->process_config(m); }
          catch (std::runtime_error& err){
            GEXCEPTION(err,"Gadget::process_config() failed\n");
            success = -1;
          }

          if (success == -1) {
            m->release();
            this->flush();
            GDEBUG("Gadget (%s) process config failed\n", this->module()->name());
            return GADGET_FAIL;

          }

          //Push this onto next gadgets queue, other gadgets may need this configuration information
          if (this->next()) {
            if (this->next()->putq(m) == -1) {
              m->release();
              GDEBUG("Gadget (%s) process config failed to put config on dowstream gadget\n", this->module()->name());
              return GADGET_FAIL;
            }
          }
          continue;
        }

        int success;
        try{ success = this->process(m); }
        catch (std::runtime_error& err){
          GEXCEPTION(err,"Gadget::process() failed\n");
          success = -1;
        }

        if (success == -1) {
          m->release();
          this->flush();
          GERROR("Gadget (%s) process failed\n", this->module()->name());
          return GADGET_FAIL;
        }
      }
      return 0;
    }

    virtual int set_parameter(const char* name, const char* val, bool trigger = true) {
      boost::shared_ptr<std::string> old_value = get_string_value(name);
      GadgetPropertyBase* p = this->find_property(name);

      if (p) {
        p->string_value(val);
      } else {
        throw std::runtime_error("Attempting to set non-registered property");
      }

      parameter_mutex_.acquire();
      parameters_[std::string(name)] = std::string(val);
      parameter_mutex_.release();

      if (trigger) {
        return parameter_changed(std::string(name), std::string(val), *old_value);
      }

      return 0;
    }

    /*
    virtual int get_bool_value(const char* name) {
    return (0 == ACE_OS::strcmp(get_string_value(name)->c_str(), "true"));
    }

    virtual int get_int_value(const char* name) {
      return ACE_OS::atoi(get_string_value(name)->c_str());
    }

    virtual double get_double_value(const char* name) {
      return ACE_OS::atof(get_string_value(name)->c_str());
    }
    */

    boost::shared_ptr<std::string> get_string_value(const char* name, unsigned int recursive = 0);

    /**
    *  This trigger function is called whenever set_parameter is called with the trigger = true;
    */
    virtual int parameter_changed(std::string name, std::string new_value, std::string old_value)
    {
      return GADGET_OK;
    }

    void print_properties()
    {
      for (std::vector<GadgetPropertyBase*>::iterator it = properties_.begin(); it != properties_.end(); it++)
      {
        GDEBUG("Parameter with name: %s\n", (*it)->name());
      }
    }

    int get_number_of_properties()
    {
      return properties_.size();
    }

    GadgetPropertyBase* get_property_by_index(size_t i)
    {
      if (i >= properties_.size()) {
        return 0;
      }
      return properties_[i];
    }

    GadgetPropertyBase* find_property(const char* name)
    {
      GadgetPropertyBase* p = 0;
      parameter_mutex_.acquire();
      for (std::vector<GadgetPropertyBase*>::iterator it = properties_.begin(); it != properties_.end(); it++) {
        if (std::string(name) == std::string((*it)->name())) {
          p = *it;
          break;
        }
      }
      parameter_mutex_.release();
      return p;
    }
    void register_property(GadgetPropertyBase* p)
    {
      parameter_mutex_.acquire();
      properties_.push_back(p);
      parameter_mutex_.release();
    }

    const char* get_gadgetron_version() {
      return gadgetron_version_.c_str();
    }

  protected:
    std::vector<GadgetPropertyBase*> properties_;

    virtual int next_step(ACE_Message_Block *m)
    {
      return this->put_next(m);//next()->putq(m);
    }

    virtual int process(ACE_Message_Block * m) = 0;

    virtual int process_config(ACE_Message_Block * m) {
      return 0;
    }

    unsigned int desired_threads_;
    unsigned int threads_;
    bool pass_on_undesired_data_;
    GadgetStreamInterface* controller_;
    ACE_Thread_Mutex parameter_mutex_;
  private:
    std::map<std::string, std::string> parameters_;
    std::string gadgetron_version_;
  };


  template <typename T> class GadgetPropertyLimits
  {
  public:
    virtual bool within_limits(T& v) = 0;
    virtual const char* limits_description() = 0;
  };

  template <typename T> class GadgetPropertyLimitsNoLimits
  : public GadgetPropertyLimits<T>
  {
  public:
    virtual bool within_limits(T& v) {
      return true;
    }

    virtual const char* limits_description() {
      return "";
    }
  };

  template <typename T> class GadgetPropertyLimitsEnumeration
  : public GadgetPropertyLimits<T>
  {
  public:
    GadgetPropertyLimitsEnumeration(std::initializer_list<T> valid_vals) {
      valid_vals_.insert(valid_vals_.end(), valid_vals.begin(), valid_vals.end());
    }

    virtual bool within_limits(T& v)
    {
      typename std::vector<T>::iterator it;
      it = find(valid_vals_.begin(), valid_vals_.end(), v);
      if (it != valid_vals_.end()) return true;
      return false;
    }

    virtual const char* limits_description()
    {
      if (!limits_desc_.size()) {
        std::stringstream strstream;
        typename std::vector<T>::iterator it;
        it = valid_vals_.begin();
        if (it != valid_vals_.end()) {
          strstream << "[";
          strstream << *it;
          it++;
          while (it != valid_vals_.end()) {
            strstream << ", " << *it;
            it++;
          }
          strstream << "]";
        }
        limits_desc_ = strstream.str();
      }
      return limits_desc_.c_str();
    }

  protected:
    std::vector<T> valid_vals_;
    std::string limits_desc_;
  };

  template <typename T> class GadgetPropertyLimitsRange
  : public GadgetPropertyLimits<T>
  {
  public:
    GadgetPropertyLimitsRange(T min_val, T max_val)
    : min_(min_val)
    , max_(max_val)
    {
    }

    virtual bool within_limits(T& v)
    {
      return ( (v >= min_) && (v <= max_) );
    }

    virtual const char* limits_description()
    {
      if (!limits_desc_.size()) {
        std::stringstream strstream;
        strstream << "[" << min_ << ":" << max_ << "]" << std::endl;
        limits_desc_ = strstream.str();
      }
      return limits_desc_.c_str();
    }

  protected:
    T min_;
    T max_;
    std::string limits_desc_;
  };
  template <typename T, typename L> class GadgetProperty
  : public GadgetPropertyBase
  {
  public:
    GadgetProperty(const char* name, const char* type_string, const char* description,
      Gadget* g, T default_value, L limits)
      : GadgetPropertyBase(name,type_string,description)
      , g_(g)
      , limits_(limits)
      {
        g_->register_property(this);
        this->value(default_value);
      }

      T value()
      {
        if (is_reference_) {
          boost::shared_ptr<std::string> val = this->g_->get_string_value(this->name());
          std::stringstream(*val) >> std::boolalpha >> value_;
        }
        return value_;
      }

      void value(T v)
      {
        value_ = v;
        std::stringstream strstream;
        strstream << std::boolalpha << v;
        strstream >> str_value_;
        is_reference_ = false;
        if (!limits_.within_limits(v)) {
          GERROR("Property: %s, value: %s, limits:%s\n", this->name(), str_value_.c_str(), this->limits_.limits_description());
          throw std::runtime_error("Value assigned outside limit range");
        }
      }

      virtual void string_value(const char* val)
      {
        GadgetPropertyBase::string_value(val);

        if (!is_reference_)
        {
          T tmp;
          std::stringstream(val) >> std::boolalpha >> tmp;
          this->value(tmp);
        }
      }


      bool operator==(const T &v) const
      {
        return this->value() == v;
      }

      virtual const char* limits_description()
      {
        return limits_.limits_description();
      }

    protected:
      T value_;
      L limits_;
      Gadget* g_;
    };
    template <typename T, typename L> class GadgetProperty<std::vector<T>,L >
    : public GadgetPropertyBase
    {
    public:
      GadgetProperty(const char* name, const char* type_string, const char* description,
        Gadget* g, std::initializer_list<T> default_value, L limits)
        : GadgetPropertyBase(name,type_string,description)
        , g_(g)
        , limits_(limits)
        {
          g_->register_property(this);
          this->value(default_value);
        }

        std::vector<T> value()
        {
          if (is_reference_) {
            boost::shared_ptr<std::string> val = this->g_->get_string_value(this->name());
            std::stringstream ss(*val);
            T chunk;
            values_ = std::vector<T>();
            while (ss >> std::boolalpha >> chunk)
            values_.push_back(chunk);
          }
          return values_;
        }

        void value(std::vector<T> v)
        {
          values_ = v;
          std::stringstream strstream;
          for (T val : values_)
          strstream << std::boolalpha << val << " ";
          strstream >> str_value_;
          is_reference_ = false;
          if (!limits_.within_limits(v)) {
            GERROR("Property: %s, value: %s, limits:%s\n", this->name(), str_value_.c_str(), this->limits_.limits_description());
            throw std::runtime_error("Value assigned outside limit range");
          }
        }

        virtual void string_value(const char* val)
        {
          GadgetPropertyBase::string_value(val);

          if (!is_reference_)
          {
            std::vector<T> tmp;
            T chunk;
            std::stringstream ss(val);
            while (ss >> std::boolalpha >> chunk )
            tmp.push_back(chunk);
            this->value(tmp);
          }
        }


        bool operator==(const std::vector<T> &v) const
        {
          return this->value() == v;

        }

        virtual const char* limits_description()
        {
          return limits_.limits_description();
        }

      protected:
        std::vector<T> values_;
        L limits_;
        Gadget* g_;
      };

      #define GADGET_PROPERTY(varname, vartype, description, defaultvalue) GadgetProperty<vartype, GadgetPropertyLimitsNoLimits<vartype> > varname{#varname,#vartype, description, this, defaultvalue, GadgetPropertyLimitsNoLimits<vartype>()}
      #define GADGET_PROPERTY_LIMITS(varname, vartype, description, defaultvalue, limitstype, ...) GadgetProperty<vartype, limitstype<vartype> > varname{#varname,#vartype, description, this, defaultvalue, limitstype<vartype>{ __VA_ARGS__ }}

      class BasicPropertyGadget : public Gadget
      {
      public:
        BasicPropertyGadget()
        : Gadget()
        {
          desired_threads_ = threads.value();
          pass_on_undesired_data_ = pass_on_undesired_data.value();
        }

        virtual ~BasicPropertyGadget() {}

      protected:
        GADGET_PROPERTY(using_cloudbus,bool,"Indicates whether the cloudbus is in use and available", false);
        GADGET_PROPERTY(pass_on_undesired_data,bool, "If true, data not matching the process function will be passed to next Gadget", false);
        GADGET_PROPERTY(threads,int, "Number of threads to run in this Gadget", 1);
        #ifdef _WIN32
        GADGET_PROPERTY(workingDirectory, std::string, "Where to store temporary files", "c:\\temp\\gadgetron\\");
        #else
        GADGET_PROPERTY(workingDirectory, std::string, "Where to store temporary files", "/tmp/gadgetron/");
        #endif // _WIN32
      };

      template <class P1> class Gadget1 : public BasicPropertyGadget
      {

      protected:
        int process(ACE_Message_Block* mb)
        {
          GadgetContainerMessage<P1>* m = AsContainerMessage<P1>(mb);

          if (!m) {
            if (!pass_on_undesired_data_) {
              GERROR("Gadget1::process, conversion of message block");
              return -1;
            } else {
              return (this->next()->putq(mb));
            }

          }

          return this->process(m);
        }

        virtual int process(GadgetContainerMessage<P1>* m) = 0;

      };

      template <class P1, class P2> class Gadget2 : public BasicPropertyGadget
      {

      protected:
        int process(ACE_Message_Block* mb)
        {

          GadgetContainerMessage<P1>* m1 = AsContainerMessage<P1>(mb);

          GadgetContainerMessage<P2>* m2 = 0;
          if (m1) {
            m2 = AsContainerMessage<P2>(m1->cont());
          }

          if (!m1 || !m2) {
            if (!pass_on_undesired_data_) {
              GERROR("%s -> %s, (%s, %s, %p, %p), (%s, %s, %p, %p)\n",
              this->module()->name(),
              "Gadget2::process, Conversion of Message Block Failed",
              typeid(GadgetContainerMessage<P1>*).name(),
              typeid(m1).name(),
              mb,
              m1,
              typeid(GadgetContainerMessage<P2>*).name(),
              typeid(m2).name(),
              mb->cont(),
              m2);
              return -1;
            } else {
              return (this->next()->putq(mb));
            }
          }

          return this->process(m1,m2);
        }

        virtual int process(GadgetContainerMessage<P1>* m1, GadgetContainerMessage<P2>* m2) = 0;

      };


      template <class P1, class P2, class P3> class Gadget3 : public BasicPropertyGadget
      {

      protected:
        int process(ACE_Message_Block* mb)
        {

          GadgetContainerMessage<P1>* m1 = AsContainerMessage<P1>(mb);

          GadgetContainerMessage<P2>* m2 = 0;
          if (m1) {
            m2 = AsContainerMessage<P2>(m1->cont());
          }

          GadgetContainerMessage<P3>* m3 = 0;
          if (m2) {
            m3 = AsContainerMessage<P3>(m2->cont());
          }

          if (!m1 || !m2 || !m3) {
            if (!pass_on_undesired_data_) {
              GERROR("%s -> %s, (%s, %s, %p), (%s, %s, %p), (%s, %s, %p)\n",
              this->module()->name(),
              "Gadget3::process, Conversion of Message Block Failed",
              typeid(GadgetContainerMessage<P1>*).name(),
              typeid(m1).name(),
              m1,
              typeid(GadgetContainerMessage<P2>*).name(),
              typeid(m2).name(),
              m2,
              typeid(GadgetContainerMessage<P3>*).name(),
              typeid(m3).name(),
              m3);
              return -1;
            } else {
              return (this->next()->putq(mb));
            }
          }

          return this->process(m1,m2,m3);
        }

        virtual int process(GadgetContainerMessage<P1>* m1, GadgetContainerMessage<P2>* m2, GadgetContainerMessage<P3>* m3) = 0;

      };

      /* Macros for handling dyamic linking */
      // #define GADGET_DECLARE(GADGET) GADGETRON_LOADABLE_DECLARE(GADGET)
      // #define GADGET_FACTORY_DECLARE(GADGET) GADGETRON_LOADABLE_FACTORY_DECLARE(Gadget,GADGET)

      #define GADGET_DECLARE(GADGET)
      #define GADGET_FACTORY_DECLARE(GADGET) GADGETRON_LOADABLE_FACTORY_DECLARE(Gadget,GADGET)

    }

    #endif //GADGET_H
