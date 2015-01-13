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

#include <stdexcept>

#define GADGET_FAIL -1
#define GADGET_OK    0

namespace Gadgetron{

    //Forward declarations
    class GadgetStreamController;

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

            int t = this->get_int_value("threads");
            if (t > 0) {
                GDEBUG("Setting number of threads of gadget %s to %d\n", this->module()->name(), t);
                this->desired_threads(t);
            }

            return this->activate( THR_NEW_LWP | THR_JOINABLE,
                this->desired_threads() );
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

        virtual void set_controller(GadgetStreamController* controller) {
            controller_ = controller;
        }

        virtual GadgetStreamController* get_controller()
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
                    GDEBUG("Gadget (%s) process failed\n", this->module()->name());
                    return GADGET_FAIL;
                }
            }
            return 0;
        }

        int set_parameter(const char* name, const char* val, bool trigger = true) {
	  boost::shared_ptr<std::string> old_value = get_string_value(name);

	    parameter_mutex_.acquire();
            parameters_[std::string(name)] = std::string(val);
	    parameter_mutex_.release();

            if (trigger) {
	      return parameter_changed(std::string(name), std::string(val), *old_value);
            }

            return 0;
        }

        int get_bool_value(const char* name) {
            return (0 == ACE_OS::strcmp(get_string_value(name)->c_str(), "true"));
        }

        int get_int_value(const char* name) {
            return ACE_OS::atoi(get_string_value(name)->c_str());
        }

        double get_double_value(const char* name) {
            return ACE_OS::atof(get_string_value(name)->c_str());
        }

	boost::shared_ptr<std::string> get_string_value(const char* name, unsigned int recursive = 0);

        /**
        *  This trigger function is called whenever set_parameter is called with the trigger = true;
        */
        virtual int parameter_changed(std::string name, std::string new_value, std::string old_value)
        {
            return GADGET_OK;
        }

	const char* get_gadgetron_version() {
	  return gadgetron_version_.c_str();
	}

    protected:
        virtual int next_step(ACE_Message_Block *m)
        {
            return this->put_next(m);//next()->putq(m);
        }

        virtual int process(ACE_Message_Block * m) = 0;

        virtual int process_config(ACE_Message_Block * m) {
            return 0;
        }

        unsigned int desired_threads_;
        bool pass_on_undesired_data_;
        GadgetStreamController* controller_;
	ACE_Thread_Mutex parameter_mutex_;
    private:
        std::map<std::string, std::string> parameters_;
	std::string gadgetron_version_;
    };


    template <class P1> class Gadget1 : public Gadget
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

    template <class P1, class P2> class Gadget2 : public Gadget
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


    template <class P1, class P2, class P3> class Gadget3 : public Gadget
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
