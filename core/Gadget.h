#pragma once

#include <boost/dll/alias.hpp>
#include <boost/shared_ptr.hpp>
#include <initializer_list>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>

#include "mrd/types.h"

#include "Channel.h"
#include "Context.h"
#include "GadgetContainerMessage.h"
#include "Node.h"
#include "log.h"

#define GADGET_FAIL -1
#define GADGET_OK    0


namespace Gadgetron {

    class GadgetPropertyBase {
    public:
        GadgetPropertyBase(const char *name, const char *type_string, const char *description)
                : name_(name), type_str_(type_string), description_(description), str_value_(""), is_reference_(false),
                  reference_gadget_(""), reference_property_("") {

        }

        virtual const char *name() const {
            return name_.c_str();
        }

        virtual const char *string_value() {
            return str_value_.c_str();
        }

        virtual void string_value(const char *value) {
            str_value_ = value;
            size_t at_pos = str_value_.find('@');
            if (at_pos != std::string::npos) {
                //There was an add sign, which means look for that parameter on another gadget
                std::string reference_property_ = str_value_.substr(0, at_pos);
                std::string reference_gadget_ = str_value_.substr(at_pos + 1);
                is_reference_ = true;
            }
        }

        virtual const char *type_string() {
            return type_str_.c_str();
        }

        virtual const char *description() {
            return description_.c_str();
        }

        virtual const char *limits_description() {
            return "";
        }

        virtual ~GadgetPropertyBase() = default;

    protected:
        std::string name_;
        std::string type_str_;
        std::string description_;
        std::string str_value_;
        bool is_reference_;
        std::string reference_gadget_;
        std::string reference_property_;
    };

    class ChannelAdaptor {
    public:
        explicit ChannelAdaptor(Core::OutputChannel& out) : channel(out) { }

        int putq(GadgetContainerMessageBase *msg) {
            if (msg) {
                channel.push_message(to_message(msg));
            }
            msg->release();
            return 0;
        }

    private:
        static Core::Message to_message(GadgetContainerMessageBase *message) {
            std::vector<std::unique_ptr<Core::MessageChunk>> messages;
            auto *current_message = message;
            while (current_message) {
                messages.emplace_back(current_message->take_message());
                current_message = current_message->cont();
            }
            return Core::Message(std::move(messages));
        }

        Core::OutputChannel& channel;
    };

    class Gadget {
    public:
        Gadget() {}

        virtual ~Gadget() {
            GDEBUG_STREAM("Shutting down Gadget (" << this->name_ << ")");
        }

        virtual int set_parameter(const char *name, const char *val, bool trigger = true) {
            std::string old_value = get_string_value(name);
            GadgetPropertyBase *p = this->find_property(name);

            if (p) {
                p->string_value(val);
            } else {
                throw std::runtime_error("Attempting to set non-registered property");
            }

            std::lock_guard<std::mutex> guard(parameter_mutex_);
            parameters_[std::string(name)] = std::string(val);

            if (trigger) {
                return parameter_changed(std::string(name), std::string(val), old_value);
            }

            return 0;
        }

        int putq(GadgetContainerMessageBase *msg) {
            return this->next_channel->putq(msg);
        }

        virtual int close(unsigned long flags = 1) { return 0; }

        void set_context(const Core::Context& context);

        std::string get_string_value(const char *name);

        /**
        *  This trigger function is called whenever set_parameter is called with the trigger = true;
        */
        virtual int parameter_changed(std::string name, std::string new_value, std::string old_value) {
            return GADGET_OK;
        }

        void print_properties() {
            for (std::vector<GadgetPropertyBase *>::iterator it = properties_.begin(); it != properties_.end(); it++) {
                GDEBUG("Parameter with name: %s\n", (*it)->name());
            }
        }

        int get_number_of_properties() {
            return properties_.size();
        }

        GadgetPropertyBase *get_property_by_index(size_t i) {
            if (i >= properties_.size()) {
                return 0;
            }
            return properties_[i];
        }

        inline GadgetPropertyBase *find_property(const char *name) {
            GadgetPropertyBase *p = 0;
            std::lock_guard<std::mutex> guard(parameter_mutex_);
            for (std::vector<GadgetPropertyBase *>::iterator it = properties_.begin(); it != properties_.end(); it++) {
                if (std::string(name) == std::string((*it)->name())) {
                    p = *it;
                    break;
                }
            }
            return p;
        }

        void register_property(GadgetPropertyBase *p) {
            std::lock_guard<std::mutex> guard(parameter_mutex_);
            properties_.push_back(p);
        }

        void next(std::shared_ptr<ChannelAdaptor> n) {
            next_channel = n;
        }

        std::shared_ptr<ChannelAdaptor> next() {
            return next_channel;
        }

        void set_name(const std::string& name) {
            name_ = name;
        }

        std::vector<GadgetPropertyBase *> properties_;

        virtual int process(GadgetContainerMessageBase *m) = 0;

        virtual int process_config(const mrd::Header& header) {
            return 0;
        }

    protected:
        std::string name_;
        std::mutex parameter_mutex_;
        Core::Context context;

    private:
        std::map<std::string, std::string> parameters_;
        std::shared_ptr<ChannelAdaptor> next_channel;
    };


    template<typename T>
    class GadgetPropertyLimits {
    public:
        virtual bool within_limits(T &v) = 0;

        virtual const char *limits_description() = 0;
        virtual ~GadgetPropertyLimits() = default;
    };

    template<typename T>
    class GadgetPropertyLimitsNoLimits
            : public GadgetPropertyLimits<T> {
    public:
        virtual bool within_limits(T &v) {
            return true;
        }

        virtual const char *limits_description() {
            return "";
        }
    };

    template<typename T>
    class GadgetPropertyLimitsEnumeration
            : public GadgetPropertyLimits<T> {
    public:
        GadgetPropertyLimitsEnumeration(std::initializer_list<T> valid_vals) {
            valid_vals_.insert(valid_vals_.end(), valid_vals.begin(), valid_vals.end());
        }

        virtual bool within_limits(T &v) {
            typename std::vector<T>::iterator it;
            it = find(valid_vals_.begin(), valid_vals_.end(), v);
            if (it != valid_vals_.end()) return true;
            return false;
        }

        virtual const char *limits_description() {
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

    template<typename T>
    class GadgetPropertyLimitsRange
            : public GadgetPropertyLimits<T> {
    public:
        GadgetPropertyLimitsRange(T min_val, T max_val)
                : min_(min_val), max_(max_val) {
        }

        virtual bool within_limits(T &v) {
            return ((v >= min_) && (v <= max_));
        }

        virtual const char *limits_description() {
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


    template<typename T>
    inline void GadgetProperty_extract_value(const char *val, T &tmp) {
        std::stringstream(val) >> std::boolalpha >> tmp;
    }

    template<>
    inline void GadgetProperty_extract_value(const char *val, std::string &tmp) {
        tmp = std::string(val);
    }

    template<typename T, typename L>
    class GadgetProperty
            : public GadgetPropertyBase {
    public:
        GadgetProperty(const char *name, const char *type_string, const char *description,
                       Gadget *g, T default_value, L limits)
                : GadgetPropertyBase(name, type_string, description), g_(g), limits_(limits) {
            g_->register_property(this);
            this->value(default_value);
        }

        T value() const {
            if (is_reference_) {
                std::string val = this->g_->get_string_value(this->name());
                std::stringstream(val) >> std::boolalpha >> value_;
            }
            return value_;
        }

        void value(T v) {
            value_ = v;
            std::stringstream strstream;
            strstream << std::boolalpha << v;
            strstream >> str_value_;
            is_reference_ = false;
            if (!limits_.within_limits(v)) {
                GERROR("Property: %s, value: %s, limits:%s\n", this->name(), str_value_.c_str(),
                       this->limits_.limits_description());
                throw std::runtime_error("Value assigned outside limit range");
            }
        }

        virtual void string_value(const char *val) {
            GadgetPropertyBase::string_value(val);

            if (!is_reference_) {
                T tmp;
//          std::stringstream(val) >> std::boolalpha >> tmp;
                GadgetProperty_extract_value(val, tmp);
                this->value(tmp);
            }
        }

        bool operator==(const T &v) const {
            return this->value() == v;
        }

        operator T() const { return this->value(); }

        virtual const char *limits_description() {
            return limits_.limits_description();
        }

    protected:
        mutable T value_;
        Gadget *g_;
        L limits_;
    };

    template<typename T, typename L>
    class GadgetProperty<std::vector<T>, L> : public GadgetPropertyBase {
    public:
        GadgetProperty(const char *name, const char *type_string, const char *description,
                       Gadget *g, std::initializer_list<T> default_value, L limits)
                : GadgetPropertyBase(name, type_string, description), g_(g), limits_(limits) {
            g_->register_property(this);
            this->value(default_value);
        }

        std::vector<T> value() const {
            if (is_reference_) {
                std::string val = this->g_->get_string_value(this->name());
                std::stringstream ss(val);
                T chunk;
                values_ = std::vector<T>();
                while (ss >> std::boolalpha >> chunk)
                    values_.push_back(chunk);
            }
            return values_;
        }

        operator std::vector<T>() const { return this->value(); }

        void value(std::vector<T> v) {
            values_ = v;
            std::stringstream strstream;
            for (T val : values_)
                strstream << std::boolalpha << val << " ";
            strstream >> str_value_;
            is_reference_ = false;
            if (!limits_.within_limits(v)) {
                GERROR("Property: %s, value: %s, limits:%s\n", this->name(), str_value_.c_str(),
                       this->limits_.limits_description());
                throw std::runtime_error("Value assigned outside limit range");
            }
        }

        virtual void string_value(const char *val) {
            GadgetPropertyBase::string_value(val);

            if (!is_reference_) {
                std::vector<T> tmp;
                T chunk;
                std::stringstream ss(val);
                while (ss >> std::boolalpha >> chunk)
                    tmp.push_back(chunk);
                this->value(tmp);
            }
        }


        bool operator==(const std::vector<T> &v) const {
            return this->value() == v;

        }

        virtual const char *limits_description() {
            return limits_.limits_description();
        }

    protected:
        mutable std::vector<T> values_;
        L limits_;
        Gadget *g_;
    };

#define GADGET_PROPERTY(varname, vartype, description, defaultvalue) GadgetProperty<vartype, GadgetPropertyLimitsNoLimits<vartype> > varname{#varname,#vartype, description, this, defaultvalue, GadgetPropertyLimitsNoLimits<vartype>()}
#define GADGET_PROPERTY_LIMITS(varname, vartype, description, defaultvalue, limitstype, ...) GadgetProperty<vartype, limitstype<vartype> > varname{#varname,#vartype, description, this, defaultvalue, limitstype<vartype>{ __VA_ARGS__ }}

    template<class P1>
    class Gadget1 : public Gadget {

    protected:
        int process(GadgetContainerMessageBase *mb) final {
            GadgetContainerMessage<P1> *m = AsContainerMessage<P1>(mb);

            if (!m) {
                return (this->next()->putq(mb));
            }

            return this->process(m);
        }

        virtual int process(GadgetContainerMessage<P1> *m) = 0;
    };

    template<class P1, class P2>
    class Gadget2 : public Gadget {

    protected:
        int process(GadgetContainerMessageBase *mb) final {
            GadgetContainerMessage<P1> *m1 = AsContainerMessage<P1>(mb);

            GadgetContainerMessage<P2> *m2 = 0;
            if (m1) {
                m2 = AsContainerMessage<P2>(m1->cont());
            }

            if (!m1 || !m2) {
                return (this->next()->putq(mb));
            }

            return this->process(m1, m2);
        }

        virtual int process(GadgetContainerMessage<P1> *m1, GadgetContainerMessage<P2> *m2) = 0;
    };


    template<class P1, class P2, class P3>
    class Gadget3 : public Gadget {

    protected:
        int process(GadgetContainerMessageBase *mb) final {

            GadgetContainerMessage<P1> *m1 = AsContainerMessage<P1>(mb);

            GadgetContainerMessage<P2> *m2 = 0;
            if (m1) {
                m2 = AsContainerMessage<P2>(m1->cont());
            }

            GadgetContainerMessage<P3> *m3 = 0;
            if (m2) {
                m3 = AsContainerMessage<P3>(m2->cont());
            }

            if (!m1 || !m2 || !m3) {
                return (this->next()->putq(mb));

            }

            return this->process(m1, m2, m3);
        }

        virtual int
        process(GadgetContainerMessage<P1> *m1, GadgetContainerMessage<P2> *m2, GadgetContainerMessage<P3> *m3) = 0;
    };

    template<class P1, class P2>
    class Gadget1Of2 : public Gadget {
    protected:
        int process(GadgetContainerMessageBase *mb) final {
            GadgetContainerMessage<P1> *m1 = nullptr;
            GadgetContainerMessage<P2> *m2 = nullptr;

            if ((m1 = AsContainerMessage<P1>(mb))) {
                return this->process(m1);
            } else if ((m2 = AsContainerMessage<P2>(mb))) {
                return this->process(m2);
            } else {
                return (this->next()->putq(mb));
            }
        }

        virtual int process(GadgetContainerMessage<P1> *m1) = 0;

        virtual int process(GadgetContainerMessage<P2> *m1) = 0;
    };


    class LegacyGadgetNode : public Core::Node {
    public:
        LegacyGadgetNode(
                std::unique_ptr<Gadget> gadget_ptr,
                const Core::Context& context,
                const std::unordered_map<std::string, std::string> &props
        );

        void process(Core::GenericInputChannel& in,
                     Core::OutputChannel& out) override;

    private:

        std::unique_ptr<Gadget> gadget;
    };
}


#define GADGET_FACTORY_DECLARE(GadgetClass)                                         \
std::unique_ptr<Gadgetron::Core::Node> legacy_gadget_factory_##GadgetClass(         \
        const Gadgetron::Core::Context &context,                                    \
        const std::string &name,                                                    \
        const std::unordered_map<std::string, std::string> &props                   \
) {                                                                                 \
    auto gadget = std::make_unique<GadgetClass>();                                  \
    gadget->set_name(name);                                                         \
    return std::make_unique<Gadgetron::LegacyGadgetNode>(                           \
            std::move(gadget),                                                      \
            context,                                                                \
            props                                                                   \
    );                                                                              \
}                                                                                   \
                                                                                    \
BOOST_DLL_ALIAS(                                                                    \
        legacy_gadget_factory_##GadgetClass,                                        \
        gadget_factory_export_##GadgetClass                                         \
)                                                                                   \

