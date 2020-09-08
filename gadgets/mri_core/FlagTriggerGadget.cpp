//
// Created by david on 24/08/2020.
//

#include "FlagTriggerGadget.h"
#include <boost/algorithm/string.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include "ChannelAlgorithms.h"
#include "io/from_string.h"
#include "mri_core_acquisition_bucket.h"
#include "mri_core_data.h"

#include <boost/config/warning_disable.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/ast/variant.hpp>

using namespace Gadgetron;
static const std::unordered_map<std::string, Gadgetron::FlagTriggerGadget::TriggerFlags> flag_map = {

    {"first_in_encode_step1", FlagTriggerGadget::TriggerFlags::first_in_encode_step1},
    {"last_in_encode_step1", FlagTriggerGadget::TriggerFlags::last_in_encode_step1},
    {"first_in_encode_step2", FlagTriggerGadget::TriggerFlags::first_in_encode_step2},
    {"last_in_encode_step2", FlagTriggerGadget::TriggerFlags::last_in_encode_step2},
    {"first_in_average", FlagTriggerGadget::TriggerFlags::first_in_average},
    {"last_in_average", FlagTriggerGadget::TriggerFlags::last_in_average},
    {"first_in_slice", FlagTriggerGadget::TriggerFlags::first_in_slice},
    {"last_in_slice", FlagTriggerGadget::TriggerFlags::last_in_slice},
    {"first_in_contrast", FlagTriggerGadget::TriggerFlags::first_in_contrast},
    {"last_in_contrast", FlagTriggerGadget::TriggerFlags::last_in_contrast},
    {"first_in_phase", FlagTriggerGadget::TriggerFlags::first_in_phase},
    {"last_in_phase", FlagTriggerGadget::TriggerFlags::last_in_phase},
    {"first_in_repetition", FlagTriggerGadget::TriggerFlags::first_in_repetition},
    {"last_in_repetition", FlagTriggerGadget::TriggerFlags::last_in_repetition},
    {"first_in_set", FlagTriggerGadget::TriggerFlags::first_in_set},
    {"last_in_set", FlagTriggerGadget::TriggerFlags::last_in_set},
    {"first_in_segment", FlagTriggerGadget::TriggerFlags::first_in_segment},
    {"last_in_segment", FlagTriggerGadget::TriggerFlags::last_in_segment},
    {"is_noise_measurement", FlagTriggerGadget::TriggerFlags::is_noise_measurement},
    {"is_parallel_calibration", FlagTriggerGadget::TriggerFlags::is_parallel_calibration},
    {"is_parallel_calibration_and_imaging", FlagTriggerGadget::TriggerFlags::is_parallel_calibration_and_imaging},
    {"is_reverse", FlagTriggerGadget::TriggerFlags::is_reverse},
    {"is_navigation_data", FlagTriggerGadget::TriggerFlags::is_navigation_data},
    {"is_phasecorr_data", FlagTriggerGadget::TriggerFlags::is_phasecorr_data},
    {"last_in_measurement", FlagTriggerGadget::TriggerFlags::last_in_measurement},
    {"is_hpfeedback_data", FlagTriggerGadget::TriggerFlags::is_hpfeedback_data},
    {"is_dummyscan_data", FlagTriggerGadget::TriggerFlags::is_dummyscan_data},
    {"is_rtfeedback_data", FlagTriggerGadget::TriggerFlags::is_rtfeedback_data},
    {"is_surfacecoilcorrectionscan_data", FlagTriggerGadget::TriggerFlags::is_surfacecoilcorrectionscan_data},

    {"compression1", FlagTriggerGadget::TriggerFlags::compression1},
    {"compression2", FlagTriggerGadget::TriggerFlags::compression2},
    {"compression3", FlagTriggerGadget::TriggerFlags::compression3},
    {"compression4", FlagTriggerGadget::TriggerFlags::compression4},
    {"user1", FlagTriggerGadget::TriggerFlags::user1},
    {"user2", FlagTriggerGadget::TriggerFlags::user2},
    {"user3", FlagTriggerGadget::TriggerFlags::user3},
    {"user4", FlagTriggerGadget::TriggerFlags::user4},
    {"user5", FlagTriggerGadget::TriggerFlags::user5},
    {"user6", FlagTriggerGadget::TriggerFlags::user6},
    {"user7", FlagTriggerGadget::TriggerFlags::user7},
    {"user8", FlagTriggerGadget::TriggerFlags::user8}};



namespace {

namespace ast {

namespace fusion = boost::fusion;
namespace x3 = boost::spirit::x3;

struct AcquisitionFlags_ : x3::symbols<FlagTriggerGadget::TriggerFlags> {
    AcquisitionFlags_() {
        for (auto [key, value] : flag_map)
            this->add(key, value);
    }
};

enum class BinaryOperator { AND, OR };

struct BinaryOperator_ : x3::symbols<BinaryOperator> {
    BinaryOperator_() {
        this->add("&&", BinaryOperator::AND)("||", BinaryOperator::OR)("AND", BinaryOperator::AND)("OR",
                                                                                                   BinaryOperator::OR);
    }
};

struct BinaryOperation;
struct NegateOperation;

struct Term {
    FlagTriggerGadget::TriggerFlags trigger;
    Term& operator=(FlagTriggerGadget::TriggerFlags flags){ trigger = flags; return *this;}
};

struct Expression : x3::variant < Term,x3::forward_ast<NegateOperation>,
    x3::forward_ast<BinaryOperation>>{
    using base_type::base_type;
    using base_type::operator=;
};

struct NegateOperation {
    int unused;
    Expression expression;
    NegateOperation& operator=(const Expression& exp) { expression = exp; return *this;}
};
struct BinaryOperation {
    Expression left;
    BinaryOperator op;
    Expression right;
};



struct eval {
    bool operator()(Term const & term, uint64_t flags) {
        return (uint64_t(1) >> static_cast<uint64_t>(term.trigger)) & flags;
    }

    bool operator()(NegateOperation const& op, uint64_t flags) { return !(*this)( (Expression const &)op, flags); }

    bool operator()(BinaryOperation const& op, uint64_t flags) {
        bool left = (*this)(op.left, flags);
        bool right = (*this)(op.right, flags);

        switch (op.op) {
        case BinaryOperator::AND:
            return left && right;
        case BinaryOperator::OR:
            return left || right;
        }
        throw std::runtime_error("This really cannot happen.");
    }

    bool operator()(Expression const& expression, uint64_t flags) {
        return boost::apply_visitor([this,flags](auto const & val){return (*this)(val, flags);},expression);
    }
};

} // namespace ast
} // namespace
BOOST_FUSION_ADAPT_STRUCT(ast::NegateOperation, unused,expression);
BOOST_FUSION_ADAPT_STRUCT(ast::BinaryOperation, left, op,right);
BOOST_FUSION_ADAPT_STRUCT(ast::Term, trigger);

namespace boolean_grammer {
namespace x3 = boost::spirit::x3;
using x3::char_;

x3::rule<class expression, ast::Expression> const expression = "expression";
x3::rule<class binop, ast::BinaryOperation> const binop = "binop";
x3::rule<class negateop, ast::NegateOperation> const negateop = "negateop";
x3::rule<class term, ast::Term> const term = "term";

ast::AcquisitionFlags_ acquisitionflags;
ast::BinaryOperator_ binaryoperator;

auto const binop_def = expression >> binaryoperator >> expression;
auto const negateop_def = x3::lit('!') >> x3::attr(1) >>
                          expression;
auto const term_def = acquisitionflags;
auto const expression_def = term | negateop | binop | x3::lit('(') >> expression >> x3::lit(')');

BOOST_SPIRIT_DEFINE(expression, binop, negateop, term);

} // namespace boolean_grammer
namespace {
template <class Iterator> ast::Expression parse_triggers(Iterator first, Iterator last) {
    namespace x3 = boost::spirit::x3;

    ast::Expression expression;
    ast::BinaryOperation bop;
    ast::NegateOperation nop;
    ast::Term term;
    x3::ascii::space_type space;
    bool r = x3::phrase_parse(first, last, boolean_grammer::expression, space, expression);

    if (!r || first != last) {
        throw std::runtime_error("Not passing");
    }
    return expression;
}
}

void Gadgetron::FlagTriggerGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {

    for (const auto& group : Core::Algorithm::buffer(in, this->predicate))
    {
        auto bucket = AcquisitionBucket();
        for (auto acq : group) {
            bucket.add_acquisition(std::move(acq));
        }
        out.push(std::move(bucket));
    }
}

Gadgetron::FlagTriggerGadget::FlagTriggerGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : ChannelGadget(context, props) {
    using namespace ranges;
    auto expression = parse_triggers(trigger_flags.begin(),trigger_flags.end());
    this->predicate = [expression](const Core::Acquisition& acq){
        auto& [head,data,traj] = acq;
        ast::eval eval;
        return eval(expression,head.flags);
    };
}

namespace Gadgetron {
GADGETRON_GADGET_EXPORT(FlagTriggerGadget);
}