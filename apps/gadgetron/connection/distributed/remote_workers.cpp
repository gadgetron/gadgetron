//
// Created by dchansen on 2/5/19.
//

#include "remote_workers.h"
#include "RemoteChannel.h"
#include <boost/asio/io_service.hpp>
#include <boost/process.hpp>
#include <boost/process/async.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi.hpp>
#include <cstdlib>
#include <string>

BOOST_FUSION_ADAPT_STRUCT(
    Gadgetron::Server::Distributed::Address,
    (std::string, ip)(std::string, port))

namespace {
using Gadgetron::Server::Distributed::Address;
using namespace Gadgetron::Core;

template <typename Iterator>
optional<std::vector<Address>> parse_remote_workers(Iterator first, Iterator last)
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    using ascii::space;
    using qi::_1;
    using qi::char_;
    using qi::double_;
    using qi::lexeme;
    using qi::phrase_parse;

    qi::rule<Iterator, std::string(), ascii::space_type> ipv6 = '[' >> lexeme[+(char_ - ']')] >> ']';
    qi::rule<Iterator, Address(), ascii::space_type> address_rule = '"' >> (lexeme[+(char_ - (lexeme[':'] | lexeme[']'] | lexeme['[']))] | ipv6) >> ':' >> lexeme[+(char_ - '"')] >> '"';
    auto result = std::vector<Address> {};
    bool r = phrase_parse(first, last,
        (
            '[' >> address_rule % ',' >> ']'
            ),
        space,
        result);

    if (first != last) // fail if we did not get a full match
        return none;
    return result;
}

}

std::vector<Gadgetron::Server::Distributed::Address> Gadgetron::Server::Distributed::get_remote_workers()
{
    namespace bp = boost::process;

    std::string remote_worker_command = std::getenv("GADGETRON_REMOTE_WORKER_COMMAND");

    boost::asio::io_service ios;
    std::future<std::string> output;
    bp::system(remote_worker_command, bp::std_out > output, bp::std_err > bp::null, ios);

    auto worker_list = output.get();

    auto maybe_result = parse_remote_workers(worker_list.begin(),worker_list.end());
    if (!maybe_result){
        std::stringstream sstream;
        sstream << "Could not parse remote worker list: " << worker_list;
      throw std::runtime_error(sstream.str());
    }
    return std::move(*maybe_result);
}
