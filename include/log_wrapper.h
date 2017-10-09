/* Copyright 2017 Baytekov Nikita */

#ifndef INCLUDE_LOG_WRAPPER_H_
#define INCLUDE_LOG_WRAPPER_H_

#include <boost/move/utility_core.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>

#define LOG_(x) BOOST_LOG_TRIVIAL(x)

namespace logging = boost::log;
namespace sinks = boost::log::sinks;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;


void init_logging() {
    logging::register_simple_formatter_factory<logging::trivial::severity_level, char>("Severity");
    logging::add_file_log(keywords::file_name = "%Y-%m-%d.log",
                          keywords::auto_flush = true,
                          keywords::format = "(%TimeStamp%)[%Severity%]: %Message%");

    logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::warning);
    logging::add_common_attributes();
}

#endif  // INCLUDE_LOG_WRAPPER_H_
