/* Copyright 2017 Baytekov Nikita */
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "../include/default_types.h"
#include "../include/log_wrapper.h"

#ifndef INCLUDE_DATASETREADER_H_
#define INCLUDE_DATASETREADER_H_

template<typename T>
class DatasetReader {
 public:
    Mat<T> read_csv(const char path[], char sep = ',');

 private:
    Vec<T> read_vals(std::string entry_str, char sep = ',');
};


template<typename T>
Mat<T> DatasetReader<T>::read_csv(const char path[], char sep) {
    std::ifstream file(path);
    std::string entry_val;

    std::getline(file, entry_val);  // get size of Y-axis from this one
    // LOG_(trace) << "First line:\"" << entry_val << "\"";
    Vec<T> vals = read_vals(entry_val, sep);
    // LOG_(trace) << "Read vector: " << vals;
    Mat<T> out_mat(vals);

    // LOG_(trace) << "Out matrix" << out_mat;
    while(file.good()) {
        std::getline(file, entry_val);
        if (entry_val.size() > 0)
            out_mat.hadd(read_vals(entry_val, sep));
        else
            LOG_(warning) << "Empty string was discovered during file read.";
    }

    LOG_(info) << "Dataset \"" << path << "\" was successfully read.";
    return out_mat;
}

template<typename T>
Vec<T> DatasetReader<T>::read_vals(std::string entry_str, char sep) {
    // LOG_(trace) << "Reading values from vector...";
    // LOG_(trace) << "Entry string: \"" << entry_str << "\"";
    Vec<T> out_vec;
    T add_val;

    std::stringstream ss;
    ss.str(entry_str);
    std::string token;
    while (std::getline(ss, token, sep).good()) {
        // LOG_(trace) << "New token:\"" << token << "\"";
        std::stringstream token_ss;
        token_ss.str(token);
        token_ss >> add_val;
        // LOG_(trace) << "Read value: " << add_val;
        // LOG_(trace) << out_vec;
        out_vec.append(add_val);
        // LOG_(trace) << "Extended vec: " << out_vec;
    }

    return out_vec;
}

#endif  // INCLUDE_DATASETREADER_H_
