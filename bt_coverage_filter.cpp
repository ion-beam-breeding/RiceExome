/**
 * Copyright (c) 2015-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file bt_coverage_filter.cpp
 * @author Hiroyuki Ichida <histfd@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 (GPL-2.0)
 * as published by the Free Software Foundation, Inc.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "histd.h"
#include "headerline.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <unistd.h>

// keywords
#define KEY_CHR         "#CHROM"
#define KEY_START       "ChrStart"
#define KEY_END         "ChrEnd"
#define KEY_FRAC        ".frac"
#define DEFAULT_COVERAGE_THRESHOLD    1.0
//-----------------------------------------------------------------------------
bool is_all_present(const hi::StringArray &elements, \
        const hi::SzArray &column_pos, float *coverage_threshold){

    double frac;
    for(hi::SzArray::const_iterator pos=column_pos.begin(); pos!=column_pos.end(); pos++){
        if(1.0 - std::atof(elements[*pos].c_str()) > 1e-10)
            return false;
    }

    return true;
}
//-----------------------------------------------------------------------------
bool is_line_specific(const hi::StringArray &elements, \
        const hi::SzArray &column_pos, float *coverage_threshold, size_t *specific_column){

    bool isLineSet=false;
    float frac;
    for(size_t pos=0; pos!=column_pos.size(); pos++){
        frac = std::atof(elements[column_pos.at(pos)].c_str());
        if(frac < *coverage_threshold && false==isLineSet){
            *specific_column = pos;
            isLineSet = true;
        }
        else if(1.0 - frac > 1e-10)
            return false;
    }

    if(isLineSet)
        return true;
    return false;
}
//-----------------------------------------------------------------------------
bool find_column_position(const hi::StringArray &header_elements, \
        const char *keyword, size_t *column_pos){

    for(hi::StringArray::const_iterator element=header_elements.begin(); \
            element!=header_elements.end(); element++){

        if(std::string::npos != element->find(keyword)){
            *column_pos = std::distance(header_elements.begin(), element);
            return true;
        }
    }

    return false;
}
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
int nArg=1;
size_t szKeyFrac=std::strlen(KEY_FRAC);

    if(argc <= nArg){
        std::cerr << USAGE_STRING << argv[0] \
                << " (-n program) (-t threshold) input_file" << ENDL;
        exit(EXIT_FAILURE);
    }

    // parse arg
    char option;
    std::string program_name="Program";
    float coverage_threshold=DEFAULT_COVERAGE_THRESHOLD;
    while ((option = getopt(argc, argv, "n:t:")) != -1){
        switch (option){
            case 'n':
                program_name = optarg;
                break;
            case 't':
                coverage_threshold = std::atof(optarg);
                break;
        }
    }
    const char *input_fn = argv[argc-1];

    // open the input file
    std::ifstream file(input_fn, std::ios::in);
    if(file.fail()){
        std::cerr << ERROR_STRING << "the input file (" \
                << input_fn << ") open failed." << ENDL;
        exit(EXIT_FAILURE);
    }

    // find header line, parse & output
    std::string header;
    bool is_header_found=false;
    while(std::getline(file, header)){
        if("##" == header.substr(0, 2)){
            std::cout << header << ENDL;
            continue;
        }
        else if('#'==header[0] && std::string::npos != header.find(KEY_CHR)){
            is_header_found = true;
            break;
        }
    }
    if(! is_header_found){
        std::cerr << ERROR_STRING << "invalid file structure. Header line (" \
                << KEY_CHR << ") can't be found." << ENDL;
        exit(EXIT_FAILURE);
    }

    // header
    std::string cmdstr = generate_cmd_string(argc, argv);
    std::stringstream inputstr;
    inputstr << "input_fn=" << input_fn;
    write_basic_header(__FILE__, __DATE__, __TIME__, cmdstr.c_str(), inputstr.str().c_str(), std::cout);
    std::cout << "#Program\tLine\t" << header.substr(header.find_first_not_of('#')) << "\tLink"<< ENDL;
    hi::StringArray header_elements;
    hi::split(header_elements, header, '\t');

    // find chr, start, and end columns in the header
    size_t chr_column_pos, start_column_pos, end_column_pos;
    if(! find_column_position(header_elements, KEY_CHR, &chr_column_pos) ||
        ! find_column_position(header_elements, KEY_START, &start_column_pos) ||
        ! find_column_position(header_elements, KEY_END, &end_column_pos) ){
        std::cerr << ERROR_STRING \
                << "invalid header structure. Can't find the keywords." << ENDL;
        exit(EXIT_FAILURE);
    }

    // determine the positions of ".frac" columns in the header
    hi::SzArray column_pos;
    hi::StringArray strain_names;
    for(hi::StringArray::iterator iter=header_elements.begin(); \
            iter!=header_elements.end(); iter++){

        ssize_t key_pos = iter->find(KEY_FRAC);
        if(std::string::npos!=key_pos && (iter->length()-szKeyFrac)==key_pos){
            size_t pos = std::distance(header_elements.begin(), iter);
            column_pos.push_back(pos);

            std::string name = header_elements[pos].substr(0, header_elements[pos].length()-std::strlen(KEY_FRAC));
            strain_names.push_back(name);
        }
    }

    // process each record
    std::string line, last_chr="", last_start="", last_end="";
    hi::StringArray elements;
    size_t specific_pos;
    while(std::getline(file, line)){
        elements.clear();
        hi::split(elements, line, '\t');

        if(elements[chr_column_pos]==last_chr \
                && elements[start_column_pos]==last_start \
                && elements[end_column_pos]==last_end){

            continue;
        }

        if(is_line_specific(elements, column_pos, &coverage_threshold, &specific_pos))
            std::cout << program_name << '\t' << strain_names[specific_pos] \
            << '\t' << line << '\t' \
            << "=HYPERLINK(\"http://localhost:60151/goto?locus=" \
            << elements[chr_column_pos] << ':' << elements[start_column_pos] \
            << '-' << elements[end_column_pos] \
            << "\", \"link\")" << ENDL;

        last_chr   = elements[chr_column_pos];
        last_start = elements[start_column_pos];
        last_end   = elements[end_column_pos];
    }

    exit(EXIT_SUCCESS);
}
