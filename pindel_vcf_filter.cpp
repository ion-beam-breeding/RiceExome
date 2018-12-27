/**
 * Copyright (c) 2015-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file pindel_vcf_filter.cpp
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

#include"histd.h"
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<string>

//------------------------------------------------------------------------------
#define POS_REF_BASE            3
#define POS_ALT_BASE            4
#define POS_GENOTYPE_START      9
#define POS_INFO_COLUMN         7
#define SVLEN_FILE_THRESHOLD    100000  // SVLEN >= (value) will be put in x_file
//------------------------------------------------------------------------------
int get_svlen(const std::string &infostr){

    ssize_t start = infostr.find("SVLEN=");
    if(std::string::npos == start)
        return -1;
    ssize_t length = infostr.find_first_of(';', start+6) - start -6;

    return std::abs(std::atoi(infostr.substr(start+6, length).c_str()));
}
//------------------------------------------------------------------------------
bool process_vcf(const char *input_fn, const char *o_fn, const char *x_fn){

    std::ifstream file(input_fn, std::ios::in);
    if(file.fail()){
        std::cerr << ERROR_STRING << "input file (" \
                << input_fn << ") open failed." << ENDL;
        return false;
    }

    // output files
    std::ofstream o_file(o_fn, std::ios::out);
    if(o_file.fail()){
        std::cerr << ERROR_STRING << "output file (" \
                << o_fn << ") open failed." << ENDL;
        return false;
    }
    std::ofstream x_file(x_fn, std::ios::out);
    if(x_file.fail()){
        std::cerr << ERROR_STRING << "large file (" \
                << x_fn << ") open failed." << ENDL;
        return false;
    }

    // process file
    std::string line;
    hi::StringArray elements;
    bool isAllRefType;
    int sv_length;
    while(std::getline(file, line)){
        if('#' == line[0]){
            o_file << line << ENDL;
            x_file << line << ENDL;
            continue;
        }

        elements.clear();
        hi::split(elements, line, '\t');

        // Skip if all strain has 0/0 genotype
        isAllRefType = true;
        for(size_t i=POS_GENOTYPE_START; i<elements.size(); ++i){
            if("0/0" != elements.at(i).substr(0,3)){
                isAllRefType = false;
                break;
            }
        }
        if(isAllRefType)
            continue;

        sv_length = get_svlen(elements.at(POS_INFO_COLUMN));
        if(0 > sv_length)
            sv_length = std::max( \
                    elements.at(POS_REF_BASE).length(), \
                    elements.at(POS_ALT_BASE).length());

        if(SVLEN_FILE_THRESHOLD <= sv_length){
            x_file << line << ENDL;
            continue;
        }

        o_file << line << ENDL;
    }

    file.close();
    o_file.close();
    x_file.close();
    return true;
}
//------------------------------------------------------------------------------
inline void print_usage(const char *cmd){
    std::cerr << USAGE_STRING << cmd \
            << " (-t large_threshold) -i [vcf_fn] -o [out_fn] -x [large_fn]" \
            << ENDL;
}
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){
size_t nArg=1;

    if(argc <= nArg){
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // parse arguments
    char option;
    std::string i_fn, o_fn, x_fn;
    bool is_i_set=false, is_o_set=false, is_x_set=false;
    int large_threshold=SVLEN_FILE_THRESHOLD;
    while ((option = getopt(argc, argv, "i:o:x:t:")) != -1){
        switch (option){
            case 'i':
                i_fn = optarg;
                is_i_set = true;
                break;
            case 'o':
                o_fn = optarg;
                is_o_set = true;
                break;
            case 'x':
                x_fn = optarg;
                is_x_set = true;
                break;
            case 't':
                large_threshold = std::atoi(optarg);
                break;
            default:
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
                break;
        }
    }
    if(! (is_i_set && is_o_set && is_x_set))
        print_usage(argv[0]);

    process_vcf(i_fn.c_str(), o_fn.c_str(), x_fn.c_str());

    exit(EXIT_SUCCESS);
}
// -----------------------------------------------------------------------------
