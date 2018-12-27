/**
 * Copyright (c) 2015-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file merge_vc.cpp
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
#include <string>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <map>

#define KEY_PROGRAM "#Program"
#define KEY_LINE    "Line"
#define KEY_CHR	    "CHROM"
#define KEY_START   "ChrStart"
#define KEY_END	    "ChrEnd"
//------------------------------------------------------------------------------
typedef std::map<std::string,std::string> MutationDB;
//------------------------------------------------------------------------------
bool determine_position(hi::StringArray &header_elements, \
        const std::string &keyword, int *column_pos){

    *column_pos = -1;
    ssize_t hit;
    for(size_t pos=0; pos<header_elements.size(); pos++){
    	hit = header_elements.at(pos).find(keyword);
    	if(std::string::npos != hit){
    	    *column_pos = pos;
    	    return true;
    	}
    }
    return false;
}
//------------------------------------------------------------------------------
inline std::string create_locus_id(const std::string &line, \
        const std::string &chr, const std::string &start, const std::string &end){

    std::stringstream sstr;
    sstr << line << '|' << chr << '_' \
    	<< std::setw(8) << std::setfill('0') << start << '-' \
    	<< std::setw(8) << std::setfill('0') << end;
    return sstr.str();
}
//------------------------------------------------------------------------------
bool read_file(const char *filename, MutationDB &mutationDB, std::string &header_line){

    std::ifstream file(filename, std::ios::in);
    if(file.fail()){
    	std::cerr << ERROR_STRING \
                << "an input file (" << filename \
                << ") can't open for reading." << ENDL;
    	return false;
    }

    // header
    std::cout << "##ORIGINAL_FILE: " << filename << ENDL;
    while(std::getline(file, header_line)){
        if('#' == header_line[0] && std::string::npos != header_line.find(KEY_CHR))
            break;
        else if('#' == header_line[0]){
            std::cout << header_line << ENDL;
            continue;
        }
    }
    hi::StringArray elements;
    hi::split(elements, header_line, '\t');

    // find column position
    int pos_program, pos_line, pos_chr, pos_start, pos_end;
    determine_position(elements, KEY_PROGRAM, &pos_program);
    determine_position(elements, KEY_LINE, &pos_line);
    determine_position(elements, KEY_CHR, &pos_chr);
    determine_position(elements, KEY_START, &pos_start);
    determine_position(elements, KEY_END, &pos_end);
    if(0 > pos_program || 0 > pos_line || 0 > pos_chr || 0 > pos_start || 0 > pos_end){
    	std::cerr << ERROR_STRING << "1 or more of critical columns can't find." << ENDL;
    	return false;
    }

    // process records
    std::string line, locus_id, prev_record;
    std::stringstream sstr;
    size_t tabPos;
    MutationDB::iterator iter;
    while(std::getline(file, line)){
    	elements.clear();
            hi::split(elements, line, '\t');

    	locus_id = create_locus_id(elements[pos_line], \
                elements[pos_chr], elements[pos_start], elements[pos_end]);
    	iter = mutationDB.find(locus_id);
    	if(mutationDB.end() == iter)
    	    mutationDB.insert(std::pair<std::string,std::string>(locus_id, line));
        else{
            tabPos = iter->second.find('\t');
            if(std::string::npos == tabPos)
                continue;
            sstr.str("");
            sstr << iter->second.substr(0, tabPos) \
                    << '/' << elements[pos_program] << iter->second.substr(tabPos);
            iter->second = sstr.str();
        }
    }

    file.close();
    return true;
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    // parse arguments
    int nArg=1;
    if(argc<=nArg){
    	std::cerr << USAGE_STRING << argv[0] << " file1 file2 file3..." << ENDL;
    	exit(EXIT_FAILURE);
    }

    // create input filelist
    std::stringstream inputstr;
    std::string header_line;
    MutationDB mutationDB;
    for(size_t i=nArg; i<argc; i++){
    	if(! read_file(argv[i], mutationDB, header_line)){
            std::cerr << WARNING_STRING << "can't open an input file (" \
                    << argv[i] << "). Skipped." << ENDL;
            continue;
        }

        if(nArg == i)
            inputstr << "input_fn=" << argv[i];
        else
            inputstr << ";input_fn_" << i-nArg+1 << "=" << argv[i];
    }

    // header
    std::string cmdstr = generate_cmd_string(argc, argv);
    write_basic_header(__FILE__, __DATE__, __TIME__, \
            cmdstr.c_str(), inputstr.str().c_str(), std::cout);
    std::cout << header_line << ENDL;

    // output results
    for(MutationDB::iterator iter=mutationDB.begin(); iter!=mutationDB.end(); iter++)
	   std::cout << iter->second << ENDL;

    exit(EXIT_SUCCESS);
}
//------------------------------------------------------------------------------
