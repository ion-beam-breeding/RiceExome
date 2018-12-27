/**
 * Copyright (c) 2015-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file genotype_filter.cpp
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

#include "genotype_filter.h"

// -----------------------------------------------------------------------------
int main(int argc, char *argv[]){
size_t nArg=1;

    if(argc <= nArg){
        std::cerr << USAGE_STRING \
                << argv[0] << " (options) [input_fn]" << ENDL;
        std::cerr \
                << " -n  Program name [" << DEFAULT_PROGRAM_NAME << "]" << ENDL \
                << " -e  Number of heterozygous lines allowed [" \
                << NUM_HETERO_ALLOWED << "]" << ENDL \
                << " -a  Use alt-rate filter [FALSE]" << ENDL \
                << " -f  Max allowed fraction of lines over alt-rate [" \
                << MAX_ALT_LINE_FRAC<< "]" << ENDL \
                << " -d  Use read depth filter [FALSE]" << ENDL \
                << " -r  Min support reads [" \
                << MIN_SUPPORT_READS << "]" << ENDL \
                << " -l  Max null alleles allowed [" \
                << NUM_NULL_ALLOWED << "]" << ENDL \
                << " -m  Multi-allelic mode [FALSE]" << ENDL \
                << " -v  Invert selection [FALSE]" << ENDL \
                << " -x  Exclude line name []" << ENDL;
        exit(EXIT_FAILURE);
    }

    // parse arguments
    char option;
    bool use_alt_rate_filter=false, use_depth_filter=false;
    bool allow_multiple_alleles=false, isInvertSelection=false;
    float maxAltLineFrac=MAX_ALT_RATE;
    int minSupportReads=MIN_SUPPORT_READS, nHeteroAllowed=NUM_HETERO_ALLOWED;
    int nNullAllowed=NUM_NULL_ALLOWED;
    StringSet ignoreNames;
    std::string program_name=DEFAULT_PROGRAM_NAME;
    while ((option = getopt(argc, argv, "n:f:r:e:mvadx:l:")) != -1){
        switch (option){
            case 'n':
                program_name = optarg;
                break;
            case 'f':
                maxAltLineFrac = std::atof(optarg);
                break;
            case 'r':
                minSupportReads = std::atoi(optarg);
                break;
            case 'e':
                nHeteroAllowed = std::atoi(optarg);
                break;
            case 'm':
                allow_multiple_alleles = true;
                break;
            case 'v':
                isInvertSelection = true;
            case 'a':
                use_alt_rate_filter = true;
                break;
            case 'd':
                use_depth_filter = true;
                break;
            case 'x':
                ignoreNames.insert(optarg);
                break;
            case 'l':
                nNullAllowed = std::atoi(optarg);
                break;
        }
    }
    const char *input_fn = argv[optind];

    // open the input file
    std::ifstream file(input_fn, std::ios::in);
    if(file.fail()){
        std::cerr << ERROR_STRING \
                << "the input file (" << input_fn << ") open failed." << ENDL;
        exit(EXIT_FAILURE);
    }

    // find header
    bool isProgDefined=false;
    int posChrom;
    std::string header_line;
    while(std::getline(file, header_line)){
        if("##" == header_line.substr(0, 2)){
            std::cout << header_line << ENDL;
            continue;
        }
        else if('#'==header_line[0]){
            if(KEY_CHROM == header_line.substr(1, std::strlen(KEY_CHROM))){
                posChrom = 0;
                break;
            }
            else if(KEY_PROGRAM == header_line.substr(1, std::strlen(KEY_PROGRAM))){
                isProgDefined = true;
                break;
            }
        }
    }

    // parse header line
    hi::StringArray header_elements;
    hi::split(header_elements, header_line, '\t');
    if(isProgDefined){
        // remove previous decisions
        find_column_position(header_elements, KEY_CHROM, &posChrom);
        header_elements.erase(header_elements.begin(), header_elements.begin()+posChrom);

        // update the header
        std::stringstream headerstr;
        for(hi::StringArray::const_iterator \
                iter=header_elements.begin(); iter!=header_elements.end(); ++iter){

            headerstr << '\t' << *iter;
        }
        header_line = headerstr.str();
    }

    // determine the positions of '_GT' columns
    IntArray gtColumns;
    if(! find_column_positions( \
            header_elements, KEY_GT, ignoreNames, gtColumns)){

        std::cerr << ERROR_STRING << "invalid header structure. \"" \
                << KEY_GT << "\" wasn't found. " << ENDL;
        exit(EXIT_FAILURE);
    }

    // determine the range of '_DP' columns
    IntArray dpColumns;
    if(use_depth_filter \
            && false==find_column_positions( \
                header_elements, KEY_DP, ignoreNames, dpColumns)){

        std::cerr << ERROR_STRING << "invalid header structure. \"" \
                << KEY_DP << "\" wasn't found. " << ENDL \
                << "Read depth filter was inactivated." << ENDL;
        use_depth_filter = false;
    }

    // determine the range of '_AD' columns
    IntArray adColumns;
    if(use_alt_rate_filter \
            && false==find_column_positions( \
                header_elements, KEY_AD, ignoreNames, adColumns)){
        std::cerr << ERROR_STRING << "invalid header structure. \"" \
                << KEY_AD << "\" wasn't found. " << ENDL \
                << "Alt rate filter was inactivated." << ENDL;
        use_alt_rate_filter = false;
    }

    // column -> strain name
    ColumnToName names;
    create_column_to_name(header_elements, gtColumns, KEY_GT, names);

    // header
    std::string cmdstr = generate_cmd_string(argc, argv);
    std::stringstream sstr;
    sstr << "input_fn=" << input_fn;
    write_basic_header(__FILE__, __DATE__, __TIME__, \
            cmdstr.c_str(), sstr.str().c_str(), std::cout);
    std::cout << '#' << KEY_PROGRAM << '\t' << "Line" << '\t' \
            <<  header_line.substr(1) << ENDL;

    // process records
    std::string line;
    hi::StringArray elements;
    int num_alleles=2;
    while(std::getline(file, line)){
        elements.clear();
        hi::split(elements, line, '\t');

        // get the max number of alleles in the records
        if(use_alt_rate_filter)
            num_alleles = get_max_alleles(elements, adColumns);

        // remove decision from previous run
        if(isProgDefined){
            program_name = elements.at(0);
            elements.erase(elements.begin(), elements.begin()+posChrom);
        }

        if(2==num_alleles){
            proc_by_two_alleles_mode( \
                elements, names, program_name, \
                gtColumns, adColumns, &nHeteroAllowed, &nNullAllowed, \
                &use_alt_rate_filter, &maxAltLineFrac, \
                &use_depth_filter, &minSupportReads, &isInvertSelection);
        }
        else if(allow_multiple_alleles && 3<=num_alleles){
            proc_by_multiple_alleles_mode( \
                elements, names, program_name, \
                gtColumns, adColumns, &nHeteroAllowed, &nNullAllowed, \
                &isInvertSelection);
        }
    }

    file.close();
    exit(EXIT_SUCCESS);
}
//------------------------------------------------------------------------------
