/**
 * Copyright (c) 2016-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file validatetarget_location.cpp
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
#include "seq.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>

//------------------------------------------------------------------------------
struct ShiftMatchData{
    int shift, nMatch;
};
//------------------------------------------------------------------------------
static bool by_nMatch(const ShiftMatchData &a, const ShiftMatchData &b){
    return(a.nMatch < b.nMatch);
}
typedef std::vector<ShiftMatchData> ShiftMatchDataArray;
//------------------------------------------------------------------------------
bool parse_record(std::string &line, TargetInfo &target, bool *isSubjectNotFound){

    hi::StringArray arr;
    hi::split(arr, line, '\t');
    int start, end;

    if(4 < arr.size()){
        target.query.chr = arr[0];
        start  = std::atoi(arr[1].c_str());
        end    = std::atoi(arr[2].c_str());
        target.query.start  = std::min(start, end);
        target.query.end    = std::max(start, end);;
        target.query.length = std::atoi(arr[3].c_str());
    }

    if(11 > arr.size()  &&  "NOT FOUND"==arr[4])
        *isSubjectNotFound = true;
    else if(11 <= arr.size()){
        target.subject.chr = arr[4];
        start  = std::atoi(arr[5].c_str());
        end    = std::atoi(arr[6].c_str());
        target.subject.start  = std::min(start, end);
        target.subject.end    = std::max(start, end);;
        target.subject.length = std::atoi(arr[7].c_str());
        target.dist_start = std::atoi(arr[8].c_str());
        target.dist_end = std::atoi(arr[9].c_str());
    }
    else{
        std::cerr << ERROR_STRING << __LINE__ << ENDL;
        return false;
    }

    return true;
}
//------------------------------------------------------------------------------
char * extract_seq(const char *ref, ChromosomeLocation &location, size_t *szQuery){

    *szQuery = std::abs(location.end - location.start + 1);
    char *seq = new char [(*szQuery)+1];
    std::memset(seq, '\0', (*szQuery)+1);
    std::strncpy(seq, &ref[location.start], *szQuery);
    return seq;
}
//------------------------------------------------------------------------------
bool is_match(const char *query, const size_t *szQuery, const char *subject, \
    const size_t *szSubject, const size_t *szMaxShiftAllowed, int *bestShift){

    size_t length = std::min(*szSubject, *szQuery);
    int nMatch=0;
    for(size_t shift=0; shift<(*szMaxShiftAllowed); ++shift){
        for(int i=0; i<length-shift; ++i){
            if(query[i]==subject[i+shift])
                ++nMatch;
        }
        if(0.95 < double(nMatch) / double(length-shift)){
            *bestShift = shift;
            return true;
        }
    }
    return false;
}
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){
size_t nArg=4;

	if(argc<=nArg){
		std::cerr << USAGE_STRING << argv[0] \
                << " [input_tab] [old_seq] [new_seq] [output_fn]" << ENDL;
		exit(EXIT_FAILURE);
	}
	const char *inFn=argv[1], *querySeqFn=argv[2];
    const char *subjectSeqFn=argv[3], *outFn=argv[4];

    // input file
    std::ifstream infile(inFn, std::ios::in);
    if(infile.fail()){
        std::cerr << ERROR_STRING << __LINE__ << ENDL;
        exit(EXIT_FAILURE);
    }

    // reference sequences
    hi::CSeq query_seqobject;
    if(RV_TRUE!=query_seqobject.open(querySeqFn)){
        std::cerr << ERROR_STRING << __LINE__ << ENDL;
        exit(EXIT_FAILURE);
    }
    hi::CSeq subject_seqobject;
    if(RV_TRUE!=subject_seqobject.open(subjectSeqFn)){
        std::cerr << ERROR_STRING << __LINE__ << ENDL;
        exit(EXIT_FAILURE);
    }

    // output
    std::ofstream outfile(outFn, std::ios::out);
    if(outfile.fail()){
        std::cerr << ERROR_STRING << __LINE__ << ENDL;
        exit(EXIT_FAILURE);
    }

    const size_t maxShiftAllowed=5;
    char *query_genome, *subject_genome, *query_seq, *subject_seq;
    size_t szQueryGenome, szSubjectGenome, szQuery, szSubject;
    int bestShift, diff;

    std::string line;
    bool isSubjectNotFound, isEstimateMatched, isFirst=true;
    int dist_start, dist_end, diff_length;
    TargetInfo target, last;
    ChromosomeLocation estimate;

    while(std::getline(infile, line)){
        isSubjectNotFound = false;
        if(! parse_record(line, target, &isSubjectNotFound)){
            std::cerr << WARNING_STRING \
                    << "parse_record() failed. line=" << line << ENDL;
            continue;
        }

        if(target.query.chr!=last.query.chr){
            if(! isFirst){
                delete[] query_genome;
                delete[] subject_genome;
            }
            query_genome = query_seqobject.read();
            subject_genome = subject_seqobject.read();
            isFirst = true;
        }

        // not found -- estimate from last record
        if(isSubjectNotFound){
            if(isFirst){
                outfile << target << "\tNOT_FOUND" << ENDL;
                continue;
            }

            estimate = target.query;
            estimate.start += last.dist_start;
            estimate.end   += last.dist_end;
            query_seq = extract_seq(query_genome, target.query, &szQuery);
            subject_seq = extract_seq(subject_genome, estimate, &szSubject);
            if(is_match(query_seq, &szQuery, subject_seq, &szSubject, &maxShiftAllowed, &bestShift)){
                target.subject.start = \
                        std::min(estimate.start + bestShift-1, estimate.end + bestShift-1);
                target.subject.end = \
                        std::max(estimate.start + bestShift-1, estimate.end + bestShift-1);
                target.dist_start = last.dist_start;
                target.dist_end = last.dist_end;
                outfile << target << "\tFILLED_FROM_LAST_SHIFT" << ENDL;
                last = target;
            }
            else{
                ChromosomeLocation empty_record;
                target.subject = empty_record;
                outfile << target << "\tNOT_FOUND" << ENDL;
            }
            delete[] subject_seq;
            delete[] query_seq;
            continue;
        }

        // perfect match -- use as is
        diff = std::abs(target.dist_start - last.dist_start) + std::abs(target.dist_end - last.dist_end);
        if(0==diff || std::abs(last.query.length - last.subject.length)==diff || isFirst){
            outfile << target << ENDL;
            last = target;
            isFirst = false;
            continue;
        }

        // try estimated location from last dist_start/end if the diff is too big
        estimate = target.query;
        estimate.start += last.dist_start;
        estimate.end   += last.dist_end;
        query_seq = extract_seq(query_genome, target.query, &szQuery);
        subject_seq = extract_seq(subject_genome, estimate, &szSubject);
        if(is_match(query_seq, &szQuery, subject_seq, &szSubject, &maxShiftAllowed, &bestShift)){
            target.subject.start = \
                    std::min(estimate.start + bestShift-1, estimate.end + bestShift-1);
            target.subject.end = \
                    std::max(estimate.start + bestShift-1, estimate.end + bestShift-1);
            target.dist_start = last.dist_start;
            target.dist_end   = last.dist_end;
            outfile << target << "\tADJUSTED_BY_LAST_SHIFT" << ENDL;
        }
        else{
            delete[] subject_seq;
            subject_seq = extract_seq(subject_genome, target.subject, &szSubject);
            if(is_match(query_seq, &szQuery, subject_seq, &szSubject, &maxShiftAllowed, &bestShift))
                outfile << target << "\tADJUSTMENT_NOT_SUCCEED" << ENDL;
            else{
                ChromosomeLocation empty_record;
                target.subject = empty_record;
                outfile << target << "\tNOT_FOUND" << ENDL;
                delete[] query_seq;
                delete[] subject_seq;
                continue;
            }
        }
        last = target;
        delete[] query_seq;
        delete[] subject_seq;
    }
    infile.close();
    outfile.close();

	exit(EXIT_SUCCESS);
}
