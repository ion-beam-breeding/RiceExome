/**
 * Copyright (c) 2016-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file histd.h
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

#ifndef EXOME_HISTD_H
#define EXOME_HISTD_H

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

// constants
#define ENDL        '\n'
#define RV_TRUE      1
#define RV_FALSE    -1
#define HI_NAMESPACE    hi

// strings
#define ERROR_STRING	"\x1b[31;1;7mError: \x1b[m"
#define WARNING_STRING	"\x1b[32;1mWarning: \x1b[m"
#define INFO_STRING		"\x1b[36;1mInfo: \x1b[m"
#define USAGE_STRING	"\x1b[36;1mUsage: \x1b[m"

//-----------------------------------------------------------------------------
namespace HI_NAMESPACE{
    // typedefs
    typedef int RETVAL;
    typedef std::vector<std::size_t> SzArray;
	typedef std::vector<std::string> StringArray;
    typedef std::vector<int> IntArray;
	typedef std::vector<double> DoubleArray;

    // functions
    void bad_alloc_exception(const char *function);
    char * FileRead(const char *file);
    bool split(StringArray &result, const std::string line, const char delimiter);
    bool split(StringArray &result, const char *line, const char *delim);
    RETVAL rmspace(char *seq);
    int toupper(char *str);
}	// End of namespace

//-----------------------------------------------------------------------------
struct ChromosomeLocation{
    std::string chr;
    int start, end, length;
    ChromosomeLocation(){
        chr = "";
        start = end = length = -1;
    }
    ChromosomeLocation & operator = (const ChromosomeLocation &b){
        this->chr = b.chr;
        this->start = b.start;
        this->end   = b.end;
        this->length = b.length;
        return *this;
    }
    friend std::ostream & operator << (std::ostream &ost, const ChromosomeLocation &data){
        ost << data.chr << '\t' << data.start << '\t' << data.end << '\t' << data.length;
		return ost;
	}
};
//-----------------------------------------------------------------------------
struct TargetInfo{
    ChromosomeLocation query, subject;
    int dist_start, dist_end;
    TargetInfo(){
        dist_start = dist_end = -1;
    }
    TargetInfo & operator = (const TargetInfo &b){
        this->query = b.query;
        this->subject = b.subject;
        this->dist_start = b.dist_start;
        this->dist_end   = b.dist_end;
        return *this;
    }
    friend std::ostream & operator << (std::ostream &ost, const TargetInfo &data){
        ost << data.query << '\t' << data.subject << '\t' << data.dist_start << '\t' << data.dist_end;
		return ost;
	}
    bool operator < (const TargetInfo &b){
        return (std::min(std::abs(this->dist_start), std::abs(this->dist_end)) < std::min(std::abs(b.dist_start), std::abs(b.dist_end)));
    }
};
typedef std::vector<TargetInfo> TargetInfoArray;
//------------------------------------------------------------------------------

#endif
