/**
 * Copyright (c) 2017-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file headerline.h
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

#ifndef EXOME_HEADERLINE_H
#define	EXOME_HEADERLINE_H

#include <cstring>
#include <time.h>

#define SIZEOF_DATESTR  32
#define SIZEOF_TIMESTR  32
// -----------------------------------------------------------------------------
std::string generate_cmd_string(int argc, char *argv[]){

    std::stringstream sstr;
    for(int i=0; i<argc; ++i){
        if(0 < i)
            sstr << ' ' << argv[i];
        else
            sstr << argv[i];
    }
    return sstr.str();
}
// -----------------------------------------------------------------------------
template <typename Ttype>
bool write_basic_header(const char *src_fn, const char *src_date, const char *src_time, \
        const char *cmd_string, const char *custom_string, Ttype &output){

    // get current date & time
    time_t t = time(NULL);
    struct tm *timestamp = localtime(&t);
    int year = timestamp->tm_year;
    if(1900 > year)
        year += 1900;
    char datestr[SIZEOF_DATESTR], timestr[SIZEOF_TIMESTR];
    std::memset(datestr, '\0', SIZEOF_DATESTR);
    std::snprintf(datestr, SIZEOF_DATESTR, \
            "%d/%02d/%02d", year, timestamp->tm_mon+1, timestamp->tm_mday);
    std::memset(timestr, '\0', SIZEOF_TIMESTR);
    std::snprintf(timestr, SIZEOF_TIMESTR, \
            "%02d:%02d:%02d", timestamp->tm_hour, timestamp->tm_min, timestamp->tm_sec);

    // write command info & parameters to the output
    output << "##INFO src=" << src_fn << ";compile_date=" << src_date \
            << ";compile_time=" << src_time << ENDL;
    output << "##INFO run_date=" << datestr << ";run_time=" << timestr << ENDL;
    if(NULL != cmd_string)
        output << "##INFO cmd='" << cmd_string << '\'' << ENDL;
    if(NULL != custom_string)
        output << "##INFO " << custom_string << ENDL;

    return true;
}
// -----------------------------------------------------------------------------

#endif
