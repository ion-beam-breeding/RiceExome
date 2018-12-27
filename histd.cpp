/**
 * Copyright (c) 2016-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file histd.cpp
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

namespace HI_NAMESPACE{

//-----------------------------------------------------------------------------
void bad_alloc_exception(const char *function){
	std::cerr << ERROR_STRING << "a fatal error occurred in " << function \
			<< ". Program execution aborted." << ENDL;
	std::exit(EXIT_FAILURE);
}
//-----------------------------------------------------------------------------
char * FileRead(const char *file){

	int fDesc = open(file, O_RDONLY);
	if(0 > fDesc)
		return NULL;

	lseek(fDesc, 0, SEEK_SET);
	long size = lseek(fDesc, 0, SEEK_END);

	char *buf = new char [size+10];
	std::memset(buf, '\0', size+10);

	lseek(fDesc, 0, SEEK_SET);
	read(fDesc, buf, size);
	close(fDesc);

	return buf;
}
//-----------------------------------------------------------------------------
bool split(StringArray &result, const std::string line, const char delimiter){

	std::stringstream sstr(line);
	std::string element;
	while(std::getline(sstr, element, delimiter))
		result.push_back(element);

	return(true);
}
//-----------------------------------------------------------------------------
bool split(StringArray &result, const char *line, const char *delim){

	size_t size = std::strlen(line);

	try{
		char *buf = new char [size+1];
		std::memset(buf, '\0', size+1);
		std::strncpy(buf, line, size);

		char *saveptr;
		char *element = strtok_r(buf, delim, &saveptr);
		while(NULL != element){
			result.push_back(element);
			element = strtok_r(NULL, delim, &saveptr);
		}
		delete[] buf;
	}
	catch(std::bad_alloc){
		return false;
	}

	return true;
}
//-----------------------------------------------------------------------------
RETVAL rmspace(char *seq){

	size_t size = std::strlen(seq);
	if(0 == size)
		return RV_FALSE;

	char *newseq = new char [size+1];
	std::memset(newseq, '\0', size+1);

	size_t newpos=0;
	for(size_t pos=0; pos<size; pos++){
		if(0 != std::isalpha(seq[pos])){
			newseq[newpos] = seq[pos];
			++newpos;
		}
	}

	std::memset(seq, '\0', size);
	std::strncpy(seq, newseq, newpos);
	delete[] newseq;

	return RV_TRUE;
}
//-----------------------------------------------------------------------------
int toupper(char *str){
size_t pos=0;

	char *chr = str;
	while ('\0' != *chr) {
		if('a'<=*chr && 'z'>=*chr) {
			++pos;
			*chr -= 0x20;
		}
		++chr;
    }
    return pos;
}
//-----------------------------------------------------------------------------
}	// End of namespace
