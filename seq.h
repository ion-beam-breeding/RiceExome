/**
 * Copyright (c) 2005-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file seq.h
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

#ifndef EXOME_SEQ_H

#include "histd.h"
#include <new>
#include <string>
#include <cstring>
#include <vector>

namespace HI_NAMESPACE{

	class CSeq{
	public:
		RETVAL open(const char *file);
		char * read(void);
		bool eof(void);
		void close(void);
		CSeq();
		CSeq(const char *file);
		~CSeq();

	private:
		char *fasta_read(void);
		void clear(void);
		char *WFile;
		char *NextStart;
		char *CurrentTitle;
		char *CurrentFn;
		size_t NextPos;
		bool IsEof;
		bool WFile_use;
		RETVAL Warning;
	};

	RETVAL clean_seq(char *buffer);
}	// End of namespace

#define EXOME_SEQ_H
#endif
