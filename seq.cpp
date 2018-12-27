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

#include "seq.h"

namespace HI_NAMESPACE{

//-----------------------------------------------------------------------------
RETVAL CSeq::open(const char *file){

	if(this->WFile_use){
		this->WFile_use = false;
		this->IsEof		= false;
		this->NextStart = NULL;
		delete[] this->WFile;
		delete[] this->CurrentTitle;
		delete[] this->CurrentFn;
	}

	this->WFile = FileRead(file);
	if(this->WFile==NULL)
		return RV_TRUE;
	this->WFile_use = true;
	this->NextStart = &WFile[0];

	try{
		std::size_t size=std::strlen(file);
		this->CurrentFn = new char [size+1];
		std::memset(this->CurrentFn, '\0', size+1);
	}
	catch(std::bad_alloc){
		bad_alloc_exception("CSeq::open()");
	}

	return RV_TRUE;
}
//-----------------------------------------------------------------------------
char * CSeq::read(void){
char *seq;

	if(false == this->WFile_use)
		return NULL;

	seq = this->fasta_read();
	clean_seq(seq);

	return seq;
}
//-----------------------------------------------------------------------------
bool CSeq::eof(void){
	return this->IsEof;
}
//-----------------------------------------------------------------------------
void CSeq::close(void){

	if(this->WFile_use==true)
		delete[] this->WFile;

	clear();
}
//-----------------------------------------------------------------------------
CSeq::CSeq(){
	this->clear();
}
//-----------------------------------------------------------------------------
CSeq::CSeq(const char *file){
	this->clear();
	this->open(file);
}
//-----------------------------------------------------------------------------
CSeq::~CSeq(){
	this->close();
}
//-----------------------------------------------------------------------------
void CSeq::clear(void){
	this->NextPos		= 0;
	this->IsEof			= false;
	this->WFile_use		= false;
}
//-----------------------------------------------------------------------------
char * CSeq::fasta_read(void){
char *sect_head, *sect_body, *sect_seq, *seq;
std::size_t size, shift;

	try{

	// ターゲットとなるセットの切り出し
	sect_head = std::strchr(this->NextStart, '>');
	if(sect_head=='\0')
		return(sect_head);
	size = std::strcspn(sect_head+1, ">");
	sect_body = new char [size+10];
	std::memset(sect_body, '\0', size+10);
	std::strncpy(sect_body, sect_head, size);

	// 次のシークエンスがあるかチェックしてフラグをセット
	this->NextStart = sect_head+size;
	if(std::strchr(this->NextStart, '>')==NULL)
		this->IsEof = true;
	else
		this->IsEof = false;

	// タイトル行の抽出
	size = std::strcspn(sect_body, "\r\n");
	if(size>0){
		shift = std::strspn(sect_body, "> ");
		this->CurrentTitle = new char [size+10];
		std::memset(this->CurrentTitle, '\0', size+10);
		std::strncpy(this->CurrentTitle, sect_body+shift, size-shift);
	}
	else{
		size = std::strlen(this->CurrentFn);
		this->CurrentTitle = new char [size+1];
		std::memset(this->CurrentTitle, '\0', size+1);
		std::strcpy(this->CurrentTitle, CurrentFn);
	}

	// シークエンスの抽出と整形
	sect_seq = sect_body+size;
	size = std::strlen(sect_seq);
	if(size<=0)
		return(NULL);

	seq = new char [size+1];
	std::memset(seq, '\0', size+1);
	std::strncpy(seq, sect_seq, size);
	}
	catch(std::bad_alloc){
		return(NULL);
	}

	return(seq);
}
//-----------------------------------------------------------------------------
RETVAL clean_seq(char *buffer){
	RETVAL ret = hi::rmspace(buffer);
	if(RV_TRUE != ret)
		return ret;
	hi::toupper(buffer);
	return ret;
}
//-----------------------------------------------------------------------------

}	// End of namespace
