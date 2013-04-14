/*! \file SplitStream.hxx
 */ 
// CIS: Random forest based segmentation project
//
// Copyright (c) 2013 Tobias Binder.
// 
// This software was developed at the University of Heidelberg by
// Tobias Binder, Bjoern Andres, Thorsten Beier and Arthur Kuehlwein.
// Enquiries shall be directed to tobias.binder@iwr.uni-heidelberg.de.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// - Redistributions in binary form must reproduce the above copyright notice, 
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// - All advertising materials mentioning features or use of this software must 
//   display the following acknowledgement: ``This product includes the CIS
//   package developed by Tobias Binder and others''.
// - The name of the author must not be used to endorse or promote products 
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef _SplitStream_H
#define _SplitStream_H

#include <iostream>
#include <vector>
#include <string>
#ifdef CHARON_PETSC
#include "petsc.h"
#endif

class SplitStream
{
public:
	std::vector<std::ostream*> streams;
	SplitStream(std::vector<std::ostream*> &streams) : streams(streams) {}
	SplitStream(std::ostream &stream1, std::ostream &stream2) : streams(2) { streams[0] =&stream1, streams[1] =&stream2; }
	SplitStream(std::ostream &stream)  : streams(1) { streams[0] =&stream; }
	SplitStream() {}

	void assign(std::ostream &stream)  { streams.resize(1); streams[0] =&stream; }
	void assign(std::ostream &stream1, std::ostream &stream2)  { streams.resize(2); streams[0] =&stream1, streams[1] =&stream2; }
	void assign(std::vector<std::ostream*> &streams)  { this->streams =streams; }

	bool isZeroRank()
	{
#ifdef CHARON_PETSC		
		int initialized =0;
		MPI_Initialized(&initialized);
		if(!initialized) //just output everything if MPI is not yet initialized
			return true;

		PetscMPIInt rank;		
		MPI_Comm_rank(PETSC_COMM_WORLD,&rank); //WARNING: this line crashes, if petsc was not initialized!
		return !rank;
#else
		return true;
#endif
	}

	~SplitStream(void) {}
};

extern SplitStream sout; //needs to be instantiated in some lib or cpp-file!

inline SplitStream& operator<<( SplitStream& strm, bool &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, int &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, long &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, float &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, double &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, const double &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, char *w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, const char *w )
{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, std::string &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, const std::string &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, char &w )
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, std::ofstream& (*w)(std::ofstream&))
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, std::ostream& (*w)(std::ostream&))
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, std::ios& (*w)(std::ios&))
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

inline SplitStream& operator<<( SplitStream& strm, std::ios_base& (*w)(std::ios_base&))
	{ if(strm.isZeroRank()) { for(unsigned int i =0; i < strm.streams.size(); ++i) (*strm.streams[i]) << w; } return strm; }

#endif //pragma once
