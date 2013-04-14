/*! \file StringTool.hxx
 *  \brief String manipulation.
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
#ifndef _StringTool_H
#define _StringTool_H

#include <string>
#include <sstream>
#include <algorithm>

class StringTool
{
public:
	StringTool(void) {}
	~StringTool(void) {}

	static std::string trimRight(const std::string &s, const std::string &t = " \t\r\n")
	{ 
		std::string d(s); 
		std::string::size_type i (d.find_last_not_of(t));
		if (i == std::string::npos)
			return "";
		else
		return d.erase (d.find_last_not_of (t) + 1) ; 
	}

	static std::string trimLeft(const std::string &s, const std::string &t = " \t\r\n") 
	{ 
		std::string d(s); 
		return d.erase (0, s.find_first_not_of (t)) ; 
	}

	static std::string trim(const std::string &s, const std::string &t = " \t\r\n")
	{ 
		std::string d(s); 
		return trimLeft (trimRight(d, t), t) ; 
	}

	static void explode(std::string str, char delimiter, std::vector<std::string> &result)
	{
		std::stringstream strm;
		char p[1024];
		strm << str;
		while(!strm.eof())
		{
			strm.getline(p, 1024, delimiter);
			result.push_back(p);
		}
	}
};

#endif //pragma once
