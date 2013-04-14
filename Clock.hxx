/*! \file Clock.hxx
 * \brief Time measuring.
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
#ifndef _Clock_H
#define _Clock_H

#include <time.h>
#include <string>

class Clock
{
private:
    clock_t start;
	clock_t time;
	bool autoReset;
public:
	Clock(bool autoReset =true) : autoReset(autoReset) {start =clock();}
	~Clock(void) {}

	inline clock_t getDiff()
		{
			time =clock()-start;
			if(autoReset)
				start =clock();
			return time;
		}
	inline std::string getString()
		{
			time =clock()-start;
			if(autoReset)
				start =clock();
			char str[255];
			double seconds =(double)time/(double)CLOCKS_PER_SEC;
			unsigned int minutes =((unsigned int)seconds)/60 % 60;
			unsigned int hours =((unsigned int)seconds)/3600 % 24;
			unsigned int days =((unsigned int)seconds)/86400;
			seconds =fmod(seconds,60);
			if(minutes)
			{
				if(hours)
				{
					if(days)
						sprintf(str, "%03d:%02d:%02d:%02.3f (dd:hh:mm:ss.sss)", days, hours, minutes, seconds);
					else
						sprintf(str, "%02d:%02d:%02.3f (hh:mm:ss.sss)", hours, minutes, seconds);
				}
				else
					sprintf(str, "%02d:%02.3f (mm:ss.sss)", minutes, seconds);
			}
			else
				sprintf(str, "%02.3f seconds", seconds);
			return str;
		}
};
#endif //pragma once
