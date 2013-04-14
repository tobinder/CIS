/*! \file FileTool.hxx
 * \brief Filepath and directory manipulation and extraction.
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
#ifndef _FileTool_H
#define _FileTool_H

#if defined(sun)         || defined(__sun)      || defined(linux)       || defined(__linux) \
 || defined(__linux__)   || defined(__CYGWIN__) || defined(BSD)         || defined(__FreeBSD__) \
 || defined(__OPENBSD__) || defined(__MACOSX__) || defined(__APPLE__)   || defined(sgi) \
 || defined(__sgi)
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>
	#define MODE ,0711
	#define CHARON_LINUX
#elif defined(_WIN32) || defined(__WIN32__)
//	#include <crtdbg.h>
    #include <direct.h>
    #include <io.h>
    #define MODE
	#define mkdir _mkdir
	#define unlink _unlink
	#define chdir _chdir
	#define getcwd _getcwd
	#define CHARON_WINDOWS
#endif

#include "StringTool.hxx"

class FileTool
{
public:
	FileTool(void) {}
	~FileTool(void) {}

	#ifdef CHARON_LINUX
	static const char slash ='/';
	#elif defined(CHARON_WINDOWS)
	static const char slash ='\\';
	#endif

	//this function creates all directories contained
	//in path provided they do not already exist
	static int makePath(char *path)
	{
		std::string p(path);
		return makePath(p);
	}
	static int makePath(std::string &path)
	{
		slashConvert(path);
		std::vector<std::string> dirs;
		StringTool::explode(path, FileTool::slash, dirs);

		std::string cdir ="";
		for(unsigned int i =0; i < dirs.size()-1; ++i)
		{
			if(i)
				cdir +=slash+dirs[i];
			else
				cdir +=dirs[i];
			makeDir(cdir);
		}
		return makeDir(cdir+slash+dirs[dirs.size()-1]);
	}

	static int makeDir(const std::string &dir)
	{
		return mkdir(dir.c_str() MODE);
	}

	static int changeDir(const std::string &dir)
	{
		return chdir(dir.c_str());
	}

	static std::string getCurrentDir()
	{
		char cwd[2048];
		getcwd((char*)cwd, 2048);
		return std::string(cwd);
	}

	static void slashConvert(std::string &src)
	{
		#ifdef CHARON_LINUX
		char wrongSlash ='\\';
		#elif defined(CHARON_WINDOWS)
		char wrongSlash ='/';
		#endif
		replace(src.begin(), src.end(), wrongSlash, (char) FileTool::slash);
	}

	static bool exists(std::string file)
	{
		std::ifstream test(file.c_str());
		return !test.fail();
	}

	static int remove(std::string file)
	{
		return unlink(file.c_str());
	}

	static int rename(std::string oldFile, std::string newFile)
	{
		return ::rename(oldFile.c_str(), newFile.c_str());
	}
};

#ifdef CHARON_LINUX
const char FileTool::slash;
#endif

#endif //pragma once
