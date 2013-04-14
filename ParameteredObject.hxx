/*! \file ParameteredObject.hxx
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
#ifndef _ParameteredObject_H
#define _ParameteredObject_H

#include <map>
#include <sstream>
#include "ParameterFile.hxx"
#include "Parameter.hxx"
#include <stdarg.h>

//todo: derive all classes from this one


template <class T>
class ParameteredObject
{
private:
	static std::map<std::string, unsigned int> genericClassNameCount;

	std::string genericName()
	{
		ParameteredObject::genericClassNameCount[className]++;
		std::stringstream str;
		str << className << ParameteredObject::genericClassNameCount[className];
		return str.str();
	}
protected:
	std::string className, instanceName;
	std::vector<AbstractParameter *> params;
	std::map<std::string, ParameteredObject<T>*> in, out;

	//this functions needs to be called by the derived class in order
	//to register all objects which can be used as inout data for
	//this object. this object can then use the data of this object
	//to generate results which can in turn be read from other objects.
	//the argument list MUST end with a 0 entry!
	void inObjects(const char *name, ...)
	{
		in[name] =0;
		va_list arguments;
		va_start(arguments, name);
		const char *nameN;
		while(nameN =va_arg(arguments, const char *))
		{
			in[nameN] =0;
		}
		va_end(arguments); //clean up argument list
	}

/*
	//this function needs to be called by the derived class in order
	//to register all parameters it contains. hence, this class can
	//automatically load and save these parameters to a file without
	//any additional coding. the argument list MUST end with a 0 entry!
	void init(std::string name, ...)
	{
		instanceName =name;
		va_list arguments;
		va_start(arguments, name);
		AbstractParameter *param;
		while((param =va_arg(arguments, AbstractParameter *)) ? true : false)
		{
			params.push_back(param);
		}
		va_end(arguments); //clean up argument list
	}
	void init(std::string name, std::vector<AbstractParameter*> &list)
	{
		instanceName = name;
		params = list;
	}
*/

	void init(std::string name, unsigned int count, AbstractParameter** arr)
	{
		instanceName = name == "" ? genericName() : name;
		params.clear();
		params.resize(count);
		for(unsigned int i =0; i < count; ++i)
			params[i] =arr[i];
	}

	void addInit(unsigned int count, AbstractParameter** arr)
	{
		for(unsigned int i =0; i < count; ++i)
			params.push_back(arr[i]);
	}
public:
	ParameteredObject(std::string className)
		: className(className)
	{
	}
	virtual ~ParameteredObject() {}

	std::string getClassName()
	{
		return className;
	}

	std::string &name()
	{
		return instanceName;
	}

	virtual void save(ParameterFile &pf, std::string name)
	{
		instanceName =name;
		for(unsigned int i =0; i < params.size(); ++i)
			params[i]->save(pf, instanceName);
	}

	virtual void load(ParameterFile &pf, std::string name)
	{
		instanceName =name;
		for(unsigned int i =0; i < params.size(); ++i)
		{
			params[i]->load(pf, instanceName);
		}
	}

	const ParameteredObject<T>* getIn(std::string name)
	{
		if(in.end() != in.find(name))
		{
			return in[name];
		}
		else
			throw getClassName()+" (" + this->name() + ".getIn()): tried to access unkown [in]-parameter '" + name + "'";
	}

	void setIn(ParameteredObject *object, std::string name)
	{
		if(in.end() != in.find(name))
		{
			in[name] =object;
		}
		else
			throw getClassName()+" (" + this->name() + ".setIn()): tried to access unkown [in]-parameter '" + name + "'";
	}
};

template<class T> std::map<std::string, unsigned int> ParameteredObject<T>::genericClassNameCount;

#endif //pragma once
