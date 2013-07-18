// IceGrain: Extraction and parameterization of grain boundary networks of ice
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
//   display the following acknowledgement: ``This product includes the IceGrain
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
#ifndef _Parameter_H
#define _Parameter_H

#include "ParameterFile.hxx"

//Includes for hdf5
#include <cgp/cgp_hdf5.hxx>


class AbstractParameter
{
protected:
	std::string description;
	std::string name;

	AbstractParameter()
	{
	}
public:
	virtual ~AbstractParameter() {}

	void assign(std::string description, std::string name)
	{
		this->description =description;
		this->name =name;
	}

	std::string getDescripton()
	{
		return description;
	}

	std::string getName()
	{
		return name;
	}

	virtual AbstractParameter& operator=(const AbstractParameter &B)
	{
		description =B.description;
		name =B.name;
		return *this;
	}

	virtual void save(ParameterFile &pf, std::string name) =0;
	virtual void load(ParameterFile &pf, std::string name) =0;
};

template <class T>
class Parameter : public AbstractParameter
{
private:
	T defaultValue;
	T value;
public:
	Parameter()
		: AbstractParameter()
	{
	}
	virtual ~Parameter() {}

	void assign(std::string description, std::string name, T defaultValue)
	{
		AbstractParameter::assign(description, name);
		this->defaultValue =defaultValue;
		value =defaultValue;
	}

	virtual Parameter<T>& operator=(const T &B)
	{
		value =B;
		return *this;
	}

	virtual Parameter<T>& operator=(const Parameter<T> &B)
	{
		AbstractParameter::operator=(B);
		defaultValue =B.defaultValue;
		value =B.value;
		return *this;
	}

	operator T() const
	{
		return value;
	}

	const T &operator()() const
	{
		return value;
	}

	virtual void save(ParameterFile &pf, std::string name)
	{
		pf.set<T>(name+"."+this->name, value);
	}

	virtual void load(ParameterFile &pf, std::string name)
	{
		value =pf.get<T>(name+"."+this->name, defaultValue);
	}
    static void writeParamToHDF5(std::string paramName, int paramStdVal, hid_t & group_id, ParameterFile * paramFile)
    {
        Parameter<int> param;
        param.assign("", paramName, paramStdVal);
        param.load(*paramFile, "config");
        
        int paramArr[1];
        paramArr[0] = param;

        hsize_t dims[2];
        dims[0] = 1;
        dims[1] = 1;

        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        hid_t dataset_id = H5Dcreate1(group_id, paramName.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, paramArr);
    }
    static void writeParamToHDF5(std::string paramName, float paramStdVal, hid_t & group_id, ParameterFile * paramFile)
    {
        Parameter<float> param;
        param.assign("", paramName, paramStdVal);
        param.load(*paramFile, "config");
        
        float paramArr[1];
        paramArr[0] = param;

        hsize_t dims[2];
        dims[0] = 1;
        dims[1] = 1;

        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        hid_t dataset_id = H5Dcreate1(group_id, paramName.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, paramArr);
    }
    static void writeParamToHDF5(std::string paramName, std::string paramStdVal, hid_t & group_id, ParameterFile * paramFile)
    {
        Parameter<std::string> param;
        param.assign("", paramName, paramStdVal);
        param.load(*paramFile, "config");
        
        std::string paramStr = param;
        const char * paramChar = paramStr.c_str();

        hsize_t dims[2];
        dims[0] = 1;
        dims[1] = 1;

        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        hid_t datatype = H5Tcopy (H5T_C_S1);
        H5Tset_size (datatype, H5T_VARIABLE);

        hid_t dataset_id = H5Dcreate1(group_id, paramName.c_str(), datatype, dataspace_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &paramChar);
    }
};

template <class T>
class ParameterList : public AbstractParameter
{
private:
	std::vector<T> value;
	std::string defaultValue;
public:
	ParameterList()
		: AbstractParameter()
	{
	}
	virtual ~ParameterList() {}

	void assign(std::string description, std::string name, std::string defaultValue)
	{
		AbstractParameter::assign(description, name);
		this->defaultValue =defaultValue;
		ParameterFile pf;
		value =pf.getList<T>("", defaultValue);
	}

	virtual ParameterList<T>& operator=(const std::vector<T> &B)
	{
		value =B;
		return *this;
	}

	virtual ParameterList<T>& operator=(const ParameterList<T> &B)
	{
		AbstractParameter::operator=(B);
		defaultValue =B.defaultValue;
		value =B.value;
		return *this;
	}

	operator std::vector<T>() const
	{
		return value;
	}

	std::vector<T> &operator()()
	{
		return value;
	}

	virtual void save(ParameterFile &pf, std::string name)
	{
		pf.set<T>(name+"."+this->name, value);
	}

	virtual void load(ParameterFile &pf, std::string name)
	{
		value =pf.getList<T>(name+"."+this->name, defaultValue);
	}
};

#endif //pragma once
