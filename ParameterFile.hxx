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
#ifndef _ParameterFile_H
#define _ParameterFile_H

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <algorithm>
#include "SplitStream.hxx"
#include "StringTool.hxx"
#include "FileTool.hxx"

class ParameterFile
{
private:
	void _set(std::string parameter, std::string value)
	{
		transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
		if(isSet(parameter))
			setParams["[mod] "+parameter] =value;
		else
		{
			setParams["      "+parameter] =value;
			parameterLines.push_back(parameter);
		}
		params[parameter] =value;
	}

	//this vector keeps a copy of all parameters in order to preserve their ordering in the actual file
	//it is used for saving the file in the same order and appending new parameters at the end of the file
	std::vector<std::string> parameterLines;

	std::map<std::string, std::string> params; //maps parameters to their respective values
	std::map<std::string, std::string> setParams; //stores modified and new parameters in a seperate map (see showSetParams)
	bool noFoundWarnings; //warn if parameter was not found?
	bool convertSlashes; //convert linux to windows slashes and vice versa? (depends on define CHARON_LINUX/CHARON_WINDOWS in StringTool)
	char delimiter; //delimiter for lists of values (default is ';')
public:
	ParameterFile(void) : noFoundWarnings(false), convertSlashes(true), delimiter(';') {}
	ParameterFile(std::string fileName) : noFoundWarnings(false), convertSlashes(true), delimiter(';') {load(fileName);}
	static std::string filepath;
	~ParameterFile(void) {}

	void setDelimiter(char delimiter)
	{
		this->delimiter =delimiter;
	}
	void setConvertSlashes(bool convertSlashes)
	{
		this->convertSlashes =convertSlashes;
	}
	bool isSet(std::string parameter)
	{
		transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
		return (params.end() != params.find(parameter));
	}

	template<class T> T get(std::string parameter, T defaultValue)
	{
		transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
		if(isSet(parameter))
		{
			std::stringstream strm;
			std::string p =params[parameter];
			strm << p;
			T v;
			strm >> v;
			return v;
		}
		else
		{
			if(noFoundWarnings) sout << "setting unkown parameter '" << parameter << "' to '" << defaultValue << "'" << std::endl;
			set<T>(parameter, defaultValue);
			return defaultValue;
		}
	}

	template<class T> std::vector<T> getList(std::string parameter, std::string defaultValue)
	{
		transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
		std::vector<T> result;
		if(isSet(parameter))
		{
			if(params[parameter] != "" && params[parameter] != "none")
			{
				std::vector<std::string> tmp;
				StringTool::explode(params[parameter], delimiter, tmp);

				result.resize(tmp.size());
				for(unsigned int i =0; i < tmp.size(); ++i)
				{
					std::stringstream strm;
					strm << tmp[i];
					T v;
					strm >> v;
					result[i] =v;
				}
			}
			return result;
		}
		else
		{
			if(noFoundWarnings) sout << "setting unkown parameter '" << parameter << "' to '" << defaultValue << "'" << std::endl;
			set(parameter, defaultValue);
			return getList<T>(parameter, defaultValue);
		}
	}

	template<class T> void set(std::string parameter, T value)
	{
		transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
		std::stringstream str;
		str << value;
		std::string v;
		str >> v;
		_set(parameter, v);
	}

	template<class T> void set(std::string parameter, std::vector<T> value)
	{
		std::stringstream str;
		if(value.size())
		{
			str << value[0];
			for(unsigned int i =1; i < value.size(); ++i)
			{
				str << ";" << value[i];
			}
			std::string res;
			str >> res;
			_set(parameter, res);
		}
	}


	std::vector<std::string> getKeyList(std::string beginsWith)
	{
		std::vector<std::string> result;
		std::map<std::string, std::string>::const_iterator i =params.begin();
		while(i != params.end())
		{
			std::string key =i->first;
			if(beginsWith.length() && key.substr(0, beginsWith.length()) == beginsWith)
				result.push_back(key);
			++i;
		}
		return result;
	}

	void showSetParams()
	{
		std::map<std::string, std::string>::const_iterator i =setParams.begin();
		int p =0;
		unsigned int maxSize =0;
		while(i != setParams.end())
		{
			std::string key =i->first;
			if(key.size() > maxSize)
				maxSize =(unsigned int)key.size();
			++i;
		}

		i =setParams.begin();
		while(i != setParams.end())
		{
			std::string key =i->first;
			std::string value =i->second;
			std::string fillUp;
			if((p++%2) != 0)
			{
				fillUp =std::string(maxSize-key.size(), ' ');
				sout << key << " " << fillUp << "  " << value << std::endl;
			}
			else
			{
				fillUp =std::string(maxSize-key.size(), '-');
				sout << key << " " << fillUp << "> " << value << std::endl;
			}
			++i;
		}
	}

	void clear()
	{
		parameterLines.clear();
		params.clear();
		setParams.clear();
	}

	void resetSetParams()
	{
		setParams.clear();
	}

	void setNotFoundWarningsOn(bool noFoundWarnings)
	{
		this->noFoundWarnings =noFoundWarnings;
	}

	bool save(std::string fileName)
	{
		std::ofstream file;
		file.open(fileName.c_str());
		if(file.fail())
		{
			throw "ERROR: Parameter file '" + fileName  + "' could not be saved. IO error.";
		}
		else
		{
			toStream(file);
			file.close();
			return true;
		}
	}

	bool load(std::string fileName)
	{
		std::ifstream file;
		file.open(fileName.c_str());
		if(!file)
		{
			std::cout << "ERROR: Parameter file '" + fileName  + "' could not be opened. IO error." << std::endl;
			return 0;
		}
		else
		{
			clear();
			fromStream(file);
			file.close();
			filepath = fileName;
			return true;
		}
	}

	//save parameters in the same order as they where inserted
	void toStream(std::ostream &strm)
	{
		//std::map<std::string, std::string>::const_iterator i =params.begin();
		//while(i != params.end())
		for(unsigned int i =0; i < parameterLines.size(); ++i)
		{
			//std::string key =i->first;
			//std::string value =i->second;
			std::string key =parameterLines[i];
			if(key != "")
			{
				std::string value =params[key];
				strm << key << "\t\t" << value << std::endl;
			}
			else
				strm << std::endl;
			//++i;
		}
	}

	void fromStream(std::istream &strm)
	{
		while(!strm.eof())
		{
			std::string key, value, line;
			if(strm.peek() == '\n')
				parameterLines.push_back(""); //preserve empty lines
			strm >> key;
			key =StringTool::trim(key);
			if(key != "")
			{
				if(key.substr(0, 1) != "#")
				{
					//sout << "reading " << key << "..." << std::endl;
					line ="";
					bool cont =true;
					do
					{
						getline(strm, value, '\n'); //allow spaces and empty std::strings as parameter value
						if(StringTool::trim(value) != "")
						{
							std::string test =value.substr(value.length()-2, value.length()-1);
							cont =(test == " \\"); //allow to continue line with this char sequence
							value =StringTool::trim(value);
							if(cont && value.size() > 2)
								value =value.substr(0, value.length()-2);
							else if(cont && value.size() <= 2)
								value ="";
							for(unsigned int i =0; i < value.size(); ++i)
							{
								if(value[i] == '#') //allow for comments after input
								{
									if(i && cont && value[i-1] == ' ') //conserve space before #
										value =StringTool::trim(value.substr(0, i))+" ";
									else
										value =StringTool::trim(value.substr(0, i));
									break;
								}
							}
							line +=value;
						}
						else
							cont =false;
					} while(cont);

					//value +=tmp;
					transform(key.begin(), key.end(), key.begin(), (int (*)(int))tolower);
					params[key] =line;
					parameterLines.push_back(key);
				}
				else
				{
					getline(strm, value, '\n'); //ignore whole line
				}
			}
		}
	}
};

template <> inline void ParameterFile::set(std::string parameter, std::string value)
{
	_set(parameter, value);
}

template <> inline bool ParameterFile::get(std::string parameter, bool defaultValue)
{
	transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
	if(isSet(parameter))
	{
		std::stringstream strm;
		std::string p =params[parameter];
		transform(p.begin(), p.end(), p.begin(), (int (*)(int))tolower);
		if(p == "true")
			return true;
		if(p == "false")
			return false;

		strm << p; //check for numbers
		bool v;
		strm >> v;
		return v;
	}
	else
	{
		if(noFoundWarnings) sout << "setting unkown parameter '" << parameter << "' to '" << defaultValue << "'" << std::endl;
		set<bool>(parameter, defaultValue);
		return defaultValue;
	}
}

template <> inline std::string ParameterFile::get(std::string parameter, std::string defaultValue)
{
	transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
	if(isSet(parameter))
	{
		std::string p =params[parameter];
		if(convertSlashes)
			FileTool::slashConvert(p);
		return p;
	}
	else
	{
		if(noFoundWarnings) sout << "setting unkown parameter '" << parameter << "' to '" << defaultValue << "'" << std::endl;
		set<std::string>(parameter, defaultValue);
		return defaultValue;
	}
}

template <> inline std::vector<std::string> ParameterFile::getList(std::string parameter, std::string defaultValue)
{
	transform(parameter.begin(), parameter.end(), parameter.begin(), (int (*)(int))tolower);
	std::vector<std::string> result;
	if(isSet(parameter))
	{
		if(params[parameter] != "" && params[parameter] != "none")
		{
			StringTool::explode(params[parameter], delimiter, result);

			if(convertSlashes)
			{
				for(unsigned int i =0; i < result.size(); ++i)
					FileTool::slashConvert(result[i]);
			}
		}
		return result;
	}
	else
	{
		if(noFoundWarnings) sout << "setting unkown parameter '" << parameter << "' to '" << defaultValue << "'" << std::endl;
		set(parameter, defaultValue);
		return getList<std::string>(parameter, defaultValue);
	}
}

#endif //pragma once
