/*! \file path_functions.h
 * \brief Filepath and filename extraction.
 *
 * This header file contains helper functions for filepath and filename extraction.
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
#pragma once
#include <string>


//a little helper funktion, an example should explain what the funktion does
//if we get an image path with a filename path_a/path_b/path_c/image.jpg we want to remove the "image.jpg"
//to get only the path of the image
std::string get_path(std::string path_and_filename)
{
    int size=path_and_filename.size();
    int position=0;
    bool found=false;

    while(size>=0&&found==false)
    {
      if (path_and_filename[size]=='/')
         {
             position=size;
             found=true;
         }
         size--;
    }

    path_and_filename.erase(position+1,path_and_filename.size()-position );
    return path_and_filename;

}

//a little helper funktion, an example should explain what the funktion does
//if we get an image path with a filename path_a/path_b/path_c/image.jpg we want to remove the "path_a/path_b/path_c/"
//to get only the name of the image ,in this case "image.jpg"
std::string get_filename(std::string path_and_filename)
{
    int size=path_and_filename.size();
    int position=0;
    bool found=false;

    while(size>=0&&found==false)
    {
      if (path_and_filename[size]=='/')
         {
             position=size;
             found=true;
         }
         size--;
    }

    if(found) path_and_filename.erase(0,position+1);
    return path_and_filename;

}

