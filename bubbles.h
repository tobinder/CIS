/*! \file bubbles.h
 * \brief Bubble area extraction.
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
#include "path_functions.h"
#include <sstream>
#include "CImg.h"
#include <fstream>
#include <vigra/localminmax.hxx>

#include "ParameteredObject.hxx"

/*! \fn extract_bubbles(std::string filepath_orig_image, std::string filepath_preproc_image = "../preprocessing/out/", std::string param_file = "parameters.txt")
 * \brief Extract bubble areas from an image and save the results into another image.
 * \param filepath_orig_image Filepath to the source image
 * \param filepath_preproc_image Filepath to the preprocessed images (*.as and *.de)
 * \param param_file Name of the parameter file
 */
void extract_bubbles(std::string filepath_orig_image, std::string filepath_preproc_image, std::string param_file = "parameters.txt")
{
    //Load original image
    std::string filename_orig_image = get_filename(filepath_orig_image);
    std::cout << "Filename of original image: " << filename_orig_image << std::endl;       

    vigra::ImageImportInfo info_orig(filepath_orig_image.c_str());
    int dim_x = info_orig.width();
    int dim_y = info_orig.height();

    //Load parameter file
    ParameterFile paramFile;
    
    if(!paramFile.load(param_file))
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

	//Load parameters
	//Bubble gray value threshold
	Parameter<int> threshold;
	threshold.assign("", "bubble_gray_threshold", 70);
	threshold.load(paramFile, "config");

	//Grain boundary lower threshold
	Parameter<int> boundary_low;
	boundary_low.assign("", "bubble_grain_lower", 80);
	boundary_low.load(paramFile, "config");
	
    vigra::FImage original_image;
    original_image.resize(dim_x, dim_y);
    importImage(info_orig, destImage(original_image));

    std::string filepath_pre_de = filepath_preproc_image;

    //Load preprocessed image (*.de.bmp)
    std::string filename_preproc_image_de = filename_orig_image;
    filename_preproc_image_de.append(".de.bmp");
    std::cout << "Filename of preprocessed image (*.de.bmp): " << filename_preproc_image_de << std::endl;  
    filepath_pre_de.append(filename_preproc_image_de);

    vigra::ImageImportInfo info_pre_de(filepath_pre_de.c_str());
    vigra::FImage preproc_image_de;
    preproc_image_de.resize(dim_x, dim_y);
    importImage(info_pre_de, destImage(preproc_image_de));   

    //Filepaths for bubble and bubble debugging images
    filepath_orig_image.resize(filepath_orig_image.size()-3);
    std::string filepath_bubbl_image = filepath_orig_image.append("b.bmp");

    //The bubble image which will be exported later on; and a debugging image
    vigra::BasicImage<bool> bubbles(dim_x, dim_y);
  
    //Initialize bubble and temporary bubble images
    for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            bubbles(x,y) = false;
        }
    }     

    for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            if(original_image(x,y) < threshold && preproc_image_de(x,y) < boundary_low)
            {                
                bubbles(x,y) = true;  				
            }            
        }
    }

    //fill inside of bubbles
    vigra::extendedLocalMinima(vigra::srcImageRange(bubbles),vigra::destImage(bubbles));

    //Export the image
    exportImage(srcImageRange(bubbles), vigra::ImageExportInfo(filepath_bubbl_image.c_str()));
    
    std::cout << "Bubble image has been written to: " << filepath_bubbl_image << std::endl;              
}
