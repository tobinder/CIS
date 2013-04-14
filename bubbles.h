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

#include "ParameteredObject.hxx"

/*! \fn extract_bubbles(std::string filepath_orig_image, std::string filepath_preproc_image = "../preprocessing/out/", std::string param_file = "parameters.txt")
 * \brief Extract bubble areas from an image and save the results into another image.
 * \param filepath_orig_image Filepath to the source image
 * \param filepath_preproc_image Filepath to the preprocessed images (*.as and *.de)
 * \param param_file Name of the parameter file
 */
void extract_bubbles(std::string filepath_orig_image, std::string filepath_preproc_image = "../preprocessing/out/", std::string param_file = "parameters.txt")
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

	//Bubble reflection threshold
	Parameter<int> ref_threshold;
	ref_threshold.assign("", "bubble_reflection_threshold", 180);
	ref_threshold.load(paramFile, "config");

	//Grain lower threshold
	Parameter<int> grain_low;
	grain_low.assign("", "bubble_grain_lower", 80);
	grain_low.load(paramFile, "config");
	
	//Grain upper threshold
	Parameter<int> grain_high;
	grain_high.assign("", "bubble_grain_upper", 130);
	grain_high.load(paramFile, "config");

	//Bubble closing radius
	Parameter<int> bubble_rad;
	bubble_rad.assign("", "bubble_closing_radius", 10);
	bubble_rad.load(paramFile, "config");

	//Grain box radius
	Parameter<int> grain_box_rad;
	grain_box_rad.assign("", "bubble_grain_radius", 10);
	grain_box_rad.load(paramFile, "config");

    vigra::FImage original_image;
    original_image.resize(dim_x, dim_y);
    importImage(info_orig, destImage(original_image));

    std::string filepath_pre_de = filepath_preproc_image;
    std::string filepath_pre_ed = filepath_preproc_image;

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
    std::string filepath_bubbl_temp_image = filepath_orig_image.append("t.bmp"); 

    //The bubble image which will be exported later on; and a debugging image
    vigra::BRGBImage bubbles(dim_x, dim_y);
    vigra::BRGBImage bubbles_temp(dim_x, dim_y);
    
    //Images for labeling non bubble areas
	vigra::BasicImage <bool> bubbles_bin(dim_x, dim_y);
	vigra::IImage label_image(dim_x, dim_y);

	//Images for labeling grain boundaries
	vigra::BasicImage <bool> grain_bin(dim_x, dim_y);
	vigra::IImage label_image2(dim_x, dim_y);

    //Define some colours
    //Rose
    vigra::RGBValue<unsigned int> rose;
    rose.setRed(255);
    rose.setBlue(255);
    
    //Green
    vigra::RGBValue<unsigned int> green;
    green.setGreen(255);  

    //Black
    vigra::RGBValue<unsigned int> black;

    //Blue
    vigra::RGBValue<unsigned int> blue;
    blue.setBlue(255);  
    
    //Red
    vigra::RGBValue<unsigned int> red;
    red.setRed(255);  
    
    //Yellow
    vigra::RGBValue<unsigned int> yellow;
    yellow.setRed(255);
    yellow.setGreen(255);          

    //Clear bubble and temporary bubble images
    std::cout << "Clearing image... " << std::endl;
    for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            bubbles(x,y) = 0;
            bubbles_temp(x,y) = 0;
        }
    }     

    //Paint areas below gray value threshold
    std::cout << "Painting areas below gray value threshold..." << std::endl;
    for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            if(original_image(x,y) < threshold)
            {                
                bubbles(x,y) = green;  				
            }            
        }
    }
    
	//Remove grain boundaries	
    std::cout << "Removing grain boundaries..." << std::endl;
    for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            if(bubbles(x,y) == green && (preproc_image_de(x,y) <= grain_high && preproc_image_de(x,y) >= grain_low))
            {                
                bubbles(x,y) = 0;
                bubbles_temp(x,y) = yellow;                    
            }          
        }
    }

	//Paint the binary image (grain boundaries)
	for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            if(bubbles_temp(x,y) == yellow)
            {                
                grain_bin(x,y) = true;                                
            }
            else
			{
				grain_bin(x,y) = false;
			}
        }
    }

	int nr_boundaries = vigra::labelImageWithBackground(vigra::srcImageRange(grain_bin), vigra::destImage(label_image2), false, 0);

	std::cout<<"Number of marked grain boundaries: "<<nr_boundaries<<std::endl;
    std::vector<int>  boundary_size(nr_boundaries,0);	
    
    //Paint the binary image (non-bubble areas)
    for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            if(bubbles(x,y) == green)
            {                
                bubbles_bin(x,y) = false;  				
            }
            else bubbles_bin(x,y) = true;				
        }
    }

	for(int x = 0; x < dim_x; x++)
    {
        for(int y = 0; y < dim_y; y++)
        {
            //Area labeling starts with 1
            if(label_image(x,y) > 0)
            {
                boundary_size[label_image(x,y)-1]++;                
            }
        }
	}    
	
	//Mark areas above the gray value threshold and index those which are coherent
	int nr_areas = vigra::labelImageWithBackground(vigra::srcImageRange(bubbles_bin), vigra::destImage(label_image), false, 0);
	
	std::cout<<"Number of marked areas: "<<nr_areas<<std::endl;
	std::string filepath_output="output.txt";
	
	std::vector<int>  area_size(nr_areas,0);
    std::vector<long> area_x_center(nr_areas,0);
    std::vector<long> area_y_center(nr_areas,0);

    for(int x = 0; x < dim_x; x++)
    {
        for(int y = 0; y < dim_y; y++)
        {
            //Area labeling starts with 1
            if(label_image(x,y) > 0)
            {
                area_size[label_image(x,y)-1]++;
                area_x_center[label_image(x,y)-1] += x;
                area_y_center[label_image(x,y)-1] += y;
            }
        }
	}
    
    //Calculate center of mass positions and average area size, considering only areas which are smaller than the closing box radius. Those areas count double.
    int avg_area_size = 0;
    for(int area = 0; area < nr_areas; area++)
    {
        if(area_size[area] < 4*bubble_rad*bubble_rad) avg_area_size = avg_area_size + 2*area_size[area];
        area_x_center[area] = area_x_center[area]/area_size[area];
        area_y_center[area] = area_y_center[area]/area_size[area];
    }
    avg_area_size = avg_area_size/nr_areas;
    
    //Write to output file (for testing only, CAN BE REMOVED)
    /*std::ofstream output_file("../preprocessing/WHUT.txt");
	cimg_library::CImg<unsigned char> img(dim_x, dim_y, 1, 3);
    img.fill(0);
    const unsigned char purple[] = {255, 0, 255};
    for(int area=0; area<nr_areas; area++)
    {
        output_file << "Area " << area+1 << " with size " << area_size[area] << " at ("<<area_x_center[area]<<","<<area_y_center[area]<<"): " << "\n";
        std::stringstream stringstr;
		stringstr << area+1;
		img.draw_text(area_x_center[area], area_y_center[area], stringstr.str().c_str(), purple); 
    }	
    output_file.close();
	img.save("../preprocessing/WHUT.bmp");*/
	
	//Close very small areas first (area size smaller than the bubble closing radius in the parameter file)
    std::vector<float> avg_gray_val(nr_areas,0.0f);      //Average gray value
    std::vector<int> ignore_pixels(nr_areas,0);          //Pixels above the threshold will be ignored
    std::vector<float> r_perc(nr_areas,0.0f);            //Those pixels' ratio on the area as a whole in percent 
    std::vector<int> min_gray_val(nr_areas,255);         //Minimal gray value
    std::vector<int> max_gray_val(nr_areas,0);           //Maximal gray value
    std::vector<int> diff(nr_areas);                     //Difference between minimal and maximal gray value
    int high = 110;                                      //The threshold. Pixels above will be ignored
	
    for(int y = 0; y < dim_y; y++)
	{
		for(int x = 0; x < dim_x; x++)
		{
            int area_id=label_image(x,y)-1;
            if(area_size[area_id] < 4*bubble_rad*bubble_rad)
			{
                if(preproc_image_de(x,y) < min_gray_val[area_id]) min_gray_val[area_id] = preproc_image_de(x,y);
                if(preproc_image_de(x,y) > max_gray_val[area_id]) max_gray_val[area_id] = preproc_image_de(x,y);          
                if(preproc_image_de(x,y) < high)
                {                        
                    avg_gray_val[area_id] = avg_gray_val[area_id] + (float)preproc_image_de(x,y);
                }
                else if(preproc_image_de(x,y) >= high)
                {                            
                    ignore_pixels[area_id]++;
                }                        
			}
		}
	}

    for(int area_id = 0; area_id < nr_areas; area_id++)
    {
        avg_gray_val[area_id] = (float)(avg_gray_val[area_id]/((float)(area_size[area_id] - ignore_pixels[area_id])));
        r_perc[area_id] = (float)((float)ignore_pixels[area_id]/(float)area_size[area_id])*100;
        diff[area_id] = max_gray_val[area_id] - min_gray_val[area_id];
        
        //Debug output
        //std::cout << "Area " << area_id+1 << " with size " << area_size[area_id] << " has an avg gray value of " << avg_gray_val << " and " << r_perc <<"% of all pixels are reflections" << std::endl;            
        //std::cout << "MIN GRAY VAL = " << min_gray_val << ", MAX GRAY VAL = " << max_gray_val <<  "=> DIFF = " << diff << "AVG_AREA_SIZE = " << avg_area_size << std::endl;
	}

    /* Close the area if at least one of the following conditions is met
     * 1) More than 45% of the pixels in the area are above the threshold
     * 2) The area size is equal to or smaller than the average area size
     * 3) The difference between the maximal and minimal gray values is equal to or higher than 45            
     */

    for(int y = 0; y < dim_y; y++)
    {
        for(int x = 0; x < dim_x; x++)
        {
            int area_id=label_image(x,y)-1;
            if(r_perc[area_id] > (float)45 || area_size[area_id] <= avg_area_size || diff[area_id] >= 45)            
            {
                bubbles_temp(x,y) = rose;
                bubbles(x,y) = green;			
            }
        }
    }
    
    //Export the image
    exportImage(srcImageRange(bubbles), vigra::ImageExportInfo(filepath_bubbl_image.c_str()));
    //These are temporary images and thus are exported for testing purposes only
    /*exportImage(srcImageRange(bubbles_temp), vigra::ImageExportInfo(filepath_bubbl_temp_image.c_str()));
    exportImage(srcImageRange(label_image), vigra::ImageExportInfo("../preprocessing/Blubb1.bmp")); 
    exportImage(srcImageRange(label_image2), vigra::ImageExportInfo("../preprocessing/Blubb.bmp"));*/ 
    
    std::cout << "Done.." << std::endl;              
}
