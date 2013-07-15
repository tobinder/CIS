/*! \file boundary_features_gray.h
 * \brief Boundary feature calculation for grayscale images.
 *
 * This header file contains functions for boundary feature calculations.
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
 
#include <math.h>
 
#include "CImg.h"
 
#include "boundary_data_structure.h"
#include "cgp_structure.hxx"
 
//functions that differs not from the color case are defined in boundary_features.h

std::vector<float> compute_single_arc_region_features_gray(int arc_index, std::vector <float> area_mean, std::vector <long> area_size, marray::Marray<unsigned int> const & two_boundings)
{
    //now we have to look in the two boundings which regions are the neigbours
    size_t n0=two_boundings(arc_index,0);
    size_t n1=two_boundings(arc_index,1);
 
    std::vector<float> r_features;
 
    r_features.push_back(fabs(area_mean[n0-1]-area_mean[n1-1]));
    r_features.push_back(area_size[n0-1]+area_size[n1-1]);
 
    return r_features;
}

std::vector<float> compute_single_arc_cross_section_features_gray(int arc_index,
                                                                  std::vector< std::vector<float> > & cross_section_mean,
                                                                  std::vector< std::vector<float> > & cross_section_stdabw)
{
    std::vector<float> c_features;

    //along the cross section the mean and stdabw gray value for each pixel has been computed and is now taken along the whole arc
    //possibly some of these features can improve the random forest classification, not tested yet
    c_features.push_back(get_mean(cross_section_mean[arc_index]));
    c_features.push_back(get_standard_deviation(cross_section_mean[arc_index]));
    c_features.push_back(get_quantile(cross_section_mean[arc_index],0.25));
    c_features.push_back(get_quantile(cross_section_mean[arc_index],0.5));
    c_features.push_back(get_quantile(cross_section_mean[arc_index],0.75));
    c_features.push_back(get_min(cross_section_mean[arc_index]));
    c_features.push_back(get_max(cross_section_mean[arc_index]));
 
    c_features.push_back(get_mean(cross_section_stdabw[arc_index]));
    c_features.push_back(get_standard_deviation(cross_section_stdabw[arc_index]));
    c_features.push_back(get_quantile(cross_section_stdabw[arc_index],0.25));
    c_features.push_back(get_quantile(cross_section_stdabw[arc_index],0.5));
    c_features.push_back(get_quantile(cross_section_stdabw[arc_index],0.75));
    c_features.push_back(get_min(cross_section_stdabw[arc_index]));
    c_features.push_back(get_max(cross_section_stdabw[arc_index]));
 
    return c_features;
 
}
 
std::vector<float> compute_single_arc_features_gray(std::vector< std::vector<point> > const & arcs,
                                                    std::vector< std::vector<float> > const & curvs,
                                                    std::vector <float> & area_mean,
                                                    std::vector <long> & area_size,
                                                    std::vector< std::vector<float> > & cross_section_mean,
                                                    std::vector< std::vector<float> > & cross_section_stdabw,
                                                    marray::Marray<unsigned int> const & two_boundings,
                                                    float * feature_array,
                                                    vigra::BImage const & main_image,
                                                    vigra::BImage const & prob_map,
                                                    int arc_index,
                                                    int dim_x,
                                                    int dim_y,
                                                    int max_features,
                                                    std::string param_file_name)
 
{
    // some important parameters should not need to be changed in the source code
    // so let's load them from a parameter file
    ParameterFile paramFile;
 
    if( !paramFile.load(param_file_name) )
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }
 
    //select boundary features for random forest by setting this boolean variables in the parameter file
	Parameter<bool> feature1;
	feature1.assign("", "feature1", true);
    feature1.load(paramFile,"config");
 
//    std::cout << "Parameter feature1: " << feature1 << std::endl;
 
	Parameter<bool> feature2;
	feature2.assign("", "feature2", true);
    feature2.load(paramFile,"config");
 
//    std::cout << "Parameter feature2: " << feature2 << std::endl;
 
	Parameter<bool> feature3;
	feature3.assign("", "feature3", true);
    feature3.load(paramFile,"config");
 
//    std::cout << "Parameter feature3: " << feature3 << std::endl;
 
	Parameter<bool> feature4;
	feature4.assign("", "feature4", true);
    feature4.load(paramFile,"config");
 
//    std::cout << "Parameter feature4: " << feature4 << std::endl;
 
	Parameter<bool> feature5;
	feature5.assign("", "feature5", true);
    feature5.load(paramFile,"config");
 
//    std::cout << "Parameter feature5: " << feature5 << std::endl;
 
	Parameter<bool> curvature;
	curvature.assign("", "curvature", true);
    curvature.load(paramFile,"config");
 
//    std::cout << "Parameter curvature: " << curvature << std::endl;
 
	Parameter<bool> probmap;
	probmap.assign("", "probmap", true);
    probmap.load(paramFile,"config");
 
//    std::cout << "Parameter probmap: " << probmap << std::endl;
 
	Parameter<bool> arcsize;
	arcsize.assign("", "arcsize", true);
    arcsize.load(paramFile,"config");
 
//    std::cout << "Parameter arcsize: " << arcsize << std::endl;
 
	Parameter<bool> region;
	region.assign("", "region", true);
    region.load(paramFile,"config");
 
//    std::cout << "Parameter region: " << region << std::endl;

	Parameter<bool> cross_section;
	cross_section.assign("", "cross_section", true);
    cross_section.load(paramFile,"config"); 
 
    std::vector<float> features;
    //We compute first the "pixel-feature based"-boundary features
    //we loop over all features
    for (int i=1;i<max_features+1;i++)
    {
        bool feature_calculation=true;
        if (i==1 && feature1==0) feature_calculation=false;
        if (i==2 && feature2==0) feature_calculation=false;
        if (i==3 && feature3==0) feature_calculation=false;
        if (i==4 && feature4==0) feature_calculation=false;
        if (i==5 && feature5==0) feature_calculation=false;
        if (feature_calculation==true)//only load pixel features when selected
        {
            //values stores the "pixel-based" features (loaded from the *.bin) for one arc (for each pixel in the arc the value)
            std::vector<float> pixel_feature_values;
            //std::cout<<"pixel_features_values.size():"<<pixel_feature_values.size()<<std::endl;
            pixel_feature_values=pixel_feature_values_on_one_arc(arcs,feature_array,arc_index,dim_x,dim_y,i,max_features);
            features.push_back(get_mean(pixel_feature_values));
            features.push_back(get_standard_deviation(pixel_feature_values));
            features.push_back(get_quantile(pixel_feature_values,0.25));
            features.push_back(get_quantile(pixel_feature_values,0.5));
            features.push_back(get_quantile(pixel_feature_values,0.75));
            features.push_back(get_min(pixel_feature_values));
            features.push_back(get_max(pixel_feature_values));
        }
    }
 
    //NOW WE WANT TO ADD THE CURVATURE FEATURES
    //curvature is a vector with the curvature for each pixel of the arc at arc_index
    if (curvature==true)
    {
        std::vector<float> curvature_values;
        curvature_values=curvs[arc_index];
        features.push_back(get_mean(curvature_values));
//        features.push_back(get_standard_deviation(curvature_values));
//        features.push_back(get_quantile(curvature_values,0.25));
//        features.push_back(get_quantile(curvature_values,0.5));
//        features.push_back(get_quantile(curvature_values,0.75));
//        features.push_back(get_min(curvature_values));
//        features.push_back(get_max(curvature_values));
    }
 
    //NOW WE WANT TO ADD THE PROBABILITY MAP FEATURES
    if (probmap==true)
    {
        std::vector<float> prob_map_values;
        prob_map_values=prob_map_values_on_one_arc(arcs,arc_index,prob_map);
        features.push_back(get_mean(prob_map_values));
//        features.push_back(get_standard_deviation(prob_map_values));
//        features.push_back(get_quantile(prob_map_values,0.25));
//        features.push_back(get_quantile(prob_map_values,0.5));
//        features.push_back(get_quantile(prob_map_values,0.75));
//        features.push_back(get_min(prob_map_values));
//        features.push_back(get_max(prob_map_values));
    }
 
    //NOW WE ADD THE SIZE OF ONE ARC AS A FEATURE
    if (arcsize==true) features.push_back( ((float)arcs[arc_index].size()) );
 
    //NOW WE COMPUTE THE REGION BASED FEATURES .. there is 1 region feature
        //-diffence of the mean of the neighbour regions
    if (region==true)
    {
       std::vector<float> region_features=compute_single_arc_region_features_gray(arc_index,area_mean,area_size,two_boundings);
 
       for(size_t rr=0;rr<region_features.size();rr++)
       {
           features.push_back(region_features[rr]);
       }
    }

    //NOW WE COMPUTE THE CROSS-SECTION BASED FEATURES
    if (cross_section==true)
    {
       std::vector<float> cross_section_features=compute_single_arc_cross_section_features_gray(arc_index,cross_section_mean,cross_section_stdabw);
 
       for(size_t rr=0;rr<cross_section_features.size();rr++)
       {
           features.push_back(cross_section_features[rr]);
       }
    }
 
    return features;
}
 
/*! \fn boundary_features_gray(std::string source_path_image,
                            std::string source_ws_image_path,
                            std::string source_path_prob_map,
                            std::string source_path_pixel_features,
                            std::string dest_path_boundary_features,
                            std::string param_file_name,
                            std::string filepath_thumbs="no")
 * \brief Calculate boundary features for an image.
 * \param source_path_image Filepath to the source image
 * \param source_ws_image_path Filepath to the source image's watershed segmentation image
 * \param source_path_prob_map Filepath to the source image's probabilty map image
 * \param source_path_pixel_features Filepath to the source image's pixel features
 * \param dest_path_boundary_features Filepath the boundary features will be written to
 * \param param_file_name Filename of the parameter file
 */ 
void boundary_features_gray(std::string source_path_image,
                            std::string source_ws_image_path,
                            std::string source_path_prob_map,
                            std::string source_path_pixel_features,
                            std::string dest_path_boundary_features,
                            std::string param_file_name,
                            std::string filepath_thumbs="no")
{
    //IMPORT RESULTS FROM HDF5 file
    //data structur
    //arcs == two cells in cartesian coordinats
    std::vector< std::vector<point> > arcs;
    std::vector<point> junctions;
    std::vector< std::vector<float> > phis;
    std::vector< std::vector<float> > curvs;
    std::vector< std::vector<float> > cross_section_mean;
    std::vector< std::vector<float> > cross_section_stdabw;
 
    marray::Marray<unsigned int> one_boundings;
    marray::Marray<unsigned int> two_boundings;
 
    vigra::BasicImage<unsigned int> ws_region_image;
    int dim_x, dim_y;
 
    std::string filename_ws=get_filename(source_path_image);
    filename_ws.append(".h5");
    std::string filepath_to_ws_region_image=source_ws_image_path;
    filepath_to_ws_region_image.append(filename_ws);
 
    load_cgp_data_structure(filepath_to_ws_region_image,
                            ws_region_image,
                            one_boundings,
                            two_boundings,
                            arcs,
                            junctions,
                            dim_x,
                            dim_y,
                            false);
 
    int nr_of_pixel_features;
 
    //WE OPEN THE IMAGE
    std::string filepath_to_image=source_path_image;
    vigra::ImageImportInfo info_img(filepath_to_image.c_str());
    dim_x=info_img.width();
    dim_y=info_img.height();
    vigra::BImage main_image(dim_x,dim_y);
    vigra::importImage(info_img, destImage(main_image));
 
    ParameterFile paramFile;
 
    if( !paramFile.load(param_file_name) )
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }

    //NOW DETERMINE IMAGE SELECTION
    bool thumb=false;
    vigra::BasicImage<bool> selection_image(dim_x,dim_y);
    for (int x=0; x<dim_x; x++)
    {
        for (int y=0; y<dim_y; y++)
        {
            selection_image(x,y)=true;
        }
    }

    if (filepath_thumbs!="no")
    {
        filepath_thumbs.append(get_filename(source_path_image));
        filepath_thumbs.append(".reduced.bmp");

        std::ifstream thumb_file(filepath_thumbs.c_str());

        if(thumb_file)
        {
            thumb_file.close();
            std::cout<<"Reduced thumb found!"<<std::endl;

            //open thumb
            vigra::ImageImportInfo info(filepath_thumbs.c_str());
            int thumb_dim_x=0.2*dim_x;
            int thumb_dim_y=0.2*dim_y;

            if (info.width()!=thumb_dim_x || info.height()!=thumb_dim_y)
                std::cout << "Error: thumb has not expected size of "<< thumb_dim_x <<" x "<<thumb_dim_y<<"!"<<std::endl;
            else
            {
                vigra::IImage thumb_grayvalues(thumb_dim_x,thumb_dim_y);
                vigra::IImage label_image(thumb_dim_x,thumb_dim_y);
                vigra::BasicImage<bool> thumb_bool(thumb_dim_x,thumb_dim_y);

                importImage(info, destImage(thumb_grayvalues));

                for (int x=0;x<thumb_dim_x;x++)
                    for (int y=0;y<thumb_dim_y;y++)
                    {
                        if (thumb_grayvalues(x,y)<255) thumb_grayvalues(x,y)=0;
                        else thumb_grayvalues(x,y)=1;
                    }

                int nr_regions = vigra::labelImageWithBackground(vigra::srcImageRange(thumb_grayvalues), vigra::destImage(label_image), 0, 0);
                std::vector<bool> regions(nr_regions,false);

                //check for regions connected to border
                for (int y=0; y<thumb_dim_y; y++)
                {
                    if (label_image(0,y)>0) regions[label_image(0,y)-1]=true;
                    if (label_image(thumb_dim_x-1,y)>0) regions[label_image(thumb_dim_x-1,y)-1]=true;
                }

                for (int x=0; x<thumb_dim_x; x++)
                {
                    if (label_image(x,0)>0) regions[label_image(x,0)-1]=true;
                    if (label_image(x,thumb_dim_y-1)>0) regions[label_image(x,thumb_dim_y-1)-1]=true;
                }

                for (int x=0;x<thumb_dim_x;x++)
                    for (int y=0;y<thumb_dim_y;y++)
                        if (label_image(x,y)>0)
                        {
                            if (regions[label_image(x,y)-1]) thumb_bool(x,y)=false;
                            else thumb_bool(x,y)=true;
                        }
                        else thumb_bool(x,y)=true;

                for (int region=0; region<nr_regions; region++)
                    if (!regions[region])
                    {
                        std::cout<<"Inside white region ignored"<<std::endl;
                    }

                resizeImageNoInterpolation(srcImageRange(thumb_bool),destImageRange(selection_image));

                for (int x=0; x<dim_x; x++)
                {
                    selection_image(x,0)=false;
                    selection_image(x,dim_y-1)=false;
                }

                for (int y=0; y<dim_y; y++)
                {
                    selection_image(0,y)=false;
                    selection_image(dim_x-1,y)=false;
                }

                thumb=true;
            }
        }
        else std::cout<<"Thumb "<<filepath_thumbs<<" not found"<<std::endl;
    }

    if (!thumb)
    {
        //to get the complete filepath "dest_path_boundary_features" is only the path of the folder
        std::string filepath_to_image_selection=dest_path_boundary_features;
        filepath_to_image_selection.append(get_filename(source_path_image));
        filepath_to_image_selection.append(".selection.dat");

        //check if selection file exists
        std::ifstream selection_file_in(filepath_to_image_selection.c_str());

        int low_x, low_y, high_x, high_y;

        if (selection_file_in.is_open())
        {
            std::cout<<"image selection loaded from file"<<std::endl;
            selection_file_in>>low_x;
            selection_file_in>>high_x;
            selection_file_in>>low_y;
            selection_file_in>>high_y;

            selection_file_in.close();
        }
        else
        {
            std::cout<<"determine image selection"<<std::endl;

            //default values
            low_x=0;
            high_x=dim_x-1;
            low_y=0;
            high_y=dim_y-1;

            //we write the resulting values to file    
            std::ofstream selection_file(filepath_to_image_selection.c_str());

            //we compare the mean gray value of blocks of this size
	        Parameter<int> blocksize;
	        blocksize.assign("", "blocksize", 100);
            blocksize.load(paramFile,"config");

            //if the gray value average differs more than this value a border is detected
	        Parameter<int> grayvalue_difference;
	        grayvalue_difference.assign("", "grayvalue_difference", 10);
            grayvalue_difference.load(paramFile,"config");

            int old_block_grayvalue=-1;
            int block_grayvalue=0;
            bool written=false;

            //start in the center and go to the left
            for (int x=dim_x/2; x>blocksize; x-=blocksize)
            {
                for (int y=0; y<dim_y; y++)
                {
                    for (int line=x; line>x-blocksize; line--)
                        block_grayvalue+=main_image(line,y);
                }
                block_grayvalue=block_grayvalue/(blocksize*dim_y);
                //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
                if (fabs(block_grayvalue-old_block_grayvalue)>grayvalue_difference && old_block_grayvalue>-1)
                {
                    std::cout<<"low x found!"<<std::endl;
                    selection_file <<x<<" ";
                    low_x=x;
                    written=true;
                    break;
                }
                old_block_grayvalue=block_grayvalue;
                block_grayvalue=0;
            }

            //if no border has been found write 0 to file
            if (!written) selection_file <<0<<" ";
            written=false;  
            old_block_grayvalue=-1;

            //start in the center and go to the right
            for (int x=dim_x/2; x<dim_x-blocksize; x+=blocksize)
            {
                for (int y=0; y<dim_y; y++)
                {
                    for (int line=x; line<x+blocksize; line++)
                        block_grayvalue+=main_image(line,y);
                }
                block_grayvalue=block_grayvalue/(blocksize*dim_y);
                //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
                if (fabs(block_grayvalue-old_block_grayvalue)>grayvalue_difference && old_block_grayvalue>-1)
                {
                    std::cout<<"high x found!"<<std::endl;
                    selection_file <<x<<" ";
                    high_x=x;
                    written=true;
                    break;
                }
                old_block_grayvalue=block_grayvalue;
                block_grayvalue=0;
            }

            if (!written) selection_file <<dim_x-1<<" ";
            written=false;  
            old_block_grayvalue=-1;

            //start in the center and go up
            for (int y=dim_y/2; y>blocksize; y-=blocksize)
            {
                for (int x=0; x<dim_x; x++)
                {
                    for (int line=y; line>y-blocksize; line--)
                        block_grayvalue+=main_image(x,line);
                }
                block_grayvalue=block_grayvalue/(blocksize*dim_x);
                //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
                if (fabs(block_grayvalue-old_block_grayvalue)>grayvalue_difference && old_block_grayvalue>-1)
                {
                    std::cout<<"low y found!"<<std::endl;
                    selection_file <<y<<" ";
                    low_y=y;
                    written=true;
                    break;
                }
                old_block_grayvalue=block_grayvalue;
                block_grayvalue=0;
            }

            if (!written) selection_file <<0<<" ";
            written=false;  
            old_block_grayvalue=-1;

            //start in the center and go down
            for (int y=dim_y/2; y<dim_y-blocksize; y+=blocksize)
            {
                for (int x=0; x<dim_x; x++)
                {
                    for (int line=y; line<y+blocksize; line++)
                        block_grayvalue+=main_image(x,line);
                }
                block_grayvalue=block_grayvalue/(blocksize*dim_x);
                //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
                if (fabs(block_grayvalue-old_block_grayvalue)>grayvalue_difference && old_block_grayvalue>-1)
                {
                    std::cout<<"high y found!"<<std::endl;
                    selection_file <<y<< "\n";
                    high_y=y;
                    written=true;
                    break;
                }
                old_block_grayvalue=block_grayvalue;
                block_grayvalue=0;
            }

            if (!written) selection_file <<dim_y-1<< "\n";

            selection_file.close();
        }

        for (int x=0; x<=low_x; x++)
        {
            for (int y=0; y<dim_y; y++)
            {
                selection_image(x,y)=false;
            }
        }

        for (int x=low_x+1; x<high_x; x++)
        {
            for (int y=0; y<=low_y; y++)
            {
                selection_image(x,y)=false;
            }

            for (int y=high_y; y<dim_y; y++)
            {
                selection_image(x,y)=false;
            }
        }

        for (int x=high_x; x<dim_x; x++)
        {
            for (int y=0; y<dim_y; y++)
            {
                selection_image(x,y)=false;
            }
        }
    }

    //fill the vector areas containing all pixels of each area
    std::vector< std::vector<point> > areas;
    //counts the number of pixels outside image selection
    std::vector<long> outside_pixels;
    int nr_areas=0;
 
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           point p;
           p.x=x;
           p.y=y;

           if (ws_region_image(x,y)>nr_areas)
           {
               nr_areas=ws_region_image(x,y);
               areas.resize(nr_areas);
               outside_pixels.resize(nr_areas,0);
           } 

           areas[ws_region_image(x,y)-1].push_back(p);
           if (!selection_image(x,y)) outside_pixels[ws_region_image(x,y)-1]++;
        }
    }

    //in order to save memory all area mean gray values are computed now instead of doing this for the neighbors of each arc
    std::vector <float> area_mean;
    area_mean.resize(nr_areas);
    std::vector <long> area_size;
    area_size.resize(nr_areas);
 
    std::cout << "calculating mean gray-values and area size..." << std::endl;
 
    //visit all areas and compute the mean gray-value
    for(int area=0; area<nr_areas; area++)
    {
        //vector to store all the pixel values of the area
        std::vector<float> region_pixels;

        //if outside ratio is very low or very high, take all area pixels
        if(outside_pixels[area]<0.1*areas[area].size() || outside_pixels[area]>0.9*areas[area].size() || areas[area].size()-outside_pixels[area]<500)
        { 
            //fill the region_pixels with values
            for(size_t i=0;i<areas[area].size();i++)
            {
                region_pixels.push_back(   (float)(main_image(areas[area][i].x, areas[area][i].y))     );
            }
     
            //mean of this area
            area_mean[area]=get_mean(region_pixels);

            //size of this area
            area_size[area]=(long)areas[area].size();
        }
        else//take only inside pixels (which are at least 500)
        {
            //fill the region_pixels with values
            for(size_t i=0;i<areas[area].size();i++)
            {
                int x=areas[area][i].x;
                int y=areas[area][i].y;

                if (selection_image(x,y)) region_pixels.push_back((float)(main_image(x,y)));
            }
     
            //mean of inside pixels of this area
            area_mean[area]=get_mean(region_pixels);

            //size of inside pixels of this area
            area_size[area]=(long)(areas[area].size()-outside_pixels[area]);
        }
    }
 
    areas.clear();
    outside_pixels.clear();
 
    std::cout << "...done" << std::endl;

	Parameter<bool> probmap;
	probmap.assign("", "probmap", true);
    probmap.load(paramFile,"config");
 
    //NOW WE HAVE TO OPEN THE PROBABILITY MAP
    vigra::BImage prob_map(dim_x,dim_y);
 
    if (probmap==true)//only if needed
    {
        std::string filepath_to_prob_map=source_path_prob_map;
        std::string image_filename=get_filename(source_path_image);
        filepath_to_prob_map.append(image_filename);
        filepath_to_image.resize(filepath_to_image.size()-4);
        vigra::ImageImportInfo info_prob(filepath_to_prob_map.c_str());
        vigra::importImage(info_prob, destImage(prob_map));
    }
 
	Parameter<bool> curvature;
	curvature.assign("", "curvature", true);
    curvature.load(paramFile,"config");

	Parameter<bool> cross_section;
	cross_section.assign("", "cross_section", true);
    cross_section.load(paramFile,"config");
 
    if (curvature || cross_section)//only if needed
    {
        // calculate angle phi of local normal vectors
        std::cout << "calculating angles" << std::endl;
 
        for(int i=0; i<arcs.size(); i++)
        {
            std::vector<float> phi;
            phis.push_back(phi);
            if(arcs[i].size() < 5 )
            {
                // if the size of an arc is smaller than 5, use a simpler approximation of the derivative
 
                //std::cout << "Short arc... " << arcs[i].size() << " ";
                double phi0 = atan2(arcs[i][1].y - arcs[i][0].y, arcs[i][1].x - arcs[i][0].x );
                phis.back().push_back( phi0 + PI/2.0f );
 
                for(int j=1; j<arcs[i].size()-1; j++)
                {
                    double phi_temp = atan2(arcs[i][j+1].y - arcs[i][j-1].y, arcs[i][j+1].x - arcs[i][j-1].x );
                    phis.back().push_back( phi_temp + PI/2.0f );
                }
 
                double phi_end = atan2(arcs[i].back().y - arcs[i][arcs[i].size()-2].y, arcs[i].back().x - arcs[i][arcs[i].size()-2].x );
                phis.back().push_back( phi_end + PI/2.0f );
                //std::cout << "done!" << std::endl;
            }
            else
            {
                //std::cout << "Long arc... " << arcs[i].size() << " ";
                calculatePhi( &(arcs[i]), &(phis[i]));
                //std::cout << "done!"<<std::endl;
            }
        }
 
        //calculate the curvature
        std::cout << "calculating curvature" << std::endl;
        for(int i=0; i<arcs.size(); i++)
        {
            std::vector<float> curv;
            curvs.push_back(curv);
            calculateCurv( &(phis[i]), &(curvs[i]) );
        }

        if (cross_section)
        {
            //calculate cross-section
            std::cout << "calculating cross-section" << std::endl;
     
            cross_section_mean.resize(arcs.size());
            cross_section_stdabw.resize(arcs.size());
     
            for(int i=0; i<arcs.size(); i++)
            {
                std::vector<std::vector<float> > cross_arc;
                calculateCrossSection( &(arcs[i]), &(phis[i]), &main_image, &(cross_arc) );

                //for each arc we have a double vector: for each arc pixel a cross section is taken
                for(int j=0; j<cross_arc.size(); j++)
                {
                    //calculate mean and standard deviation for the cross section at every arc pixel
                    cross_section_mean[i].push_back(get_mean(cross_arc[j]));
                    cross_section_stdabw[i].push_back(get_standard_deviation(cross_arc[j]));
                }
            }
        } 

        std::cout<<"...done"<<std::endl;
    }
 
	Parameter<bool> feature1;
	feature1.assign("", "feature1", true);
    feature1.load(paramFile,"config");
 
	Parameter<bool> feature2;
	feature2.assign("", "feature2", true);
    feature2.load(paramFile,"config");
 
	Parameter<bool> feature3;
	feature3.assign("", "feature3", true);
    feature3.load(paramFile,"config");
 
	Parameter<bool> feature4;
	feature4.assign("", "feature4", true);
    feature4.load(paramFile,"config");
 
	Parameter<bool> feature5;
	feature5.assign("", "feature5", true);
    feature5.load(paramFile,"config");
 
    if (feature1==true || feature2==true || feature3==true || feature4==true || feature5==true)//only if needed
    {
        //now we load the feature file
        //OPEN THE FEATURE FILE
        std::cout<<"open pixel-feature file..."<<std::endl;
        std::string only_filename_of_the_image=get_filename(source_path_image);
        std::string filepath_to_pixel_features=source_path_pixel_features;
        filepath_to_pixel_features.append(only_filename_of_the_image);
        filepath_to_pixel_features.append(".bin");
        std::cout<<filepath_to_pixel_features<<std::endl;
 
        float temp_image;//variable to determine the nr of features
 
        FILE *fp_unknown_features;
        fp_unknown_features =fopen(filepath_to_pixel_features.c_str(),"rb");
 
        //check if file exists
        if(fp_unknown_features==NULL)
        {
            std::cout<<"Error: Feature File "<<filepath_to_pixel_features<<" is NOT existend"<<std::endl;
            exit(-1);
        }
        else
        {
            //determine the size of the featurefile, thus the nr of features
            nr_of_pixel_features = -1;
            while(!feof(fp_unknown_features))
            {
                for(int y=0;y<dim_y;y++)
                {
                    for(int x=0;x<dim_x;x++)
                    {
                        fread(&temp_image,sizeof(float),1,fp_unknown_features);
                    }
                }
                nr_of_pixel_features++;
            }
        }
 
        fclose(fp_unknown_features);
 
        //check if nr of features was read in correctly
        if(nr_of_pixel_features<1)
        {
            std::cout<<"Error: Nr of Features was not read in correctly from "<<filepath_to_pixel_features<<std::endl;
            exit(-1);
        }
 
        std::cout<<"nr of pixel features: "<<nr_of_pixel_features<<std::endl;
 
        //reopen featurefile, now we store the features in an array
        fp_unknown_features =fopen(filepath_to_pixel_features.c_str(),"rb");
 
        //an array to store all the features!
        float * pixel_features_array=new float[dim_x*dim_y*nr_of_pixel_features];
 
        size_t read_in=fread(pixel_features_array,sizeof(float),dim_x*dim_y*nr_of_pixel_features,fp_unknown_features);
        if(read_in!=dim_x*dim_y*nr_of_pixel_features)
        {
            std::cout<<"boundary_features(..) pixel_feature_file size is wrong"<<std::endl;
            exit(-1);
        }
 
        std::cout<<"...done"<<std::endl;
 
        //now we have to open/create the binary file in which the features of one image are stored
        //to do that we have to use the path "dest n_path_boundary_features"
        //later we MUST KNOW how many features are in the binary file!
        //  => the first float we store in the binary file HAS TO BE the size of data (=boundary_features*arc.size())
 
        //to get the complete filepath "dest_path_boundary_features" is only the path of the folder
        std::string filepath_to_boundary_features=dest_path_boundary_features;
        filepath_to_boundary_features.append(get_filename(source_path_image));
        filepath_to_boundary_features.append(".bin");
 
        //open a filepointer
        std::cout<<"open boundary-feature file..."<<std::endl;
        std::cout<<filepath_to_boundary_features<<std::endl;
        FILE * fp;
        fp =fopen(filepath_to_boundary_features.c_str(),"w+b");
        //IS THIS SAVE???? (to make -1 )    //has been dropped
 
        //first entry in the file should be the nr of arcs
        float float_arc_size=(float)arcs.size();
        fwrite(&float_arc_size,sizeof(float),1,fp);
        std::cout<<"...done"<<std::endl;

        //now we compute one feature and see what size the arc_feature_vector is
        //and write this as second entry!
        std::vector<float> arc_feature;
        arc_feature=compute_single_arc_features_gray(arcs,curvs,area_mean,area_size,cross_section_mean,cross_section_stdabw,two_boundings,
                                                     pixel_features_array,main_image,prob_map,0,dim_x,dim_y,nr_of_pixel_features,param_file_name);
        float float_size=(float)arc_feature.size();
 
        fwrite(&float_size,sizeof(float),1,fp);
 
        std::cout<<"compute boundary-features..."<<std::endl;
        std::cout<<"nr of boundary features: "<<arc_feature.size()<<std::endl;

        //now we compute the feature for each arc and write each value of the feature vector of each arc in a binary file
        for(int i=0;i<(int)arcs.size();i++)
        {
            arc_feature=compute_single_arc_features_gray(arcs,curvs,area_mean,area_size,cross_section_mean,cross_section_stdabw,two_boundings,
                                                         pixel_features_array,main_image,prob_map,i,dim_x,dim_y,nr_of_pixel_features,param_file_name);
 
            //std::cout<<"size arc_feature:"<<arc_feature.size()<<std::endl;
            for(int jj=0;jj<(int)arc_feature.size();jj++)
            {
              float value_to_write=arc_feature[jj];
              //write the value into the filepointer
              fwrite(&value_to_write,sizeof(float),1,fp);
            }
        }
        std::cout<<"...done"<<std::endl;

        fclose(fp);
 
        delete pixel_features_array;
 
        arcs.clear();
    }
    else
    {
        //an pseudo array of pixel features!
        float * pixel_features_array=new float[1];
        nr_of_pixel_features = 0;
 
        //now we have to open/create the binary file in which the features of one image are stored
        //to do that we have to use the path "dest n_path_boundary_features"
        //later we MUST KNOW how many features are in the binary file!
        //  => the first float we store in the binary file HAS TO BE the size of data (=boundary_features*arc.size())
 
        //to get the complete filepath "dest_path_boundary_features" is only the path of the folder
        std::string filepath_to_boundary_features=dest_path_boundary_features;
        filepath_to_boundary_features.append(get_filename(source_path_image));
        filepath_to_boundary_features.append(".bin");
 
        //open a filepointer
        std::cout<<"open boundary-feature file..."<<std::endl;
        std::cout<<filepath_to_boundary_features<<std::endl;
        FILE * fp;
        fp =fopen(filepath_to_boundary_features.c_str(),"w+b");
        //IS THIS SAVE???? (to make -1 )    //has been dropped
 
        //first entry in the file should be the nr of arcs
        float float_arc_size=(float)arcs.size();
        fwrite(&float_arc_size,sizeof(float),1,fp);
        std::cout<<"...done"<<std::endl;
 
        //now we compute one feature and see what size the arc_feature_vector is
        //and write this as second entry!
        std::vector<float> arc_feature;
 
        arc_feature=compute_single_arc_features_gray(arcs,curvs,area_mean,area_size,cross_section_mean,cross_section_stdabw,two_boundings,pixel_features_array,
                                                     main_image,prob_map,0,dim_x,dim_y,nr_of_pixel_features,param_file_name);
        float float_size=(float)arc_feature.size();
 
        fwrite(&float_size,sizeof(float),1,fp);
 
        std::cout<<"compute boundary-features..."<<std::endl;
        std::cout<<"nr of boundary features: "<<arc_feature.size()<<std::endl;

        //now we compute the feature for each arc and write each value of the feature vector of each arc in a binary file
        for(int i=0;i<(int)arcs.size();i++)
        {
            arc_feature=compute_single_arc_features_gray(arcs,curvs,area_mean,area_size,cross_section_mean,cross_section_stdabw,two_boundings,pixel_features_array,
                                                         main_image,prob_map,i,dim_x,dim_y,nr_of_pixel_features,param_file_name);
 
            //std::cout<<"size arc_feature:"<<arc_feature.size()<<std::endl;
            for(int jj=0;jj<(int)arc_feature.size();jj++)
            {
              float value_to_write=arc_feature[jj];
              //write the value into the filepointer
              fwrite(&value_to_write,sizeof(float),1,fp);
            }
        }
        std::cout<<"...done"<<std::endl;
 
        fclose(fp);
 
        /*
        //TESTING
        vigra::FImage curv(dim_x,dim_y);
 
        //loop over all arcs, boundaries of bubble/grain boundaries to output image
        for(int a=0;a<(int)arcs.size();a++)
        {
            //TEST CURV VALUES
            std::vector<point>  this_arc;
            this_arc=arcs[a];
            //now we loop over the points in this arc
            for(int p=0;p<(int)this_arc.size();p++)
            {
                int x=this_arc[p].x;
                int y=this_arc[p].y;
                if (get_mean(curvs[a])>0.12f) curv(x,y)=1;        
                else curv(x,y)=0;
            }
        }
        exportImage(srcImageRange(curv), vigra::ImageExportInfo("curvs.bmp"));
        */
 
        delete pixel_features_array;
 
        arcs.clear();
    }
}
