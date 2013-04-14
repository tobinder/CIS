/*! \file boundary_features.h
 * \brief Boundary feature calculation for colour images.
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

//we need a funktion that computes the index in the complete feature file
//we want wo give a pixel(x,y) and we want say wich feature
// the return should just be a index to the corresponding array
int get_feature_array_index(int x,int y,int dim_x,int dim_y,int feature_nr,int max_features)
{
 //the features are stored in a binary file: max_features floats pro pixel
 //the order of the file is:
 //  LUV....LUV     (21 times 3 pixels (=63) )


 //first we have to know in wich "max_feature"-block we want to use (=>means wich for
 //wich pixel we want to compute the feature:
 // block_index=dim_x*y+x

 int block_index=dim_x*y+x;
 int array_index=block_index*max_features+feature_nr;
 return array_index;

}

//IS USED TO SAVE THE VALUES OF ONe "PIXEL BASED FEATURE"  OF ONE ARC IN A STD::VECTOR
std::vector<float> pixel_feature_values_on_one_arc(std::vector< std::vector<point> > const & arcs,
                                                    float * feature_array,
                                                    int arc_index,
                                                    int dim_x,
                                                    int dim_y,
                                                    int feature_nr,
                                                    int max_features)
{
   std::vector<float> value_vector;

   for(int p=0; p<(int)arcs[arc_index].size();p++ )
   {
        int x_pixel=arcs[arc_index][p].x;
        int y_pixel=arcs[arc_index][p].y;

        int array_index=get_feature_array_index(x_pixel,y_pixel,dim_x,dim_y,feature_nr,max_features);
        float value=feature_array[array_index];
        value_vector.push_back(value);
   }

   return value_vector;
}

std::vector<float> prob_map_values_on_one_arc(std::vector< std::vector<point> > const & arcs,
                                                    int arc_index,
                                                    vigra::BImage const & prob_map)
{
    std::vector<float> value_vector;

    for(int p=0; p<(int)arcs[arc_index].size();p++ )
    {
        int x_pixel=arcs[arc_index][p].x;
        int y_pixel=arcs[arc_index][p].y;

        float value=(float)prob_map(x_pixel,y_pixel);
        value_vector.push_back(value);
    }

    return value_vector;
}

//quantile
float get_quantile(std::vector<float> values,float q)
{

    int size=values.size();
    if(size!=1)
    {
        //sort the "values"-vector
        sort( values.begin(), values.end() );

        //f_index is the index of the lower quartil
        float f_index=((float)size-1)*q;
        //if we make f_index-floor(f_index) we know how to round
        //the quantile
        float round_value=f_index-floor(f_index);

        if(round_value==0)
        {
            return values[(int)f_index];
        }
        else
        {
            float return_value=0;
            int   low_index=floor(f_index);
            return_value=((1-round_value)*values[low_index]+ (round_value)*values[low_index+1]);
            return return_value;
        }
    }
    else
    {
        return values[1];
    }
}

//mean
float get_mean(std::vector<float> values)
{
    float mean=0;
    for(int i=0;i<(int)values.size();i++)
    {
        mean=mean+values[i];
    }

    return (mean/values.size());
}

//standard devitation
float get_standard_deviation(std::vector<float> values)
{
    float standard_deviation=0;

    //if arc is only 1 pixel long we can return 0 as standard deviation
    if(values.size()==1)
    {
        return 0;
    }
    else
    {
        float mean=get_mean(values);
        for(int i=0;i<(int)values.size();i++)
        {
            standard_deviation=standard_deviation + ((values[i]-mean)*(values[i]-mean));
        }
        standard_deviation=sqrt(standard_deviation*(1.0/((float)values.size()-1.0)));
        return standard_deviation;
    }

}

float get_min(std::vector<float> values)
{
    std::vector<float>::iterator pos;
    // find  minimum elements
    pos = min_element (values.begin(), values.end());
    return  *pos;
}

float get_max(std::vector<float> values)
{
    std::vector<float>::iterator pos;
    // find  maximum elements
    pos = max_element (values.begin(), values.end());
    return  *pos;
}

std::vector<float> compute_single_arc_region_features(int arc_index,
                                                      vigra::BRGBImage const & main_image,
                                                      std::vector < std::vector<point> > & areas,
                                                      marray::Marray<unsigned int> const & two_boundings,
                                                      std::vector<std::vector<point> > const & arcs)
{
    //first we need 6 std::vector<float> to store all the pixel values of two neigbour regions on each of the 3 colors channels
    std::vector<float> a0 , a1 , b0 , b1, c0 , c1 ;

    //now we have to look in the two boundings which regions are the neigbours
    size_t n0=two_boundings(arc_index,0);
    size_t n1=two_boundings(arc_index,1);


    //now we must fill a std::vector with the caartesian coordinates of the neigbours
    std::vector<point> c0_coord ,c1_coord;

    std::vector<point> area0 = areas[n0-1];
    std::vector<point> area1 = areas[n1-1];

    for(size_t k=0; k<area0.size(); ++k)
    {
        point p0;
        p0.x=area0[k].x;
        p0.y=area0[k].y;
        c0_coord.push_back(p0);
    }
    for(size_t k=0; k<area1.size(); ++k)
    {
        point p1;
        p1.x=area1[k].x;
        p1.y=area1[k].y;
        c1_coord.push_back(p1);
    }

    //now we can fill a0 ,....,c1 with values

    for(size_t i=0;i<c0_coord.size();i++)
    {
        a0.push_back(   (float)(main_image(c0_coord[i].x, c0_coord[i].y)[0])     );
        b0.push_back(   (float)(main_image(c0_coord[i].x, c0_coord[i].y)[1])     );
        c0.push_back(   (float)(main_image(c0_coord[i].x, c0_coord[i].y)[2])     );
    }

    for(size_t i=0;i<c1_coord.size();i++)
    {
        a1.push_back(   (float)(main_image(c1_coord[i].x, c1_coord[i].y)[0])     );
        b1.push_back(   (float)(main_image(c1_coord[i].x, c1_coord[i].y)[1])     );
        c1.push_back(   (float)(main_image(c1_coord[i].x, c1_coord[i].y)[2])     );
    }

    float mean_a0,mean_a1,mean_b0,mean_b1,mean_c0,mean_c1;
    float std_a0,std_a1,std_b0,std_b1,std_c0,std_c1;
    float size_0,size_1;

    mean_a0=get_mean(a0);
    mean_a1=get_mean(a1);
    mean_b0=get_mean(b0);
    mean_b1=get_mean(b1);
    mean_c0=get_mean(c0);
    mean_c1=get_mean(c1);

    std_a0=get_standard_deviation(a0);
    std_a1=get_standard_deviation(a1);
    std_b0=get_standard_deviation(b0);
    std_b1=get_standard_deviation(b1);
    std_c0=get_standard_deviation(c0);
    std_c1=get_standard_deviation(c1);

    size_0=(float)c0_coord.size();
    size_1=(float)c1_coord.size();

    std::vector<float> r_features;

    r_features.push_back(fabs(mean_a0-mean_a1)); //sehr gut
    r_features.push_back(fabs(mean_b0-mean_b1));
    r_features.push_back(fabs(mean_c0-mean_c1));

/*
    r_features.push_back(fabs(std_a0-std_a1));
    r_features.push_back(fabs(std_b0-std_b1));
    r_features.push_back(fabs(std_c0-std_c1));

    r_features.push_back(size_0+size_1);    //nicht gut!
    r_features.push_back(fabs(size_0-size_1));
*/

/*
    if(arc_index==100)
    {

        //DEBUG

        vigra::BRGBImage test_it(main_image.width(),main_image.height());

        //NOW WE DISPLAY ONE ARC ,lets say arg 100



        for(int iii=0;iii<c0_coord.size();iii++)
        {
            test_it(c0_coord[iii].x,c0_coord[iii].y)[0]=0;
            test_it(c0_coord[iii].x,c0_coord[iii].y)[1]=0;
            test_it(c0_coord[iii].x,c0_coord[iii].y)[2]=255;

        }

        for(int iii=0;iii<c1_coord.size();iii++)
        {
            test_it(c1_coord[iii].x,c1_coord[iii].y)[0]=0;
            test_it(c1_coord[iii].x,c1_coord[iii].y)[1]=255;
            test_it(c1_coord[iii].x,c1_coord[iii].y)[2]=0;

        }

        for(int iii=0;iii<arcs[100].size();iii++)
        {
            test_it(arcs[100][iii].x,arcs[100][iii].y)[0]=255;
            test_it(arcs[100][iii].x,arcs[100][iii].y)[1]=0;
            test_it(arcs[100][iii].x,arcs[100][iii].y)[2]=0;
        }

        exportImage(srcImageRange(test_it), vigra::ImageExportInfo("der_test.bmp"));

    }

    if(arc_index==200)
    {

        //DEBUG

        vigra::BRGBImage test_it(main_image.width(),main_image.height());

        //NOW WE DISPLAY ONE ARC ,lets say arg 100



        for(int iii=0;iii<c0_coord.size();iii++)
        {
            test_it(c0_coord[iii].x,c0_coord[iii].y)[0]=0;
            test_it(c0_coord[iii].x,c0_coord[iii].y)[1]=0;
            test_it(c0_coord[iii].x,c0_coord[iii].y)[2]=255;

        }

        for(int iii=0;iii<c1_coord.size();iii++)
        {
            test_it(c1_coord[iii].x,c1_coord[iii].y)[0]=0;
            test_it(c1_coord[iii].x,c1_coord[iii].y)[1]=255;
            test_it(c1_coord[iii].x,c1_coord[iii].y)[2]=0;

        }


        for(int iii=0;iii<arcs[200].size();iii++)
        {
            test_it(arcs[200][iii].x,arcs[200][iii].y)[0]=255;
            test_it(arcs[200][iii].x,arcs[200][iii].y)[1]=0;
            test_it(arcs[200][iii].x,arcs[200][iii].y)[2]=0;
        }

        exportImage(srcImageRange(test_it), vigra::ImageExportInfo("der_test2.bmp"));

    }
    */

    return r_features;
}

std::vector<float> compute_single_arc_features(std::vector< std::vector<point> > const & arcs,
                                               std::vector< std::vector<float> > const & curvs,
                                               std::vector < std::vector<point> > & areas,
                                               marray::Marray<unsigned int> const & two_boundings,
                                               float * feature_array,
                                               vigra::BRGBImage const & main_image,
                                               vigra::BImage const & prob_map,
                                               int arc_index,
                                               int dim_x,
                                               int dim_y,
                                               int max_features)
{
    std::vector<float> features;
    //We compute first the "pixe-feature based"-boundary features
    //we got "max_featues" and we take alway feature 0,1,2   9,10,11  18,19,20
    //this is to ensure that we take always one scale (the smalles) and take all 3 color channels
    //so first we loop over all the  (max_features/3) possible "pixel-features"
    //to get the right numbers ( 0,1,2   9,10,11  18,19,20....54,55,56 ,...,63,64,65) we loop twice

    int loop_size=max_features/(3*3);  //3 channels ,3 scales...,

    for (int i=0;i<loop_size;i++)
    {
        for(int j=0;j<3;j++)
        {
            //is used to get only the first scale of the pixel based features
            //feature 0,1,2 9,10,11 18,19,20.......52,53 54 (in the case of 63 features)

            int feature_nr= 3*3*i+j;

            //values stores the "pixel-based" features (loaded from the *.bin) for one arc (for each pixel in the arc the value of
            std::vector<float> pixel_feature_values;
           // std::cout<<"pixel_features_values.size():"<<pixel_feature_values.size()<<std::endl;
            pixel_feature_values=pixel_feature_values_on_one_arc(arcs,feature_array,arc_index,dim_x,dim_y,feature_nr,max_features);
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
    //curvature is a vector with the curvature foreach pixel of the arc at arc_index
    std::vector<float> curvature_values;
    curvature_values=curvs[arc_index];
    features.push_back(get_mean(curvature_values));
    features.push_back(get_standard_deviation(curvature_values));
    features.push_back(get_quantile(curvature_values,0.25));
    features.push_back(get_quantile(curvature_values,0.5));
    features.push_back(get_quantile(curvature_values,0.75));
    features.push_back(get_min(curvature_values));
    features.push_back(get_max(curvature_values));

    //NOW WE WANT TO ADD THE PROBABILITY MAP FEATURES
    std::vector<float> prob_map_values;
    prob_map_values=prob_map_values_on_one_arc(arcs,arc_index,prob_map);
    features.push_back(get_mean(prob_map_values));  //SEHR GUT  127
    features.push_back(get_standard_deviation(prob_map_values));
    features.push_back(get_quantile(prob_map_values,0.25)); //SEHR GUT  129
    features.push_back(get_quantile(prob_map_values,0.5));  //SEHR GUT2 130
    features.push_back(get_quantile(prob_map_values,0.75));  //SEHR GUT  131
    features.push_back(get_min(prob_map_values));
    features.push_back(get_max(prob_map_values));   //SEHR GUT

    //NOW WE ADD THE SIZE OF ONE ARC AS A FEATURE
    features.push_back(  ((float)arcs[arc_index].size())    );

    //NOW WE COMPUTE THE REGION BASED FEATURES .. there is 1 region feature
        //-diffence of the mean of the neighbour regions
    std::vector<float> region_features=compute_single_arc_region_features(arc_index,main_image,areas,two_boundings,arcs);

    for(size_t rr=0;rr<region_features.size();rr++)
    {
        features.push_back(region_features[rr]);
    }

    return features;
}

/*! \fn boundary_features(std::string source_path_image,
                       std::string source_ws_image_path,
                       std::string source_path_prob_map,
                       std::string source_path_pixel_features,
                       std::string dest_path_boundary_features)
 * \brief Calculate boundary features for an image.
 * \param source_path_image Filepath to the source image
 * \param source_ws_image_path Filepath to the source image's watershed segmentation image
 * \param source_path_prob_map Filepath to the source image's probabilty map image
 * \param source_path_pixel_features Filepath to the source image's pixel features
 * \param dest_path_boundary_features Filepath the boundary features will be written to
 */
void boundary_features(std::string source_path_image,
                       std::string source_ws_image_path,
                       std::string source_path_prob_map,
                       std::string source_path_pixel_features,
                       std::string dest_path_boundary_features)
{
    //IMPORT RESULTS FROM HDF5 file
    //data structur
    //arcs == two cells in cartesian coordinats
    std::vector< std::vector<point> > arcs;
    std::vector<point> junctions;
    std::vector< std::vector<float> > phis;
    std::vector< std::vector<float> > curvs;

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

    std::vector < std::vector<point> > areas;
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
           } 
           areas[ws_region_image(x,y)-1].push_back(p);
        }
    }

    int nr_of_pixel_features;

    //WE OPEN THE IMAGE,AND CONVERT IT INTO AN "LAB"mode and save it into an BRGB Image
    std::string filepath_to_image=source_path_image;
    vigra::ImageImportInfo info_img(filepath_to_image.c_str());
    dim_x=info_img.width();
    dim_y=info_img.height();
    vigra::BRGBImage main_image(dim_x,dim_y);
    vigra::importImage(info_img, destImage(main_image));
    vigra::transformImage(vigra::srcImageRange(main_image),destImage(main_image),vigra::RGBPrime2LabFunctor<float>() );

    //NOW WE HAVE TO OPEN THE PROBABILITY MAP
    vigra::BImage prob_map(dim_x,dim_y);
    std::string filepath_to_prob_map=source_path_prob_map;
    std::string image_filename=get_filename(source_path_image);
    filepath_to_prob_map.append(image_filename);
    filepath_to_image.resize(filepath_to_image.size()-4);
    vigra::ImageImportInfo info_prob(filepath_to_prob_map.c_str());
    vigra::importImage(info_prob, destImage(prob_map));

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
            calculatePhi( &(arcs[i]), &(phis[i]), 15 );
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

    //now we load the feature file
    //OPEN THE FEATURE FILE
    std::cout<<"open pixel-feature file"<<std::endl;
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

    std::cout<<"Nr of features: "<<nr_of_pixel_features<<std::endl;
    
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

    //now we have to open/create the binary file in which the features of one image are stored
    //to do that we have to use the path "dest_path_boundary_features"
    //later we MUST KNOW how many features are in the binary file!
    //  => the first float we store in the binary file HAS TO BE the size of data (=boundary_features*arc.size())

    //to get the the complete filepath  "dest_path_boundary_features" is only the path of the folder

    std::string filepath_to_boundary_features=dest_path_boundary_features;
    filepath_to_boundary_features.append(get_filename(source_path_image));
    filepath_to_boundary_features.append(".bin");

    //open a filepointer

    std::cout<<filepath_to_boundary_features<<std::endl;
    FILE * fp;
    fp =fopen(filepath_to_boundary_features.c_str(),"w+b");
    //first entry in the file should be the arc_size

//IS THIS SAVE???? (to make -1 )    //has been dropped

    float float_arc_size=(float)arcs.size();
    fwrite(&float_arc_size,sizeof(float),1,fp);
    std::cout<<"..done"<<std::endl;

    //now we compute one feature and see what size the  arc_feature_vector is
    //and write this as second entry!

    std::vector<float> arc_feature;
    arc_feature=compute_single_arc_features(arcs,curvs,areas,two_boundings,pixel_features_array,main_image,prob_map,0,dim_x,dim_y,nr_of_pixel_features);
    float float_size=(float)arc_feature.size();

    fwrite(&float_size,sizeof(float),1,fp);

    std::cout<<"compute boundary-features..."<<std::endl;
    //now we compute the feature for each arc and write each value of the feature vector of each arc in a binary file
    for(int i=0;i<(int)arcs.size();i++)
    {
        arc_feature=compute_single_arc_features(arcs,curvs,areas,two_boundings,pixel_features_array,main_image,prob_map,0,dim_x,dim_y,nr_of_pixel_features);

        //std::cout<<"size arg_feature:"<<arc_feature.size()<<std::endl;
        for(int jj=0;jj<(int)arc_feature.size();jj++)
        {
         // std::cout<<arc_feature[i]<<std::endl;
          float value_to_write=arc_feature[jj];
          //write the value into the filepointer
          fwrite(& value_to_write,sizeof(float),1,fp);
        }
    }
    std::cout<<"done"<<std::endl;

    fclose(fp);

    delete pixel_features_array;

    arcs.clear();
}


void create_boundary_classification_set(std::string filepath_to_boundary_training_file,
                                        std::string path_to_boundary_features,              
                                        std::string filepath_to_classification_file,
                                        bool overwrite = false)    
{
    bool found_training_file=false;
    bool found_feature_file=false;

    float  nr_of_arcs=0;
    float  nr_of_features=0;
    float  nr_of_old_features=0;

    //here is some dirty stuff,the last line is readed twice(eof cant look in future),so we define a variable, what the last wirte was
    //and we check if the last read hast already been written, and i use -1 as start value because there cant be a regular arc nr -1
    int last_read_arc=-1;

    //OPEN THE FEATURE FILE:
    //FIRST WE HAVE A LOOK IF THERE IS A BOUNDARY_FEATURE_FILE
    //WITHOUT THIS FEATURE FILE WE CANT DO ANYTHING

    //we have to create the filepath to boundary_features (until now it is only the path to the folder)
    std::string filepath_to_boundary_features=path_to_boundary_features;

    filepath_to_boundary_features.append(get_filename(filepath_to_boundary_training_file));
    filepath_to_boundary_features.resize(filepath_to_boundary_features.size()-3);
    filepath_to_boundary_features.append("bin");

    std::string filepath_to_classification_info=filepath_to_classification_file;
    filepath_to_classification_info.append(".info.bin");
    filepath_to_classification_file.append(".bin");

    FILE *boundary_features_file;
    boundary_features_file=fopen(filepath_to_boundary_features.c_str(),"rb");

    //CHECK IF IS EXISTEND
    if(boundary_features_file==NULL)
    {
        found_feature_file=false;
        std::cout<<"create_boundary_classification_set() could not open "<<filepath_to_boundary_features<<std::endl;
        //exit(-1);
    }
    else
    {
        found_feature_file=true;

        //NOW WE CAN READ OUT THE FIRST 2 FLOATS
        //first float is the number of arcs in this image  =>is stored in nr_of_arcs
        //second float is the number of features           =>is stored in nr_of_features

        int read_in=fread(&nr_of_arcs,sizeof(float),1,boundary_features_file);
        if(read_in!=1)
        {
         std::cout<<"error!  create_boundary_classification() was not able to read in nr_of_arcs"<<std::endl;
         exit(-1);
        }

        read_in=fread(&nr_of_features,sizeof(float),1,boundary_features_file);
        if(read_in!=1)
        {
         std::cout<<"error!  create_boundary_classification() was not able to read in nr_of_features"<<std::endl;
         exit(-1);
        }

        std::cout<<"Nr of Arcs:"<<nr_of_arcs<<" Nr of Boundary Features:"<<nr_of_features<<std::endl;
    }

    if(found_feature_file==true)
    {
        //Now we can create an pointer to an array to store all the features
        float * boundary_features_array=new float[(int)(nr_of_arcs*nr_of_features)];

        //IF Feature file is existend we can read in the boundary features
        int read_in=fread(boundary_features_array,sizeof(float),(int)(nr_of_arcs*nr_of_features),boundary_features_file);
        if(read_in!=nr_of_arcs*nr_of_features)
        {
            std::cout<<"error!  create_boundary_classification() was not able to read in "<<nr_of_arcs*nr_of_features<<" floats"<<std::endl;
            std::cout<<"NR OF FEATURES WE COULD READ IN: "<<read_in<<std::endl;
            exit(-1);
        }
        else
        {
            std::cout<<"feature array has been loaded"<<std::endl;
        }

        fclose(boundary_features_file);

        //Now we know how many arcs are in the image, so we can create and std::vector vector_training_data
        //the index of the vector is the index of the arc, the value to the index is the label of arc
        //if the label=10 this means that there is no training data to this arc

        std::vector<int> vector_training_labels;
        vector_training_labels.resize((int)nr_of_arcs,10);

        std::ifstream training_file(filepath_to_boundary_training_file.c_str());
        if(!training_file)
        {
            found_training_file=false;
        }
        else  //FILE IS EXISTEND
        {
            found_training_file=true;
            std::cout<<"found training data"<<std::endl;
            //FILE IS NOT EMPTY //TODO check this

            int i=0;
            //temp file is used to get one line out of the training file
            std::vector<int> temp;
            temp.resize(2);
            while(!training_file.eof())
            {
                //the vector of the training data is
                training_file>>temp[0];
                training_file>>temp[1];

                if(temp[0]>nr_of_arcs-1)
                {
                    std::cout<<"error in create_boundary_classification_set(..)"<<std::endl;
                    std::cout<<"trainingsdata and image have different nr of arcs:"<<std::endl;
                    std::cout<<"found labeling for arc nr "<<temp[0]<<" in trainingsdata."<<std::endl;
                    exit(-1);
                }
                else if(last_read_arc!=temp[0])
                {
                    if(temp[1]>nr_of_classes-1 || temp[1]<0)
                    {
                        std::cout<<"error in create_boundary_classification_set(..)"<<std::endl;
                        std::cout<<"label "<<temp[1]<<" is not valid. Nr of classes: "<<nr_of_classes<<std::endl;
                        exit(-1);
                    }

                    vector_training_labels[temp[0]]=temp[1];
                    //HERE WE SET LAST READ TO THE LAST READED /WRITTEN POINTS TO AVOID THE LAST LINE
                    //IS READED TWICE
                    last_read_arc=temp[0];
                    i++;
                }
            }
            training_file.close();
            std::cout<<"training data has been loaded"<<std::endl;
        }

        if(found_training_file==false)
        {
            std::cout<<"no training data found in: "<<filepath_to_boundary_training_file<<std::endl;
        }
        else
        {
            std::cout<<"Found feature and training file"<<std::endl;
            int nr_of_old_data=0;
            //we must have a look if there is already an existing classification set! IF there is one we have to read out the fist element =>this should be
            //the nr of training data wich is already stored in the classification set.
            //so what we do:we first try to open the class. set in "rb" mode (means only binary reading)

            //info file
            FILE * info;
            info=fopen(filepath_to_classification_info.c_str(),"rb");

            if(info==NULL)
            {
                std::cout<<"No classification set until now"<<std::endl;
            }
            else
            {
                std::cout<<"Found classification info file "<<filepath_to_classification_info.c_str()<<std::endl;

                //now we must read in how many data is already in the file
                float float_nr_of_old_data;

                int read_in=fread(&float_nr_of_old_data,sizeof(float),1,info);
                if(read_in!=1)
                {
                    std::cout<<"error!  create_boundary_classification() was not able to read in the size of the classification info file"<<std::endl;
                    exit(-1);
                }
                else
                {
                    nr_of_old_data=(int)float_nr_of_old_data;
                }

                read_in=fread(&nr_of_old_features,sizeof(float),1,info);
                if(read_in!=1)
                {
                    std::cout<<"error!  create_boundary_classification() was not able to read in the old nr of features"<<std::endl;
                    exit(-1);
                }
                else if (nr_of_old_features!=nr_of_features)
                {
                    std::cout<<"Nr of old boundary features: "<<nr_of_old_features<<std::endl;
                    std::cout<<"Old classification file contains information with different nr of boundary features"<<std::endl;
                    std::cout<<"Old classification file is deleted."<<std::endl;
                    remove(filepath_to_classification_file.c_str());
                    nr_of_old_data=0;
                }
                else if (overwrite)
                {
                    std::cout<<"Old classification file is deleted."<<std::endl;
                    remove(filepath_to_classification_file.c_str());
                    nr_of_old_data=0;
                }

                fclose(info);
            }

            //now we know how many feature data is already stored in the feature file

            //now we can open the classification set file and write the values into it
            //creates a new file if no file exists
            FILE * classification_set_file;
            classification_set_file=fopen(filepath_to_classification_file.c_str(),"a+b");

            std::vector<int> class_count;
            class_count.resize(10,0);
            for(int i=0;i<(int)vector_training_labels.size();i++)
            {
                //if vector_training_labels[i] ==10 this means that we dont have training
                //data for this arc

                if(vector_training_labels[i]!=10)
                {
                    //vector to store the label and the feature values
                    std::vector<float> label_and_features;

                    //first entry in labels_and_features should be the label
                    label_and_features.push_back((float)vector_training_labels[i]);

                    //we want so know how many data is from wich class
                    class_count[vector_training_labels[i]]=class_count[vector_training_labels[i]]+1;

                    //j is the index the boundary feature array
                    for(int j=(i)*(int)nr_of_features;j<(i+1)*(int)nr_of_features;j++ )
                    {
                        float feature_value=boundary_features_array[j];
                        label_and_features.push_back(feature_value);
                    }

                    //NOW WE GOT THE FEATURES OF ONE TRAINING DATA AND WE CAN WRITE IT TO THE FILE
                    for(int k=0;k<(int)label_and_features.size();k++)
                    {
                        float value=label_and_features[k];
                        fwrite(&value,sizeof(float),1,classification_set_file);
                    }
                }
            }

            //NOW ALL THE TRAINING DATA HAS BEEN ADDED TO THE FILE
            fclose(classification_set_file);

            //Now we musst look know how many data has been added to the file
            int nr_of_new_data=0;

            //nr_of_classes is global const. int
            for(int i=1;i<nr_of_classes;i++)
            {
                if(class_count[i]!=class_count[0])
                {
                    std::cout<<"unbalanced data in "<<get_filename(filepath_to_boundary_training_file)<<std::endl;
                    //exit(-1);
                }
            }

            for(int ii=0;ii<nr_of_classes;ii++)
            {
                std::cout<<"Nr of training data class "<<ii<<": "<<class_count[ii]<<std::endl;
                nr_of_new_data=nr_of_new_data+class_count[ii];
            }

            float all_data=(float)(nr_of_old_data+nr_of_new_data);
            std::cout<<"Nr of OLD training data:"<<nr_of_old_data<<std::endl;
            std::cout<<"Nr of NEW training data:"<<nr_of_new_data<<std::endl;
            std::cout<<"Nr of ALL training data:"<<all_data<<std::endl;

            info=fopen(filepath_to_classification_info.c_str(),"w+b");

            fwrite(&all_data,sizeof(float),1,info);
            fwrite(&nr_of_features,sizeof(float),1,info);

            fclose(info);

            //and we can delete the boundary feature array
            delete boundary_features_array;
        }

    }//endif

    else
    {
        std::cout<<"nothing to do here"<<std::endl;
    }
}
