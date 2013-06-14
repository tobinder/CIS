/*! \file pixel_probabilities.h
 *  \brief Pixel probability calculations.
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

#include <stdio.h>
#include <omp.h>

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <math.h>
#include <sys/time.h>

#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/impex.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/random_forest.hxx>

#include <vigra/hdf5impex.hxx>
#include <vigra/random_forest_hdf5_impex.hxx>
#include <vigra/hdf5impex.hxx>

#include <vigra/matrix.hxx>

/*! \fn train_pixel_rf(std::string filepath_to_classification_set_file,std::string dest_filepath_to_random_forest_file,std::string param_file)
 * \brief Train a pixel random forest.
 * \param filepath_to_classification_set Filepath to the classification set
 * \param dest_filepath_to_random_forest_file Filepath to the target random forest file
 * \param param_file Name of the parameter file
 */
void train_pixel_rf(std::string filepath_to_classification_set_file,std::string dest_filepath_to_random_forest_file,std::string param_file)
{
    ParameterFile paramFile;

    if( !paramFile.load(param_file) )
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }

    Parameter<int> nr_of_trees;
    nr_of_trees.assign("", "nr_of_trees", 256);
    nr_of_trees.load(paramFile,"config");

    int nr_of_training_data_0;
    int nr_of_training_data_1;
    int nr_of_features;
    int nr_of_training_data;

    //first we load the *.info.bin file to see how many training data is in the *.0.bin and *.1.bin
    //and how many features have been used

    FILE *info;
    std::string filepath_to_classification_set_file_info=filepath_to_classification_set_file;
    filepath_to_classification_set_file_info.append(".info.bin");
    info =fopen(filepath_to_classification_set_file_info.c_str(),"rb");
    if(info==NULL)  //if the file is NOT existend
    {
        std::cout<<filepath_to_classification_set_file_info<<" is NOT existend"<<std::endl;
        exit(-1);
    }
    else        //if it is existend
    {

      int read_out=0;

      read_out=read_out+fread(&nr_of_training_data_0,sizeof(int),1,info);
      read_out=read_out+fread(&nr_of_training_data_1,sizeof(int),1 ,info);
      read_out=read_out+fread(&nr_of_features,sizeof(int),1 ,info);

      nr_of_training_data =nr_of_training_data_0+nr_of_training_data_1;
      fclose(info);

      std::cout<<"Train Random Forest with:"<<std::endl;
      std::cout<<"Nr of Training data with label=0: "<<nr_of_training_data_0<<std::endl;
      std::cout<<"Nr of Training data with label=1: "<<nr_of_training_data_1<<std::endl;
      std::cout<<"Nr of Features: "<<nr_of_features<<std::endl;

      if(read_out!=3)
      {
            std::cout<<filepath_to_classification_set_file_info<<" is wrong!"<<std::endl;
            exit(-1);
      }
    }

    //some filepointers to get the training data
    FILE *fp0;
    FILE *fp1;
    std::string filepath_to_classification_set_file_0=filepath_to_classification_set_file;
    std::string filepath_to_classification_set_file_1=filepath_to_classification_set_file;
    fp0 =fopen(filepath_to_classification_set_file_0.append(".0.bin").c_str(),"rb");   //rb means read (only) binary,no writing!
    fp1 =fopen(filepath_to_classification_set_file_1.append(".1.bin").c_str(),"rb");   //rb means read (only) binary,no writing!

    int nr_to_use=std::min(nr_of_training_data_0, nr_of_training_data_0);
//    int nr_to_use;
//    std::cout<<"Nr of Training data to be used: ";
//  std::cin>>nr_to_use;
    std::cout<<"Used Training data: "<<nr_to_use<<std::endl;

    nr_of_training_data_0=nr_to_use;
    nr_of_training_data_1=nr_to_use;
    nr_of_training_data=nr_to_use*2;

    float * feature_array_0=new float[(nr_of_training_data_0*nr_of_features)];
    float * feature_array_1=new float[(nr_of_training_data_1*nr_of_features)];
    float * feature_array=new float[(nr_of_training_data_0*nr_of_features+nr_of_training_data_1*nr_of_features)];
      int *   label_array=new int[(nr_of_training_data_0+nr_of_training_data_1)];

    //if these 2 files are existend we can read them out
    if(fp0!=NULL && fp1!=NULL)
    {
        //we read in the floats
        int read_in_0=fread(feature_array_0,sizeof(float),nr_of_features*nr_of_training_data_0 ,fp0);
        //check if all floats habe been loaded/readed
        if(read_in_0!=nr_of_training_data_0*nr_of_features)
        {
            std::cout<<"Error: "<<filepath_to_classification_set_file_0<<" has NOT been readed in correct "<<std::endl;
            exit(-1);
        }
        std::cout<<filepath_to_classification_set_file_0<<" has been loaded"<<std::endl;
        //we read in the floats
        int read_in_1=fread(feature_array_1,sizeof(float),nr_of_features*nr_of_training_data_1 ,fp1);
        //check if all floats habe been loaded/readed

        if(read_in_1!=nr_of_training_data_1*nr_of_features)
        {
            std::cout<<"Error: "<<filepath_to_classification_set_file_1<<" has NOT been readed in correct "<<std::endl;
            exit(-1);
        }
        std::cout<<filepath_to_classification_set_file_1<<" has been loaded"<<std::endl;

        //merge part 1
        for(int i=0;i<nr_of_training_data_0*nr_of_features;i++ )
        {
            feature_array[i]=feature_array_0[i];
        }
        //merge part 2
        for(int i=0;i<nr_of_training_data_1*nr_of_features;i++ )
        {
            feature_array[i+nr_of_training_data_0*nr_of_features]=feature_array_1[i];
        }

        //now we have to insert the labes into the label_array
        for(int i=0;i<nr_of_training_data;i++)
        {
            if(i<nr_of_training_data_0)
            {
                label_array[i]=0;
            }
            if(i>=nr_of_training_data_0)
            {
                label_array[i]=1;
            }
        }
    }

    if(fp0==NULL || fp1==NULL)
    {
    std::cout<<filepath_to_classification_set_file_0<<" or "<<filepath_to_classification_set_file_1<<"is NOT existend"<<std::endl;
    exit(-1);
    }

    //CLosing files
    fclose(fp0);
    fclose(fp1);
    delete feature_array_0;
    delete feature_array_1;
    //Delete feature array_0 and feature_array_1

    vigra::MultiArray<2, float> features(vigra::MultiArrayShape<2>::type(nr_of_training_data, nr_of_features),feature_array);
    for(int x=0;x<nr_of_training_data;x++)
    {
        for(int y=0;y<nr_of_features;y++)
        {
        features(x,y)=feature_array[nr_of_features*x+y];
        }
    }

    vigra::MultiArray<2, int> labels(vigra::MultiArrayShape<2>::type(nr_of_training_data, 1),label_array);

    vigra::RandomForest<> rf(vigra::RandomForestOptions().tree_count(nr_of_trees));

    std::cout<<"Train Random Forest... "<<std::endl;
    size_t numNans = 0;
    size_t numInfs = 0;
    for(int i=0; i<features.shape(0); ++i) {
        for(int j=0; j<features.shape(1); ++j) {
            if(std::isnan(features(i,j))) {
                ++numNans;
                std::stringstream s; s << "feature(" << i << "," << j << ") = " << features(i,j);
                std::cout << s.str() << std::endl;
                //throw std::runtime_error(s.str().c_str());
                features(i,j) = 0;
            }
            if(std::isinf(features(i,j))) {
                ++numInfs;
                std::stringstream s; s << "feature(" << i << "," << j << ") = " << features(i,j);
                std::cout << s.str() << std::endl;
                //throw std::runtime_error("isinf");
            }
        }
    }
    std::cout << "numNans= " << numNans << ", numInfs = " << numInfs << std::endl;
    //return;

    for(int i=0; i<labels.shape(0); ++i) {
        if(!(labels(i) == 0 || labels(i) == 1)) { throw std::runtime_error("not in 01"); }
    }

    //vigra::rf::visitors::VariableImportanceVisitor var_imp;
    //vigra::rf::visitors::OOB_Error oob;

    rf.learn(features,labels);
    //rf.learn(features,labels, create_visitor(var_imp));
    //rf.learn(features,labels, create_visitor(oob));
    //std::cout<<"oob error: "<<oob.oob_breiman<<std::endl;
    std::cout<<"...training done "<<std::endl;

    //create an variable importance output file
    //std::string filepath_to_variable_importance=dest_filepath_to_random_forest_file;
    //filepath_to_variable_importance.append("_var.hdf5");
    //var_imp.variable_importance_;
    //var_imp.save(filepath_to_variable_importance, get_filename(filepath_to_variable_importance));

    //EXPORT THE RANDOM FORREST
    std::string dest_filepath_to_random_forest_info_file=dest_filepath_to_random_forest_file;
    dest_filepath_to_random_forest_file.append("_rf.hdf5");

    //if the features have changed the old file "_rf.hdf5" can cause errors, so delete it in any case!
    remove(dest_filepath_to_random_forest_file.c_str());
    std::cout<<"Old Random Forest File deleted."<<std::endl;

    std::cout<<"Export the Random Forest..."<<std::endl;
    rf_export_HDF5(rf, dest_filepath_to_random_forest_file.c_str());
    std::cout<<"...done "<<std::endl;

    //EXPORT an other "info" file,in this info file is written how many
    //features has been used to train the random forest
    FILE *info_rf;
    dest_filepath_to_random_forest_info_file.append(".info.bin");
    info_rf =fopen(dest_filepath_to_random_forest_info_file.c_str(),"wb");

    fwrite(&nr_of_features,sizeof(int),1,info_rf);
    fclose(info_rf);

    //delete the arrays
    delete feature_array;
    delete label_array;
}

/*! \fn extract_pixel_probabilities(std::string fn,std::string path_to_feature_files,std::string filepath_to_random_forest_file,std::string dest_path_probability_maps,std::string param_file)
 * \brief Extract pixel probabilities.
 * \param fn Filepath to the image file
 * \param path_to_feature_files Path to the feature file
 * \param filepath_to_random_forest_file Filepath to the random forest file
 * \param dest_path_probability_maps Destination path of the probability maps
 * \param param_file Name of the parameter file
 */
//fn is the filepath to the image file
void extract_pixel_probabilities(std::string fn,std::string path_to_feature_files,std::string filepath_to_random_forest_file,std::string dest_path_probability_maps,std::string param_file)
{
    int nr_of_features;
    //First we open the random_forest file and the little info file wich is saved in the same folder as the random_forest file
    //The info file just contains the information how many features has been used to train the Random forest
    std::string filepath_to_random_forest_info_file=filepath_to_random_forest_file;
    filepath_to_random_forest_file.append("_rf.hdf5");
    filepath_to_random_forest_info_file.append(".info.bin");

    FILE *info;
    info =fopen(filepath_to_random_forest_info_file.c_str(),"rb");

    if(info==NULL) //file is not existend
    {
        std::cout<<"Random Forest Info File "<<filepath_to_random_forest_info_file<<" is NOT existend"<<std::endl;
        exit(-1);
    }
    else     //file is existend
    {
        int read_out=0;
        read_out=read_out+fread(&nr_of_features,sizeof(int),1 ,info);
        if(read_out!=1)
        {
            std::cout<<filepath_to_random_forest_info_file<<" is wrong!"<<std::endl;
            exit(-1);
        }
    }

    fclose (info);    

    vigra::RandomForest<> rf;

    rf_import_HDF5(rf,filepath_to_random_forest_file.c_str());

    //NOW FOR UNKNOWN IMAGES/FEATURES
    vigra::ImageImportInfo info_image(fn.c_str());

    int dim_x=info_image.width();
    int dim_y=info_image.height();
    int nr_of_pixel=dim_x*dim_y;

    std::cout<<"Image:"<<fn<<std::endl;
    std::cout<<"nr_of_pixel:"<<nr_of_pixel<<std::endl;

    //Check if there is a bubble image
    std::string fn_bubbles = fn;
    fn_bubbles.resize(fn_bubbles.size()-3);
    fn_bubbles.append("b.bmp");
    
    FILE *info_bubbles;
    info_bubbles = fopen(fn_bubbles.c_str(), "rb");

    //Load bubble image
    vigra::BasicImage<bool> bubbles(dim_x, dim_y);

    bool bubble_img = false;

    if(info_bubbles == NULL)
    {
        std::cout << "No bubble image found" << std::endl;
        //Set the bubble image to black then
        for(int y = 0; y < dim_y; y++)
        {
            for(int x = 0; x < dim_x; x++)
            {
                bubbles(x,y) = false;
            }
        }
    }
    else
    {
        std::cout << "Bubble image found" << std::endl;
        vigra::ImageImportInfo info_b_image(fn_bubbles.c_str());
        importImage(info_b_image, destImage(bubbles));
        fclose(info_bubbles);
        bubble_img = true;
    }

    //if there is a threshold defined we set all pixels above the threshold to high
    //no boundary-probability
    ParameterFile paramFile;

    if( !paramFile.load(param_file) )
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    //threshold 0 means without threshold
    //threshold below 10 means calculate threshold with given factor
    Parameter<int> gray_threshold;
    gray_threshold.assign("", "gray_threshold", 0);
    gray_threshold.load(paramFile,"config");
/*
    //Add 2 to the gray threshold if there exists a bubble image
    if(bubble_img)
    {                
        gray_threshold = gray_threshold + 2;
    }
*/
    //std::cout << "Parameter gray_threshold: " << gray_threshold << std::endl;

    Parameter<int> light_threshold;
    light_threshold.assign("", "light_threshold", 0);
    light_threshold.load(paramFile,"config");

    //std::cout << "Parameter light_threshold: " << light_threshold << std::endl;

    //threshold 0 means without threshold, feature file not loaded
    Parameter<float> feature_threshold;
    feature_threshold.assign("", "feature_threshold", 0.0);
    feature_threshold.load(paramFile,"config");

    //std::cout << "Parameter feature_threshold: " << feature_threshold << std::endl;

    Parameter<float> no_boundary_probability;
    no_boundary_probability.assign("", "no_boundary_probability", 0.95f);
    no_boundary_probability.load(paramFile,"config");

    //std::cout << "Parameter no_boundary_probability: " << no_boundary_probability << std::endl;

    //OPEN THE FEATURE FILE
    std::string only_filename_of_the_image=get_filename(fn);
    std::string filepath_to_feature_files=path_to_feature_files;
    filepath_to_feature_files.append(only_filename_of_the_image);
    filepath_to_feature_files.append(".bin");

    FILE *fp_unknown_features;
    fp_unknown_features =fopen(filepath_to_feature_files.c_str(),"rb");
    //create an array
    //float * unknown_features_array=new float[nr_of_pixel*nr_of_features];
    float **unknown_features_array = new float*[nr_of_features];
    for(int feature=0;feature<nr_of_features;feature++)
    {
        unknown_features_array[feature] = new float[nr_of_pixel];
    }

    if(fp_unknown_features==NULL)
    {
        std::cout<<"Error: Feature File "<<filepath_to_feature_files<<" is NOT existend"<<std::endl;
        exit(-1);
    }
    else
    {
        for (int x=0; x<dim_x; x++)
        {
            float * array = new float[dim_y*nr_of_features];

            int read_in=fread(array,sizeof(float),dim_y*nr_of_features,fp_unknown_features);
            if(read_in!=dim_y*nr_of_features)
            {
                std::cout<<" extract_pixel_probabilities() error, in unknown features"<<std::endl;
                exit(-1);
            }

            for (int y=0;y<dim_y;y++)
            {
                int z=dim_y*x+y;

                for(int feature=0;feature<nr_of_features;feature++)
                {
                    //unknown_features_array[nr_of_features*z+feature]=array[nr_of_features*y+feature];
                    unknown_features_array[feature][z]=array[nr_of_features*y+feature];
                }
            }

            delete array;
        }
    
        //for(size_t jj=0;jj<nr_of_pixel*nr_of_features;jj++)
        for(size_t jj=0;jj<nr_of_pixel;jj++)
        {
            //if (std::isnan(unknown_features_array[jj]))
            //{
            //    unknown_features_array[jj] = 0;
            //}
            for(int feature=0;feature<nr_of_features;feature++)
            {
                if (std::isnan(unknown_features_array[feature][jj]))
                {
                    unknown_features_array[feature][jj] = 0;
                }
            }
        }
    }

    fclose(fp_unknown_features);

    std::cout<<"Probability prediction...."<<std::endl;

    timeval start, end;
    gettimeofday(&start, 0);

    dest_path_probability_maps.append(param_file_name);
    dest_path_probability_maps.append(only_filename_of_the_image);
    //dest_path_probability_maps.append(".0.png");

    if (gray_threshold>0 || feature_threshold>0.0f)//threshold
    {
        //OPEN THE IMAGE
        vigra::FImage image;
        image.resize(dim_x,dim_y);
        importImage(info_image, destImage(image));
        
        vigra::FImage probability_image_0(dim_x,dim_y);
        vigra::IImage ws_label(dim_x,dim_y);

        vigra::BasicImage<bool> temp_img(dim_x,dim_y);
        for(int x=0;x<dim_x;x++)
        {
            for (int y=0;y<dim_y;y++)
            {
                temp_img(x,y)=false;
                probability_image_0(x,y)=1.0f;
            }
        }

        if (gray_threshold>0)//gray value threshold
        {
            if (gray_threshold<=10)//calculate threshold
            {
                std::cout<<"Calculate grayvalue threshold with factor "<<gray_threshold<<std::endl;

                //load NEEM stitching information
                std::ifstream stitch_info_file("/export/disk01/tobinder/Data/NEEM-LASM/List_StitchedImages.txt");
                std::ifstream temp_stitch_info_file("/export/disk01/tobinder/Data/NEEM-LASM/List_StitchedImages.txt");
                int pos_right=dim_x;

                //string is just for testing stuff
                std::string teststring;
                temp_stitch_info_file>>teststring;
                if(stitch_info_file) //file exists
                {
                    if(teststring.size()!=0)//file is not empty
                    {
                        std::string name1, name2, path_stitched;
                        int x, y;
                        bool found=false;

                        while(!stitch_info_file.eof() && !found)
                        {
                            stitch_info_file>>name1>>name2>>x>>y>>path_stitched;
                            if (get_filename(path_stitched)==only_filename_of_the_image)
                                {
                                std::cout<<"Found stitching information: "<<x<<" "<<y<<std::endl;
                                pos_right=4021+x/2;
                                std::cout<<"Right image starts at: "<<pos_right<<std::endl;
                                found=true;
                                }
                        }
                    }

                    stitch_info_file.close();
                    temp_stitch_info_file.close();
                }

                if (pos_right==dim_x)
                {
                    float grayvalue_histogram[256];
                    for (int i=0;i<256;i++) grayvalue_histogram[i]=0;
             
                    for (int y=0;y<dim_y;y++)
                    {
                        for (int x=0;x<dim_x;x++)
                        {
                            int value=image(x,y);
                            grayvalue_histogram[value]++;
                        }
                    }
             
                    for (int i=0;i<256;i++)
                    {
                        grayvalue_histogram[i]=grayvalue_histogram[i]/(float)(dim_x*dim_y/1000);
             
                        float average=0.0f;
                        for (int j=0;j<=i;j++) average+=(grayvalue_histogram[j]/(i+1));
             
                        if (grayvalue_histogram[i]>average*(float)gray_threshold && grayvalue_histogram[i]>gray_threshold && i>50)
                        {
                            gray_threshold=i;
                            break;
                        }
                    }
             
                    if (gray_threshold<50)
                    {
                        std::cout<<"no gray value threshold found, set threshold to 50"<<std::endl;
                        gray_threshold=50;
                    }
                    else std::cout<<"gray value threshold calculated to: "<<gray_threshold<<std::endl;  
                }
                else
                {
                    int initial=gray_threshold;

                    float grayvalue_histogram[256];
                    for (int i=0;i<256;i++) grayvalue_histogram[i]=0;
             
                    for (int y=0;y<dim_y;y++)
                    {
                        for (int x=0;x<pos_right;x++)
                        {
                            int value=image(x,y);
                            grayvalue_histogram[value]++;
                        }
                    }
             
                    for (int i=0;i<256;i++)
                    {
                        grayvalue_histogram[i]=grayvalue_histogram[i]/(float)(dim_x*dim_y/1000);
             
                        float average=0.0f;
                        for (int j=0;j<=i;j++) average+=(grayvalue_histogram[j]/(i+1));
             
                        if (grayvalue_histogram[i]>average*(float)gray_threshold && grayvalue_histogram[i]>gray_threshold && i>50)
                        {
                            gray_threshold=i;
                            break;
                        }
                    }
             
                    std::cout<<"left gray value threshold calculated to: "<<gray_threshold<<std::endl;  

                    int left_threshold=gray_threshold;
                    gray_threshold=initial;

                    for (int i=0;i<256;i++) grayvalue_histogram[i]=0;
             
                    for (int y=0;y<dim_y;y++)
                    {
                        for (int x=pos_right;x<dim_x;x++)
                        {
                            int value=image(x,y);
                            grayvalue_histogram[value]++;
                        }
                    }
             
                    for (int i=0;i<256;i++)
                    {
                        grayvalue_histogram[i]=grayvalue_histogram[i]/(float)(dim_x*dim_y/1000);
             
                        float average=0.0f;
                        for (int j=0;j<=i;j++) average+=(grayvalue_histogram[j]/(i+1));
             
                        if (grayvalue_histogram[i]>average*(float)gray_threshold && grayvalue_histogram[i]>gray_threshold && i>50)
                        {
                            gray_threshold=i;
                            break;
                        }
                    }
             
                    std::cout<<"right gray value threshold calculated to: "<<gray_threshold<<std::endl;  
                    if (left_threshold>gray_threshold) gray_threshold=left_threshold;
                    std::cout<<"maximal gray value threshold: "<<gray_threshold<<std::endl;  
                }
            }
            else std::cout<<"use defined gray value threshold: "<<gray_threshold<<std::endl;  

            if (light_threshold==0) light_threshold=256;//no light threshold
            else std::cout<<"use defined light threshold: "<<light_threshold<<std::endl;  

            for(int x=0;x<dim_x;x++)
            {
                for (int y=0;y<dim_y;y++)
                {
                    if (image(x,y)<gray_threshold || image(x,y)>light_threshold)
                    {
                        temp_img(x,y)=true;
                    }
                }
            }
        }

        // due to morphologic closing some pixels are filtered even if both thresholds overlap.
        // to guarantuee that no pixels are filtered uncomment here
        if (gray_threshold>=light_threshold)
        {
            //gray_threshold=0;
        }

        if (feature_threshold>0.0f)//threshold
        {
            //OPEN THE DECISION FEATURE FILE
            filepath_to_feature_files.resize(filepath_to_feature_files.size()-4);
            filepath_to_feature_files.append(".decision.bin");

            FILE *fp_decision_feature;
            fp_decision_feature =fopen(filepath_to_feature_files.c_str(),"rb");
            //create an array
            float * decision_feature_array=new float[nr_of_pixel];

            if(fp_decision_feature==NULL)
            {
                std::cout<<"Error: Decision Feature File "<<filepath_to_feature_files<<" is NOT existend"<<std::endl;
                exit(-1);
            }
            else
            {
                int read_in=fread(decision_feature_array,sizeof(float),nr_of_pixel,fp_decision_feature);
                if(read_in!=nr_of_pixel)
                {
                    std::cout<<" extract_pixel_probabilities() error, in unknown features"<<std::endl;
                    exit(-1);
                }
            }
    
            for(size_t jj=0;jj<nr_of_pixel;jj++)
            {
                if (std::isnan(decision_feature_array[jj]))
                {
                    decision_feature_array[jj] = 0;
                }
            }

            fclose(fp_decision_feature);
            std::cout<<"use defined feature threshold: "<<feature_threshold<<std::endl;  

            for(int x=0;x<dim_x;x++)
            {
                for (int y=0;y<dim_y;y++)
                {
                    int z=dim_x*y+x;

                    if (decision_feature_array[z]>feature_threshold)
                    {
                        temp_img(x,y)=true;
                    }
                }
            }

            delete decision_feature_array;
        }

        if (gray_threshold>0)
        {
            //exportImage(srcImageRange(temp_img), vigra::ImageExportInfo("filtered_pixels1.bmp"));

            float LoGSigma=1.0f;
            float sigma=0.1f;
            int r=5;

            std::cout<<"calculate laplace of gaussian"<<std::endl;

            //We want to do a morphological closing operation and need two float buffer-images for that
            vigra::FImage temp_img_float2(dim_x,dim_y);

            //Do a laplacian of gaussian to highligth gradients
            vigra::laplacianOfGaussian(srcImageRange(image), destImage(temp_img_float2), LoGSigma);

            //If there exists a bubble image, set outlines of bubble areas to true and inside to false
            if(bubble_img == true)
            {
                for(int x = 0; x < dim_x; x++)
                {
                    for(int y = 0; y < dim_y; y++)
                    {
                        if(bubbles(x,y) == true)
                        {                
                            temp_img(x,y) = false;
                        }
                    }
                }
                
                //Traverse along x-axis first
                for(int x = 0; x < dim_x; x++)
                {
                    for(int y = 0; y < dim_y; y++)
                    {
                        if(bubbles(x,y) == false && x+1 < dim_x)
                        {
                            if(bubbles(x+1,y) == true)
                            {
                                temp_img(x+1,y) = true;
                                temp_img_float2(x+1,y) = 5.0f;
                                probability_image_0(x+1,y) = 0.0f;
                                if(x+2 < dim_x && y+1 < dim_y)
                                {
                                    temp_img(x+2,y) = true;
                                    temp_img(x+1,y+1) = true;
                                    temp_img_float2(x+2,y) = 5.0f;
                                    temp_img_float2(x+1,y+1) = 5.0f;
                                    probability_image_0(x+2,y) = 0.0f;
                                    probability_image_0(x+1,y+1) = 0.0f;
                                }
                                if(y-1 >= 0)
                                {
                                    temp_img(x,y) = true;
                                    temp_img(x+1,y-1) = true;
                                    temp_img_float2(x,y) = 5.0f;
                                    temp_img_float2(x+1,y-1) = 5.0f;
                                    probability_image_0(x,y) = 0.0f;
                                    probability_image_0(x+1,y-1) = 0.0f;
                                }
                            }
                        }
                        if(bubbles(x,y) == true && x+1 < dim_x)
                        {
                            if(bubbles(x+1,y) == false)
                            {
                                temp_img(x,y) = true;
                                temp_img_float2(x,y) = 5.0f;
                                probability_image_0(x,y) = 0.0f;
                                if(x+1 < dim_x && y+1 < dim_y)
                                {
                                    temp_img(x+1,y) = true;
                                    temp_img(x,y+1) = true;
                                    temp_img_float2(x+1,y) = 5.0f;
                                    temp_img_float2(x,y+1) = 5.0f;
                                    probability_image_0(x+1,y) = 0.0f;
                                    probability_image_0(x,y+1) = 0.0f;
                                }
                                if(x-1 >= 0 && y-1 >= 0)
                                {
                                    temp_img(x-1,y) = true;
                                    temp_img(x,y-1) = true;
                                    temp_img_float2(x-1,y) = 5.0f;
                                    temp_img_float2(x,y-1) = 5.0f;
                                    probability_image_0(x-1,y) = 0.0f;
                                    probability_image_0(x,y-1) = 0.0f;
                                }
                            }
                        }
                    }
                }

                //Traverse along y-axis 
                for(int y = 0; y < dim_y; y++)
                {
                    for(int x = 0; x < dim_x; x++)
                    {
                        if(bubbles(x,y) == false && y+1 < dim_y)
                        {
                            if(bubbles(x,y+1) == true)
                            {
                                temp_img(x,y+1) = true;
                                temp_img_float2(x,y+1) = 5.0f;
                                probability_image_0(x,y+1) = 0.0f;
                                if(x+1 < dim_x && y+2 < dim_y)
                                {
                                    temp_img(x+1,y+1) = true;
                                    temp_img(x,y+2) = true;
                                    temp_img_float2(x+1,y+1) = 5.0f;
                                    temp_img_float2(x,y+2) = 5.0f;
                                    probability_image_0(x+1,y+1) = 0.0f;
                                    probability_image_0(x,y+2) = 0.0f;
                                }
                                if(y-1 >= 0)
                                {
                                    temp_img(x-1,y+1) = true;
                                    temp_img(x,y) = true;
                                    temp_img_float2(x-1,y+1) = 5.0f;
                                    temp_img_float2(x,y) = 5.0f;
                                    probability_image_0(x-1,y+1) = 0.0f;
                                    probability_image_0(x,y) = 0.0f;
                                }
                            }
                        }
                        if(bubbles(x,y) == true && y+1 < dim_y)
                        {
                            if(bubbles(x,y+1) == false)
                            {
                                temp_img(x,y) = true;
                                temp_img_float2(x,y) = 5.0f;
                                probability_image_0(x,y) = 0.0f;
                                if(x+1 < dim_x && y+2 < dim_y)
                                {
                                    temp_img(x+1,y+1) = true;
                                    temp_img(x,y+2) = true;
                                    temp_img_float2(x+1,y+1) = 5.0f;
                                    temp_img_float2(x,y+2) = 5.0f;
                                    probability_image_0(x+1,y+1) = 0.0f;
                                    probability_image_0(x,y+2) = 0.0f;
                                }
                                if(y-1 >= 0)
                                {
                                    temp_img(x-1,y+1) = true;
                                    temp_img(x,y) = true;
                                    temp_img_float2(x-1,y+1) = 5.0f;
                                    temp_img_float2(x,y) = 5.0f;
                                    probability_image_0(x-1,y+1) = 0.0f;
                                    probability_image_0(x,y) = 0.0f;
                                }
                            }
                        }
                    }
                }
                //exportImage(srcImageRange(temp_img), vigra::ImageExportInfo("filtered_pixels_bubble.bmp"));
            }

            std::cout<<"find connected filtered pixels"<<std::endl;

            vigra::IImage labels(dim_x,dim_y);
            int nr_labels;

            // find connected regions
            nr_labels=vigra::labelImageWithBackground(srcImageRange(temp_img), destImage(labels), true, 0);
            
            int * label_size=new int[nr_labels+1];

            for(int i=0; i<nr_labels+1; i++) label_size[i]=0;

            for(int x=0;x<dim_x;x++)
            {
                for (int y=0;y<dim_y;y++)
                {
                    label_size[labels(x,y)]++;
                }
            }

            for(int x=0;x<dim_x;x++)
            {
                for (int y=0;y<dim_y;y++)
                {
                    if (label_size[labels(x,y)]<1000) temp_img(x,y)=false;
                }
            }

            delete label_size;

            //exportImage(srcImageRange(temp_img), vigra::ImageExportInfo("filtered_pixels2.bmp"));

            std::cout<<"close filtered pixels"<<std::endl;

            //Copy from buffer
            for(int i=0; i<temp_img.width(); i++)
                for(int j=0; j<temp_img.height(); j++)
                {
                    if (temp_img(i,j)==false) temp_img_float2(i,j) = 0;
                }

            vigra::FImage temp_img_float1(temp_img_float2);
            //exportImage(srcImageRange(temp_img_float2), vigra::ImageExportInfo("filtered_pixels_laplace.bmp"));

            // Disc Dilation with radius r
            for(int i=0; i<temp_img.width(); i++)
                for(int j=0; j<temp_img.height(); j++)
                    for(int k = -r; k<=r; k++)
                        for(int l= -r; l<=r; l++)
                        {
                            if( k*k + l*l <= r*r )
                            {
                                if( k+i >= 0 && k+i < temp_img_float2.width() && j+l >=0 && j+l < temp_img_float2.height() )
                                {
                                    if( temp_img_float2( k+i, j+l ) > temp_img_float1( i, j ) )
                                        temp_img_float1( i, j) = temp_img_float2( k+i, j+l );
                                }
                            }
                        }

            //Copy from buffer
            for(int i=0; i<temp_img.width(); i++)
                for(int j=0; j<temp_img.height(); j++)
                    temp_img_float2(i,j) = temp_img_float1(i,j);

            // Disc Erosion with radius r
            for(int i=0; i<temp_img.width(); i++)
                for(int j=0; j<temp_img.height(); j++)
                    for(int k = -r; k<=r; k++)
                        for(int l= -r; l<=r; l++)
                        {
                            if( k*k + l*l <= r*r )
                            {
                                if( k+i >= 0 && k+i < temp_img_float2.width() && j+l >=0 && j+l < temp_img_float2.height() )
                                {
                                    if( temp_img_float2( k+i, j+l ) < temp_img_float1( i, j ) )
                                        temp_img_float1( i, j) = temp_img_float2( k+i, j+l );
                                }
                            }
                        }

            //Copy from buffer
            for(int i=0; i<temp_img.width(); i++)
                for(int j=0; j<temp_img.height(); j++)
                {
                    if (temp_img_float1(i,j)>0.0f) temp_img(i,j) = true;
                    else temp_img(i,j) = false;
                    temp_img_float1(i,j)=0.0f;
                }
        }

        //exportImage(srcImageRange(temp_img), vigra::ImageExportInfo("filtered_pixels3.bmp"));
        /*
        vigra::IImage labels_(dim_x,dim_y);
        int nr_labels_;

        // find connected regions
        nr_labels_=vigra::labelImageWithBackground(srcImageRange(temp_img), destImage(labels_), true, 0);

        int * label_size_=new int[nr_labels_+1];

        for(int i=0; i<nr_labels_+1; i++) label_size_[i]=0;

        for(int x=0;x<dim_x;x++)
        {
            for (int y=0;y<dim_y;y++)
            {
                label_size_[labels_(x,y)]++;
            }
        }

        for(int x=0;x<dim_x;x++)
        {
            for (int y=0;y<dim_y;y++)
            {
                if (label_size_[labels_(x,y)]<1000) temp_img(x,y)=false;
            }
        }

        delete label_size_;

        for(int i=0; i<temp_img.width(); i++)
            for(int j=0; j<temp_img.height(); j++)
            {
                bool not_marked_found=false;
                for(int k = -3*r; k<=3*r; k++)
                    for(int l= -3*r; l<=3*r; l++)
                    {
                        if( k*k + l*l > r*r )
                        {
                            if( k+i >= 0 && k+i < temp_img_float2.width() && j+l >=0 && j+l < temp_img_float2.height() )
                            {
                                if(temp_img(k+i,j+l) == false) not_marked_found=true;
                            }
                        }
                    }
                if (!not_marked_found) temp_img_float1(i,j)=1.0f;               
            }

        //Copy from buffer
        for(int i=0; i<temp_img.width(); i++)
            for(int j=0; j<temp_img.height(); j++)
            {
                if (temp_img_float1(i,j)==1.0f) temp_img(i,j) = false;
            }

        exportImage(srcImageRange(temp_img), vigra::ImageExportInfo("filtered_pixels4.bmp"));
        */
        /*
        for(int x=0;x<dim_x;x++)
        {
            for (int y=0;y<dim_y;y++)
            {
                if (temp_img(x,y)==false)
                {
                    image(x,y)=255;
                }
            }
        }
        */
        //exportImage(srcImageRange(image), vigra::ImageExportInfo("thresholded_image.bmp"));

        std::cout<<"start parallelized prediction"<<std::endl;

        #pragma omp parallel
        {
            float * pixel_unknown_features_array=new float[nr_of_features];
            vigra::MultiArray<2, float> pixel_unknown_features(vigra::MultiArrayShape<2>::type(1, nr_of_features),pixel_unknown_features_array);

            #pragma omp for
            for(int x=0;x<dim_x;x++)
            {
                for (int y=0;y<dim_y;y++)
                {
                    int z=dim_x*y+x;

                    if (temp_img(x,y)==true && probability_image_0(x,y)==0.0f) ws_label(x,y)=0;
                    else if (temp_img(x,y)==true)
                    {
                        for(int feature=0;feature<nr_of_features;feature++)
                        {
                            //pixel_unknown_features(0,feature)=unknown_features_array[nr_of_features*z+feature];
                            pixel_unknown_features(0,feature)=unknown_features_array[feature][z];
                        }

                        vigra::MultiArray<2, double> pixel_unknown_probability(vigra::MultiArrayShape<2>::type(1,2));

                        rf.predictProbabilities(pixel_unknown_features,pixel_unknown_probability);

                        probability_image_0(x,y)=(float)pixel_unknown_probability(0,0);
                        ws_label(x,y)=0;
                    }
                    else
                    {
                        probability_image_0(x,y)=no_boundary_probability;
                        ws_label(x,y)=1;
                    }
                }
            }

            delete pixel_unknown_features_array;  

        }//end of parallelisation

        std::cout<<"end parallelized prediction"<<std::endl;
        
        exportImage(srcImageRange(probability_image_0), vigra::ImageExportInfo(dest_path_probability_maps.c_str()));
        dest_path_probability_maps.resize(dest_path_probability_maps.size()-3);
        dest_path_probability_maps.append("0.bmp");
        if (gray_threshold>0) exportImage(srcImageRange(ws_label), vigra::ImageExportInfo(dest_path_probability_maps.c_str()));
    }
    else//no threshold
    {
        std::cout<<"start parallelized prediction"<<std::endl;

        vigra::FImage probability_image_0(dim_x,dim_y);

        #pragma omp parallel
        {
            float * pixel_unknown_features_array=new float[nr_of_features];
            vigra::MultiArray<2, float> pixel_unknown_features(vigra::MultiArrayShape<2>::type(1, nr_of_features),pixel_unknown_features_array);

            #pragma omp for
            for(int x=0;x<dim_x;x++)
            {
                for (int y=0;y<dim_y;y++)
                {
                    int z=dim_x*y+x;

                    for(int feature=0;feature<nr_of_features;feature++)
                    {
                        //pixel_unknown_features(0,feature)=unknown_features_array[nr_of_features*z+feature];
                        pixel_unknown_features(0,feature)=unknown_features_array[feature][z];
                    }

                    vigra::MultiArray<2, double> pixel_unknown_probability(vigra::MultiArrayShape<2>::type(1,2));

                    rf.predictProbabilities(pixel_unknown_features,pixel_unknown_probability);

                    probability_image_0(x,y)=(float)pixel_unknown_probability(0,0);
                }
            }

            delete pixel_unknown_features_array;  

        }//end of parallelisation

        std::cout<<"end parallelized prediction"<<std::endl;

        /* If there exists a bubble image, fill the bubble areas on the probability image
         * and draw their outlines on the watershed label image as well as the probability image.
         */
        if(bubble_img == true)
        {                    
            for(int x = 0; x < dim_x; x++)
            {
                for(int y = 0; y < dim_y; y++)
                {
                    if(bubbles(x,y) == true)
                    {                
                        probability_image_0(x,y) = 1.0f;                       
                    }
                }
            }
            
            //Traverse along x-axis first
            for(int x = 0; x < dim_x; x++)
            {
                for(int y = 0; y < dim_y; y++)
                {
                    if(bubbles(x,y) == false && x+1 < dim_x)
                    {
                        if(bubbles(x+1,y) == true)
                        {
                            probability_image_0(x+1,y) = 0.0f;                            
                            if(x+2 < dim_x && y+1 < dim_y)
                            {
                                probability_image_0(x+2,y) = 0.0f;
                                probability_image_0(x+1,y+1) = 0.0f;
                            }
                            if(y-1 >= 0)
                            {
                                probability_image_0(x,y) = 0.0f;
                                probability_image_0(x+1,y-1) = 0.0f;
                            }
                        }
                    }
                    if(bubbles(x,y) == true && x+1 < dim_x)
                    {
                        if(bubbles(x+1,y) == false)
                        {
                            probability_image_0(x,y) = 0.0f;                    
                            if(x+1 < dim_x && y+1 < dim_y)
                            {
                                probability_image_0(x+1,y) = 0.0f;
                                probability_image_0(x,y+1) = 0.0f;
                            }
                            if(x-1 >= 0 && y-1 >= 0)
                            {
                                probability_image_0(x-1,y) = 0.0f;
                                probability_image_0(x,y-1) = 0.0f;
                            }
                        }
                    }
                }
            }

            //Traverse along y-axis 
            for(int y = 0; y < dim_y; y++)
            {
                for(int x = 0; x < dim_x; x++)
                {
                    if(bubbles(x,y) == false && y+1 < dim_y)
                    {
                        if(bubbles(x,y+1) == true)
                        {
                            probability_image_0(x,y+1) = 0.0f;                            
                            if(x+1 < dim_x && y+2 < dim_y)
                            {
                                probability_image_0(x+1,y+1) = 0.0f;
                                probability_image_0(x,y+2) = 0.0f;
                            }
                            if(y-1 >= 0)
                            {
                                probability_image_0(x-1,y+1) = 0.0f;
                                probability_image_0(x,y) = 0.0f;
                            }
                        }
                    }
                    if(bubbles(x,y) == true && y+1 < dim_y)
                    {
                        if(bubbles(x,y+1) == false)
                        {
                            probability_image_0(x,y) = 0.0f;                            
                            if(x+1 < dim_x && y+2 < dim_y)
                            {
                                probability_image_0(x+1,y+1) = 0.0f;
                                probability_image_0(x,y+2) = 0.0f;                               
                            }
                            if(y-1 >= 0)
                            {
                                probability_image_0(x-1,y+1) = 0.0f;
                                probability_image_0(x,y) = 0.0f;                                
                            }
                        }
                    }
                }
            }
        }
        exportImage(srcImageRange(probability_image_0), vigra::ImageExportInfo(dest_path_probability_maps.c_str()));

        vigra::BasicImage<bool> temp_img(dim_x,dim_y);

        if(bubble_img) std::cout<<"threshold on probability map, combine with bubble image"<<std::endl;
        else std::cout<<"threshold on probability map"<<std::endl;

     	for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
           		if(probability_image_0(x,y)<0.78f) temp_img(x,y)=true;//TODO
                else temp_img(x,y)=false;
            }
        }

        std::cout<<"find connected filtered pixels"<<std::endl;

        vigra::IImage labels(dim_x,dim_y);
        int nr_labels;

        // find connected regions
        nr_labels=vigra::labelImageWithBackground(srcImageRange(temp_img), destImage(labels), true, 0);

        int * label_size=new int[nr_labels+1];
        int * bubble_pixels=new int[nr_labels+1];

        for(int i=0; i<nr_labels+1; i++)
        {
            label_size[i]=0;
            bubble_pixels[i]=0;
        }

        for(int x=0;x<dim_x;x++)
        {
            for (int y=0;y<dim_y;y++)
            {
                label_size[labels(x,y)]++;
                if (bubbles(x+1,y)) bubble_pixels[labels(x,y)]++;
            }
        }

        for(int x=0;x<dim_x;x++)
        {
            for (int y=0;y<dim_y;y++)
            {
                if (label_size[labels(x,y)]<1000 && bubble_pixels[labels(x,y)]<100) temp_img(x,y)=false;
                //if (label_size[labels(x,y)]<1000) std::cout<<bubble_pixels[labels(x,y)]<<std::endl;
            }
        }

        delete label_size;

        std::cout<<"close filtered pixels"<<std::endl;

        vigra::FImage temp_img_float2(dim_x,dim_y);
        int r=5;

        //Copy from buffer
        for(int i=0; i<temp_img.width(); i++)
            for(int j=0; j<temp_img.height(); j++)
            {
                if (temp_img(i,j)==false) temp_img_float2(i,j) = 0.0f;
                else temp_img_float2(i,j) = 5.0f;
            }

        vigra::FImage temp_img_float1(temp_img_float2);

        // Disc Dilation with radius r
        for(int i=0; i<temp_img.width(); i++)
            for(int j=0; j<temp_img.height(); j++)
                for(int k = -r; k<=r; k++)
                    for(int l= -r; l<=r; l++)
                    {
                        if( k*k + l*l <= r*r )
                        {
                            if( k+i >= 0 && k+i < temp_img_float2.width() && j+l >=0 && j+l < temp_img_float2.height() )
                            {
                                if( temp_img_float2( k+i, j+l ) > temp_img_float1( i, j ) )
                                    temp_img_float1( i, j) = temp_img_float2( k+i, j+l );
                            }
                        }
                    }

        //Copy from buffer
        for(int i=0; i<temp_img.width(); i++)
            for(int j=0; j<temp_img.height(); j++)
                temp_img_float2(i,j) = temp_img_float1(i,j);

        // Disc Erosion with radius r
        for(int i=0; i<temp_img.width(); i++)
            for(int j=0; j<temp_img.height(); j++)
                for(int k = -r; k<=r; k++)
                    for(int l= -r; l<=r; l++)
                    {
                        if( k*k + l*l <= r*r )
                        {
                            if( k+i >= 0 && k+i < temp_img_float2.width() && j+l >=0 && j+l < temp_img_float2.height() )
                            {
                                if( temp_img_float2( k+i, j+l ) < temp_img_float1( i, j ) )
                                    temp_img_float1( i, j) = temp_img_float2( k+i, j+l );
                            }
                        }
                    }

        //Copy from buffer
        for(int i=0; i<temp_img.width(); i++)
            for(int j=0; j<temp_img.height(); j++)
            {
                if (temp_img_float1(i,j)>0.0f) temp_img(i,j) = 0;
                else temp_img(i,j) = 1;
            }

        dest_path_probability_maps.resize(dest_path_probability_maps.size()-3);
        dest_path_probability_maps.append("0.bmp");
        exportImage(srcImageRange(temp_img), vigra::ImageExportInfo(dest_path_probability_maps.c_str()));
    }

    gettimeofday(&end, 0);
    int sec=0;
    int msec=(int)((end.tv_usec-start.tv_usec)/1000);
    if (msec<0)
    {
        msec=msec+1000;
        sec=-1;
    }
    int min=(int)((end.tv_sec-start.tv_sec)/60);
    sec=sec+end.tv_sec-start.tv_sec-60*min;
    if (sec<0)
    {
        min=min-1;
        sec=sec+60;
    }

    std::cout << "Laufzeit: "<<min<<" min "<<sec<<','<<msec<<" s"<< std::endl;

    std::cout<<"....done"<<std::endl;

    //delete unknown_features_array;
    for(int feature=0;feature<nr_of_features;feature++)
    {
        delete [] unknown_features_array[feature];
    }
    delete [] unknown_features_array;
}
