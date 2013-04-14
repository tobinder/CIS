/*! \file pixel_complete.h
 * \brief All pixel functions
 *
 * This header file contains functions that perform all pixel functions.
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

/*! \fn do_pixel_complete(ParameterFile paramFile,std::string filepath_image,std::string param_file,bool *default_rf_trained, int size, int rank)
 * \brief All pixel functions.
 * \param paramFile Parameter file
 * \param filepath_image Filepath to the image
 * \param param_file Name of the parameter file
 * \param default_rf_trained Has the default random forest been trained?
 * \param size Size
 * \param rank Rank
 */
void do_pixel_complete(ParameterFile paramFile,std::string filepath_image,std::string param_file,bool *default_rf_trained, int size, int rank)
{
    Parameter<int> threshold;
   	threshold.assign("", "threshold", 20);
    threshold.load(paramFile,"config");

   	Parameter<double> scale;
   	scale.assign("", "scale", 5);
    scale.load(paramFile,"config");

    Parameter<int> equal_tolerance;
   	equal_tolerance.assign("", "equal_tolerance", 1);
    equal_tolerance.load(paramFile,"config");

   	Parameter<std::string> fp_as_image;
   	fp_as_image.assign("", "filepath_as_image", "as/");
    fp_as_image.load(paramFile,"config");

   	Parameter<std::string> path_pixel_label;
   	path_pixel_label.assign("", "path_pixel_label", "pixel-labels/");
    path_pixel_label.load(paramFile,"config");

   	Parameter<std::string> p_pixel_feature;
   	p_pixel_feature.assign("", "path_pixel_feature", "pixel-features/");
    p_pixel_feature.load(paramFile,"config");
    std::string path_pixel_feature=p_pixel_feature;

   	Parameter<std::string> path_pixel_classification;
   	path_pixel_classification.assign("", "path_pixel_classification", "pixel-classification/");
    path_pixel_classification.load(paramFile,"config");

   	Parameter<std::string> path_pixel_rf;
   	path_pixel_rf.assign("", "path_pixel_rf", "pixel-classification/random-forests/");
    path_pixel_rf.load(paramFile,"config");

   	Parameter<std::string> path_prob_map;
   	path_prob_map.assign("", "path_prob_map", "pixel-classification/probability-maps/");
    path_prob_map.load(paramFile,"config");

   	Parameter<std::string> p_watershed;
   	p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
    p_watershed.load(paramFile,"config");
    std::string path_watershed=p_watershed;

   	Parameter<std::string> p_boundary_feature;
   	p_boundary_feature.assign("", "path_boundary_feature", "boundary-features/");
    p_boundary_feature.load(paramFile,"config");
    std::string path_boundary_feature=p_boundary_feature;

   	Parameter<std::string> default_pixel_label;
   	default_pixel_label.assign("", "default_pixel_label", "");
    default_pixel_label.load(paramFile,"config");

   	Parameter<std::string> fp_default_image;
   	fp_default_image.assign("", "filepath_default_image", "not_found");
    fp_default_image.load(paramFile,"config");
    std::string filepath_default_image = fp_default_image;

   	Parameter<std::string> filepath_default_as_image;
   	filepath_default_as_image.assign("", "filepath_default_as_image", "");
    filepath_default_as_image.load(paramFile,"config");

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "") param_file_name.resize(param_file_name.size()-4);
    else param_file = "parameters.txt";

    std::string filepath_default_pixel_classification=path_pixel_classification;
    filepath_default_pixel_classification.append(default_pixel_label);
    std::string filepath_default_pixel_rf=path_pixel_rf;
    filepath_default_pixel_rf.append(default_pixel_label);

    std::string filename_of_image=get_filename(filepath_image);
    std::string filepath_as_image=fp_as_image;
    filepath_as_image.append(filename_of_image);
    filepath_as_image.append(".as.bmp");
    
    std::string filepath_pixel_label=path_pixel_label;
    filepath_pixel_label.append(filename_of_image);
    filepath_pixel_label.append(".dat");

    path_pixel_feature.append(param_file_name);
    path_boundary_feature.append(param_file_name);
    
    std::string filepath_pixel_classification=path_pixel_classification;
    filepath_pixel_classification.append(filename_of_image);
    filepath_pixel_classification.resize(filepath_pixel_classification.size()-4);
    filepath_pixel_classification.append(param_file_name);
    
    std::string filepath_pixel_rf=path_pixel_rf;
    filepath_pixel_rf.append(filename_of_image);
    filepath_pixel_rf.resize(filepath_pixel_rf.size()-4);
    filepath_pixel_rf.append(param_file_name);
    
    std::string filepath_prob_map=path_prob_map;
    filepath_prob_map.append(param_file_name);
    filepath_prob_map.append(filename_of_image);

    std::string filepath_boundary_feature=path_boundary_feature;
    filepath_boundary_feature.append(filename_of_image);
    filepath_boundary_feature.append(".bin");

//for (int i=0;i<10;i++)
//{
    //Check if current image has default labels and features already are computed
    if ((filepath_image==filepath_default_image) && (*default_rf_trained==true))
    {
        std::cout<<"Features of this image already are computed"<<std::endl;
        extract_pixel_probabilities(filepath_image.c_str(),path_pixel_feature,filepath_default_pixel_rf,path_prob_map,param_file);
    }
    else
    {
        //Check if pixel labels for this image exist
        FILE *fp_test;
        fp_test=fopen(filepath_pixel_label.c_str(),"r");

        //Check if pixel labels for this image are default labels
        if (filepath_image==filepath_default_image) *default_rf_trained=true;
        
        //Check if grayscale or color
        vigra::ImageImportInfo info(filepath_image.c_str());
        
        //Check if classification set for this image exist
        FILE *fp_set;
        std::string test = filepath_pixel_classification;
        test.append(".info.bin");
        fp_set=fopen(test.c_str(),"r");

        if(fp_test!=NULL) //labels for this picture
        {
            if(fp_set!=NULL) //classification set for this picture
            {
                std::cout<<"Classification set found"<<std::endl;            
                fclose(fp_set);
            } 
            else
            {
                std::cout<<"Pixel labels found"<<std::endl;
                if(info.isGrayscale())
                {
                    FeatureExtractorGray pixel_features(filepath_image.c_str(),filepath_as_image,path_pixel_feature);
                    //compute features on the whole image and save to file
                    int nr_of_features = pixel_features.extract(true,param_file);
                    //write features of the training data to classification-set
                    pixel_features.extract_training_data(path_pixel_label,filepath_pixel_classification,nr_of_features);
                }
                else
                {
                    FeatureExtractor pixel_features(filepath_image.c_str(),filepath_as_image,path_pixel_feature);
                    //compute features on the whole image and save to file
                    int nr_of_features = pixel_features.extract(true,param_file);
                    //write features of the training data to classification-set
                    pixel_features.extract_training_data(path_pixel_label,filepath_pixel_classification,nr_of_features);
                }
            }

            std::cout<<"Train random forest for this image"<<std::endl;
            //train rf with pixel-classification-set    
            train_pixel_rf(filepath_pixel_classification,filepath_pixel_rf,param_file);
        }
        else //no labels for this picture
        {
            //Check if default classification set exist
            FILE *fp_default_set=NULL;
            std::string test = filepath_default_pixel_classification;
            test.append(".info.bin");
            //fp_default_set=fopen(test.c_str(),"r");

            std::cout<<"No Pixel labels found"<<std::endl;
            if(filepath_default_image == "not_found")
            {
                std::cout<<"Error: No default pixel labels found!"<<std::endl;
                exit(-1);
            }
            if(filepath_default_image == filepath_image)
            {
                std::cout<<"Error: Check label path!"<<std::endl;
                exit(-1);
            }
            if(info.isGrayscale())
            {
                //Check if pixel features have been calculated
                FILE *fp_feature_test;
                std::string filepath_pixel_features=path_pixel_feature;
                filepath_pixel_features.append(filename_of_image);
                filepath_pixel_features.append(".bin");
                fp_feature_test=fopen(filepath_pixel_features.c_str(),"r");

                if(fp_feature_test!=NULL) //features exist
                {
                    std::cout<<"Features of this image already are computed"<<std::endl;
                    fclose(fp_feature_test);
                }
                else
                {
                    FeatureExtractorGray pixel_features(filepath_image.c_str(),filepath_as_image,path_pixel_feature);
                    //compute features on the whole image and save to file
                    pixel_features.extract(true,param_file);
                }
        
                if (*default_rf_trained==false && fp_default_set==NULL)
                {   //compute features on the image corresponding to the default labels
                    FeatureExtractorGray pixel_features_default(filepath_default_image.c_str(),filepath_default_as_image,path_pixel_feature);
                    int nr_of_features = pixel_features_default.extract(true,param_file);
                    //create default classification-set
                    std::cout<<"Compute default pixel classification-set"<<std::endl;
                    pixel_features_default.extract_training_data(path_pixel_label,filepath_default_pixel_classification,nr_of_features);
                }
            
            }
            else//color
            {
                //Check if pixel features have been calculated
                FILE *fp_feature_test;
                std::string filepath_pixel_features=path_pixel_feature;
                filepath_pixel_features.append(filename_of_image);
                filepath_pixel_features.append(".bin");
                fp_feature_test=fopen(filepath_pixel_features.c_str(),"r");

                if(fp_feature_test!=NULL) //features exist
                {
                    std::cout<<"Features of this image already are computed"<<std::endl;
                    fclose(fp_feature_test);
                }
                else
                {
                    FeatureExtractor pixel_features(filepath_image.c_str(),filepath_as_image,path_pixel_feature);
                    //compute features on the whole image and save to file
                    pixel_features.extract(true,param_file);
                }
        
                if (*default_rf_trained==false && fp_default_set==NULL)
                {   //compute features on the image corresponding to the default labels
                    FeatureExtractor pixel_features_default(filepath_default_image.c_str(),filepath_default_as_image,path_pixel_feature);
                    int nr_of_features = pixel_features_default.extract(true,param_file);
                    //create default classification-set
                    std::cout<<"Compute default pixel classification-set"<<std::endl;
                    pixel_features_default.extract_training_data(path_pixel_label,filepath_default_pixel_classification,nr_of_features);
                }
            }
        
            if (*default_rf_trained==false)
            {
                if(fp_default_set!=NULL) //classification set for this picture
                {
                    std::cout<<"default classification-set found"<<std::endl;            
                    fclose(fp_default_set);
                } 

                std::cout<<"Train default random forest"<<std::endl;
                train_pixel_rf(filepath_default_pixel_classification,filepath_default_pixel_rf,param_file);
                *default_rf_trained=true;
            } else std::cout<<"Random forest already trained"<<std::endl;
        
        }//end of else/no labels for this picture
        
        //rf-prediction for the whole image, create a probability-map
        if(fp_test!=NULL)
        {
            extract_pixel_probabilities(filepath_image.c_str(),path_pixel_feature,filepath_pixel_rf,path_prob_map,param_file);
            fclose(fp_test);
        }
        else extract_pixel_probabilities(filepath_image.c_str(),path_pixel_feature,filepath_default_pixel_rf,path_prob_map,param_file);
    } 
    //compute watershed-regions
    compute_watershed_regions(filepath_prob_map,path_watershed,threshold,scale,equal_tolerance,size,rank);

    //compute boundary features
    boundary_features_gray(filepath_image.c_str(),path_watershed,path_prob_map,path_pixel_feature,path_boundary_feature,param_file);
    //}
}

/*! \fn do_pixel_complete(ParameterFile paramFile,
                       std::string filepath_image,
                       std::string param_file,
                       std::string filepath_as_image,
                       std::string path_pixel_feature,
                       std::string path_prob_map,
                       std::string path_watershed,
                       std::string path_boundary_feature,
                       int size,
                       int rank)
 * \brief All pixel functions 
 * \param paramFile Parameter file
 * \param filepath_image Filepath to the image
 * \param param_file Name of the parameter file
 * \param filepath_as_image Filepath to the preprocessed image
 * \param path_pixel_feature Path to the pixel features
 * \param path_prob_map Path to the probability map
 * \param path_watershed Path to the watershed segmentation
 * \param path_boundary_feature Path to the boundary features
 * \param size Size
 * \param rank Rank
 */
void do_pixel_complete(ParameterFile paramFile,
                       std::string filepath_image,
                       std::string param_file,
                       std::string filepath_as_image,
                       std::string path_pixel_feature,
                       std::string path_prob_map,
                       std::string path_watershed,
                       std::string path_boundary_feature,
                       int size,
                       int rank)
{
    Parameter<int> threshold;
   	threshold.assign("", "threshold", 20);
    threshold.load(paramFile,"config");

   	Parameter<double> scale;
   	scale.assign("", "scale", 5);
    scale.load(paramFile,"config");

    Parameter<int> equal_tolerance;
   	equal_tolerance.assign("", "equal_tolerance", 1);
    equal_tolerance.load(paramFile,"config");

    Parameter<std::string> filepath_default_pixel_rf;

    /* In order to avoid changing the signature of 'extract_boundary_probabilities', the path to the folders
     * containing the random forests and boundary features is stored in the third to last argument if the option
     * 'Use default Random Forest' is unchecked in the GUI, in which case 'FALSE' is appended to the end of the argument string.
     * The strings containing the two paths are seperated by a space character.
     */
    std::string arg_path_rf = path_pixel_feature;
    std::string filepath_boundary_features;
    if(arg_path_rf.compare(arg_path_rf.size()-5, 5, "FALSE") == 0)
    {            
        arg_path_rf.resize(arg_path_rf.size()-5);            
        size_t found;
        std::string space = " ";
        found = arg_path_rf.find(space);
        if(found != std::string::npos)
        {
            path_pixel_feature = arg_path_rf;
            path_pixel_feature = path_pixel_feature.substr(found+1, arg_path_rf.size()-found);
            std::cout << path_pixel_feature << std::endl;

            filepath_default_pixel_rf = arg_path_rf.substr(0, found);
        }                
    }
    else if(arg_path_rf.compare(arg_path_rf.size()-4, 4, "TRUE") == 0)
    {            
        arg_path_rf.resize(arg_path_rf.size()-4);
        filepath_default_pixel_rf.assign("", "filepath_default_pixel_rf", "NEEMcomplete");
        filepath_default_pixel_rf.load(paramFile,"config");        
        
        path_pixel_feature = arg_path_rf;
    }
    else std::cout << "This was not supposed to happen" << std::endl;

    std::cout << "Using random forest file " << filepath_default_pixel_rf() << std::endl;
        
    Parameter<std::string> path_thumbs;
    path_thumbs.assign("", "path_thumbs", "no");
    path_thumbs.load(paramFile,"config");

    std::string filename_of_image=get_filename(filepath_image);
    filepath_as_image.append(filename_of_image);
    filepath_as_image.append(".as.bmp");
    
    std::string filepath_prob_map=path_prob_map;
    filepath_prob_map.append(filename_of_image);

    std::string filepath_boundary_feature=path_boundary_feature;
    filepath_boundary_feature.append(filename_of_image);
    filepath_boundary_feature.append(".bin");

    //Check if pixel features have been calculated
    FILE *fp_feature_test;
    std::string filepath_pixel_features=path_pixel_feature;
    filepath_pixel_features.append(filename_of_image);
    filepath_pixel_features.append(".bin");
    fp_feature_test=fopen(filepath_pixel_features.c_str(),"r");

    if(fp_feature_test!=NULL) //features exist
    {
        std::cout<<"Features of this image are already computed"<<std::endl;
        fclose(fp_feature_test);
    }
    else
    {
        //Check if grayscale or color
        vigra::ImageImportInfo info(filepath_image.c_str());

        if(info.isGrayscale())
        {
            FeatureExtractorGray pixel_features(filepath_image.c_str(),filepath_as_image,path_pixel_feature);
            //compute features on the whole image and save to file
            pixel_features.extract(true,param_file);
        }
        else//color
        {
            FeatureExtractor pixel_features(filepath_image.c_str(),filepath_as_image,path_pixel_feature);
            //compute features on the whole image and save to file
            pixel_features.extract(true,param_file);
        }
    }

    //rf-prediction for the whole image, create a probability-map
    extract_pixel_probabilities(filepath_image.c_str(),path_pixel_feature,filepath_default_pixel_rf(),path_prob_map,param_file);
 
    //compute watershed-regions
    compute_watershed_regions(filepath_prob_map,path_watershed,threshold,scale,equal_tolerance,size,rank);

    //compute boundary features
    boundary_features_gray(filepath_image.c_str(),path_watershed,path_prob_map,path_pixel_feature,path_boundary_feature,param_file,path_thumbs);
}
