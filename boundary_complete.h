/*! \file boundary_complete.h
 * \brief All boundary functions.
 *
 * This header file contains all boundary functions in one function.
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

/*! \fn do_boundary_complete(ParameterFile paramFile,std::string filepath_image,std::string param_file,bool *default_rf_trained)
 * \brief Perform all boundary functions.
 * \param paramFile Parameter file
 * \param filepath_image Filepath to the image
 * \param param_file Filepath to the parameter file
 * \param default_rf_trained Has the default random forest been trained?
 */
void do_boundary_complete(ParameterFile paramFile,std::string filepath_image,std::string param_file,bool *default_rf_trained)
{
   	Parameter<std::string> p_watershed;
   	p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
    p_watershed.load(paramFile,"config");
    std::string path_watershed=p_watershed;

   	Parameter<std::string> path_prob_map;
   	path_prob_map.assign("", "path_prob_map", "pixel-classification/probability-maps/");
    path_prob_map.load(paramFile,"config");

   	Parameter<std::string> p_pixel_feature;
   	p_pixel_feature.assign("", "path_pixel_feature", "pixel-features/");
    p_pixel_feature.load(paramFile,"config");
    std::string path_pixel_feature=p_pixel_feature;

   	Parameter<std::string> p_boundary_feature;
   	p_boundary_feature.assign("", "path_boundary_feature", "boundary-features/");
    p_boundary_feature.load(paramFile,"config");
    std::string path_boundary_feature=p_boundary_feature;

   	Parameter<std::string> path_boundary_label;
   	path_boundary_label.assign("", "path_boundary_label", "boundary-labels/");
    path_boundary_label.load(paramFile,"config");

   	Parameter<std::string> path_default_prob_map;
   	path_default_prob_map.assign("", "path_default_prob_map", "pixel-classification/probability-maps/");
    path_default_prob_map.load(paramFile,"config");

   	Parameter<std::string> path_default_watershed;
   	path_default_watershed.assign("", "path_default_watershed", "pixel-classification/watershed-segmentation/");
    path_default_watershed.load(paramFile,"config");

   	Parameter<std::string> default_boundary_label;
   	default_boundary_label.assign("", "default_boundary_label", "");
    default_boundary_label.load(paramFile,"config");

   	Parameter<std::string> fp_default_image;
   	fp_default_image.assign("", "filepath_default_image", "not_found");
    fp_default_image.load(paramFile,"config");
    std::string filepath_default_image = fp_default_image;

   	Parameter<std::string> path_boundary_classification;
   	path_boundary_classification.assign("", "path_boundary_classification", "boundary-classification/trainingsets/");
    path_boundary_classification.load(paramFile,"config");

   	Parameter<std::string> path_boundary_rf;
   	path_boundary_rf.assign("", "path_boundary_rf", "boundary-classification/random-forests/");
    path_boundary_rf.load(paramFile,"config");

   	Parameter<std::string> path_rf_predictions;
   	path_rf_predictions.assign("", "path_rf_predictions", "rf-predictions/");
    path_rf_predictions.load(paramFile,"config");

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "") param_file_name.resize(param_file_name.size()-4);
    else param_file = "parameters.txt";

    //path_pixel_feature.append(param_file_name);
    //path_watershed.append(param_file_name);
    path_boundary_feature.append(param_file_name);

    std::string filepath_default_boundary_label=path_boundary_label;
    filepath_default_boundary_label.append(default_boundary_label);
    filepath_default_boundary_label.append(".bmp.dat");

    std::string filepath_default_boundary_classification=path_boundary_classification;
    filepath_default_boundary_classification.append(default_boundary_label);
    filepath_default_boundary_classification.append(param_file_name);

    std::string filepath_default_boundary_rf=path_boundary_rf;
    filepath_default_boundary_rf.append(default_boundary_label);
    filepath_default_boundary_rf.append(param_file_name);

    std::string filename_of_image=get_filename(filepath_image);
    
    std::string filepath_boundary_feature=path_boundary_feature;
    filepath_boundary_feature.append(filename_of_image);
    filepath_boundary_feature.append(".bin");

    std::string filepath_boundary_label=path_boundary_label;
    filepath_boundary_label.append(filename_of_image);
    filepath_boundary_label.resize(filepath_boundary_label.size()-3);
    filepath_boundary_label.append("bmp.dat");
    
    std::string filepath_boundary_classification=path_boundary_classification;
    filepath_boundary_classification.append(filename_of_image);
    filepath_boundary_classification.resize(filepath_boundary_classification.size()-4);
    filepath_boundary_classification.append(param_file_name);
    
    std::string filepath_boundary_rf=path_boundary_rf;
    filepath_boundary_rf.append(filename_of_image);
    filepath_boundary_rf.resize(filepath_boundary_rf.size()-4);
    filepath_boundary_rf.append(param_file_name);

//for (int i=0;i<2;i++)
for (int i=0;i<1;i++)
{
    //Check if current image has default labels and features already are computed
    if ((filepath_image==filepath_default_image) && (*default_rf_trained==true))
    {
        std::cout<<"Features of this image already are computed"<<std::endl;
        extract_boundary_probabilities(filepath_boundary_feature,path_watershed,filepath_default_boundary_rf,path_rf_predictions,path_rf_predictions,
                                       param_file);
    }
    else
    {
        //Check if boundary labels for this image exist
        FILE *fp_test;
        fp_test=fopen(filepath_boundary_label.c_str(),"r");

        //Check if boundary labels for this image are default labels
        if (filepath_image==filepath_default_image) *default_rf_trained=true;
        
        //Check if grayscale or color
        vigra::ImageImportInfo info(filepath_image.c_str());
        
        //Check if classification set for this image exist
        FILE *fp_set;
        std::string test = filepath_boundary_classification;
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
                std::cout<<"Boundary labels found"<<std::endl;
                if(info.isGrayscale())
                {
                    //compute features on the whole image and save to file
                    boundary_features_gray(filepath_image.c_str(),path_watershed,path_prob_map,path_pixel_feature,path_boundary_feature,param_file);
                    //write features of the training data to classification-set
                    create_boundary_classification_set(filepath_boundary_label,path_boundary_feature,filepath_boundary_classification);
                }
                else
                {
                    //compute features on the whole image and save to file
                    boundary_features(filepath_image.c_str(),path_watershed,path_prob_map,path_pixel_feature,path_boundary_feature);
                    //write features of the training data to classification-set
                    create_boundary_classification_set(filepath_boundary_label,path_boundary_feature,filepath_boundary_classification);
                }
            }
        
            std::cout<<"Train random forest for this image"<<std::endl;
            //train rf with boundary-classification-set    
            train_boundary_rf(filepath_boundary_classification,filepath_boundary_rf,param_file);
        }
        else //no labels for this picture
        {
            if (fp_set!=NULL) fclose(fp_set);
            //Check if default classification set exist
            FILE *fp_default_set=NULL;
            std::string test = filepath_default_boundary_classification;
            test.append(".info.bin");
            fp_default_set=fopen(test.c_str(),"r");

            std::cout<<"No Boundary labels found"<<std::endl;
            if(filepath_default_image == "not_found")
            {
                std::cout<<"Error: No default boundary labels found!"<<std::endl;
                exit(-1);
            }
            if(filepath_default_image == filepath_image)
            {
                std::cout<<"Error: Check label path!"<<std::endl;
                exit(-1);
            }
            if(info.isGrayscale())
            {
                //compute features on the whole image and save to file
                if (i==0) boundary_features_gray(filepath_image.c_str(),path_watershed,path_prob_map,path_pixel_feature,path_boundary_feature,param_file);
        
                if (*default_rf_trained==false && fp_default_set==NULL)
                {   //compute features on the image corresponding to the default labels
                    boundary_features_gray(filepath_default_image.c_str(),path_default_watershed,path_default_prob_map,path_pixel_feature,
                        path_boundary_feature,param_file);
                }
            
            }
            else//color
            {
                //compute features on the whole image and save to file
                boundary_features(filepath_image.c_str(),path_watershed,path_prob_map,path_pixel_feature,path_boundary_feature);
        
                if (*default_rf_trained==false && fp_default_set==NULL)
                {   //compute features on the image corresponding to the default labels
                    boundary_features(filepath_default_image.c_str(),path_default_watershed,path_default_prob_map,path_pixel_feature,path_boundary_feature);
                }
            }
        
            if (*default_rf_trained==false)
            {
                if(fp_default_set!=NULL) //classification set for this picture
                {
                    std::cout<<"default classification-set found"<<std::endl;            
                    fclose(fp_default_set);
                } 
                else
                {
                    //create default classification-set
                    std::cout<<"Compute default boundary classification-set"<<std::endl;
                    create_boundary_classification_set(filepath_default_boundary_label,path_boundary_feature,filepath_default_boundary_classification);
                }

                std::cout<<"Train default random forest"<<std::endl;
                train_boundary_rf(filepath_default_boundary_classification,filepath_default_boundary_rf,param_file);
                *default_rf_trained=true;
            } else std::cout<<"Random forest already computed"<<std::endl;
        
        }//end of else/no labels for this picture
        
        //rf-prediction for the whole image
        if(fp_test!=NULL)
        {
            extract_boundary_probabilities(filepath_boundary_feature,path_watershed,filepath_boundary_rf,path_rf_predictions,path_rf_predictions,param_file);
            fclose(fp_test);
        }
        else extract_boundary_probabilities(filepath_boundary_feature,path_watershed,filepath_default_boundary_rf,path_rf_predictions,path_rf_predictions,
                                            param_file);
    } 
    *default_rf_trained=false;
    }
}
