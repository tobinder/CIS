/*! \file junction_features.h
 * \brief Junction feature calculations.
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

#include "CImg.h"

#include "cgp_structure.hxx"

/*! \fn compute_single_junction_features(std::vector< std::vector<point> > arcs,
                                                    std::vector <point> junctions,
                                                    marray::Marray<unsigned int> const & one_boundings,
                                                    int junction_index,
                                                    std::vector<float> & angles)
 * \brief Compute single junction features.
 * \param arcs Vector containing the segmentation arcs
 * \param junctions Vector containing the segmentation junctions
 * \param one_boundings Marray representing one boundings
 * \param junction_index Junction index
 * \param angles Vector containing the angles
 * \return Float vector containing the features.
 */
std::vector<float> compute_single_junction_features(std::vector< std::vector<point> > arcs,
                                                    std::vector <point> junctions,
                                                    marray::Marray<unsigned int> const & one_boundings,
                                                    int junction_index,
                                                    std::vector<float> & angles)
{
    std::vector<float> features;

    //calculate angles in junctions

    //HERE WE HAVE TO LOOK UP THE COORDINATS OF THE 1-Cell
    double x_center=junctions[junction_index].x;
    double y_center=junctions[junction_index].y;

    size_t nr_of_boundarys_crossing;
    std::vector<size_t> index_of_boundings;
    if(one_boundings(junction_index,3) == 0)
    {
        index_of_boundings.resize(3);
        nr_of_boundarys_crossing = 3;
    }
    else
    {
        index_of_boundings.resize(4);
        nr_of_boundarys_crossing = 4;
    }

    // now we look up the boundary indexes and write them into the index_of_boundings vector
    for(size_t x=0; x<nr_of_boundarys_crossing; x++)
    {
        index_of_boundings[x]=one_boundings(junction_index,x);
    }

    angles.resize(nr_of_boundarys_crossing);//unsorted angles

    int min_arc_length=1000;

    for(size_t i=0;i<nr_of_boundarys_crossing;i++)
    {
        if (arcs[index_of_boundings[i]-1].size()<min_arc_length) min_arc_length=arcs[index_of_boundings[i]-1].size();
    }

    size_t nr_of_arcs_to_average=std::min(100,min_arc_length/3);

    if(nr_of_arcs_to_average>=20)
    {    
        //we have to loop over 3 or over 4 boudarys (2-cells)
        for(size_t i=0;i<nr_of_boundarys_crossing;i++)
        {
            //here we loop over the points of the boundary
            //we have to computed the weighted "center"of pixels
            double weights=0;
            double sum_x=0;
            double sum_y=0;
            double av_x;
            double av_y;

            if (fabs(x_center-arcs[index_of_boundings[i]-1][0].x)<2 && fabs(y_center-arcs[index_of_boundings[i]-1][0].y)<2)//first point in arc is close to junction
            {
                for(size_t k=0; k<nr_of_arcs_to_average; k++)
                {
                    double xx=arcs[index_of_boundings[i]-1][k].x;
                    double yy=arcs[index_of_boundings[i]-1][k].y;
                    double weight_factor=1;
                    sum_x=sum_x+(xx*weight_factor);
                    sum_y=sum_y+(yy*weight_factor);
                    weights=weights+weight_factor;
                }

                // normalize
                av_x=sum_x/weights;
                av_y=sum_y/weights;

                size_t k=0;
                while ((av_x-x_center)==0 && k<arcs[index_of_boundings[i]-1].size())
                {
                    av_x=x_center+((arcs[index_of_boundings[i]-1][k].x-x_center)/(k+1));
                    k++;
                }

                k=0;
                while ((av_y-y_center)==0 && k<arcs[index_of_boundings[i]-1].size())
                {
                    av_y=y_center+((arcs[index_of_boundings[i]-1][k].y-y_center)/(k+1));
                    k++;
                }

            }
            else if (fabs(x_center-arcs[index_of_boundings[i]-1][arcs[index_of_boundings[i]-1].size()-1].x)<2 &&
                     fabs(y_center-arcs[index_of_boundings[i]-1][arcs[index_of_boundings[i]-1].size()-1].y)<2)//last point in arc is close to junction
            {
                for(size_t k=arcs[index_of_boundings[i]-1].size()-nr_of_arcs_to_average; k<arcs[index_of_boundings[i]-1].size(); k++)
                {
                    double xx=arcs[index_of_boundings[i]-1][k].x;
                    double yy=arcs[index_of_boundings[i]-1][k].y;
                    double weight_factor=1;
                    sum_x=sum_x+(xx*weight_factor);
                    sum_y=sum_y+(yy*weight_factor);
                    weights=weights+weight_factor;
                }

                // normalize
                av_x=sum_x/weights;
                av_y=sum_y/weights;

                size_t k=0;
                while ((av_x-x_center)==0 && k<arcs[index_of_boundings[i]-1].size())
                {
                    av_x=x_center+((arcs[index_of_boundings[i]-1][arcs[index_of_boundings[i]-1].size()-k].x-x_center)/(k+1));
                    k--;
                }

                k=0;
                while ((av_y-y_center)==0 && k<arcs[index_of_boundings[i]-1].size())
                {
                    av_y=y_center+((arcs[index_of_boundings[i]-1][arcs[index_of_boundings[i]-1].size()-k].y-y_center)/(k+1));
                    k--;
                }

            }
            else
            {
                std::cout<<"Error: Arcs are not sorted correctly, no point close to junction found!"<<std::endl;
            }

            angles[i]=(atan2((av_y-y_center),(av_x-x_center)));
        }

        std::vector<float> angles_s(nr_of_boundarys_crossing);//sorted angles
        for(size_t i=0;i<nr_of_boundarys_crossing;i++) angles_s[i]=angles[i];

        std::sort(angles_s.begin(),angles_s.end());//sort by angles

        if(angles_s.size()==3)
        {
            float temp_0 = angles_s[0];
            angles_s[0] = angles_s[1] - angles_s[0];
            angles_s[1] = angles_s[2] - angles_s[1];
            angles_s[2] = 2.0f * 3.141592654f - angles_s[2] + temp_0;
        }
        if(angles_s.size()==4)
        {
            float temp_0 = angles_s[0];
            angles_s[0]=angles_s[1]-angles_s[0];
            angles_s[1]=angles_s[2]-angles_s[1];
            angles_s[2]=angles_s[3]-angles_s[2];
            angles_s[3]=2.0f * 3.141592654f - angles_s[3] + temp_0;
        }

        features.push_back(get_standard_deviation(angles_s));
        //features.push_back(get_min(angles_s));
        //features.push_back(get_max(angles_s));
        features.push_back(get_max(angles_s)-get_min(angles_s));
    }
    else
    {
        features.push_back(0.0f);
        features.push_back(0.0f);
    }

    return features;
}

/*! \fn  junction_features_and_classification_set(std::string source_path_image,
                                              std::string source_ws_image_path,
                                              std::string filepath_to_junction_training_file,
                                              std::string filepath_to_classification_file)
 * \brief Compute junction features and a classification set from a given image.
 * \param source_path_image Path to the source image
 * \param source_ws_image_path Path to the watershed segmentation image
 * \param filepath_to_junction_training_file Filepath to the junction training file
 * \param filepath_to_classification_file Filepath to the classification file
 */
void junction_features_and_classification_set(std::string source_path_image,
                                              std::string source_ws_image_path,
                                              std::string filepath_to_junction_training_file,
                                              std::string filepath_to_classification_file)
{
    //IMPORT RESULTS FROM HDF5 file
    //data structur
    //arcs == two cells in cartesian coordinats
    std::vector< std::vector<point> > arcs;
    std::vector<point> junctions;

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

    //determine the number of features
    std::vector<float> junction_feature;
    std::vector<float> angles;
    junction_feature=compute_single_junction_features(arcs,junctions,one_boundings,0,angles);
    angles.clear();

    float nr_of_features=(float)junction_feature.size();

    int nr_of_junctions=junctions.size();

    //Now we can create an pointer to an array to store all the features
    float * junction_features_array=new float[(int)(nr_of_junctions*nr_of_features)];

    std::cout<<"compute junction-features..."<<std::endl;
    std::cout<<"nr of junctions:"<<nr_of_junctions<<" nr of junction features:"<<nr_of_features<<std::endl;

    for(int i=0;i<(int)junctions.size();i++)
    {
        junction_feature=compute_single_junction_features(arcs,junctions,one_boundings,i,angles);
        angles.clear();//not needed for further calculations

        for(int jj=0;jj<(int)junction_feature.size();jj++)
        {
            junction_features_array[(int)(i*nr_of_features+jj)]=junction_feature[jj];
        }
    }
    std::cout<<"...done"<<std::endl;

    //here is some dirty stuff,the last line is readed twice(eof cant look in future),so we define a variable, what the last wirte was
    //and we check if the last read hast already been written, and i use -1 as start value because there cant be a regular junction nr -1
    int last_read_junctions=-1;

    std::string filepath_to_classification_info=filepath_to_classification_file;
    filepath_to_classification_info.append(".info.bin");
    filepath_to_classification_file.append(".bin");

    //Now we know how many junctions are in the image, so we can create and std::vector vector_training_data
    //the index of the vector is the index of the junction, the value to the index is the label of junction
    //if the label=10 this means that there is no training data to this junction

    std::vector<int> vector_training_labels;
    vector_training_labels.resize((int)nr_of_junctions,10);

    std::ifstream training_file(filepath_to_junction_training_file.c_str());
    if(!training_file)
    {
        std::cout<<"no training data found in: "<<filepath_to_junction_training_file<<std::endl;
        exit(-1);
    }
    else  //FILE IS EXISTEND
    {
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

            if(temp[0]>nr_of_junctions-1)
            {
                std::cout<<"error in create_junction_classification_set(..)"<<std::endl;
                std::cout<<"trainingsdata and image have different nr of junctions:"<<std::endl;
                std::cout<<"found labeling for junction nr "<<temp[0]<<" in trainingsdata."<<std::endl;
                exit(-1);
            }
            else if(last_read_junctions!=temp[0])
            {
                if(temp[1]>2-1 || temp[1]<0)
                {
                    std::cout<<"error in create_junction_classification_set(..)"<<std::endl;
                    std::cout<<"label "<<temp[1]<<" is not valid. Nr of classes: "<<2<<std::endl;
                    exit(-1);
                }

                vector_training_labels[temp[0]]=temp[1];
                //HERE WE SET LAST READ TO THE LAST READED /WRITTEN JUNCTIONS TO AVOID THE LAST LINE
                //IS READED TWICE
                last_read_junctions=temp[0];
                i++;
            }
        }

        training_file.close();
        std::cout<<"training data has been loaded"<<std::endl;

        //creates a new classifaction set file
        FILE * classification_set_file;
        classification_set_file=fopen(filepath_to_classification_file.c_str(),"w+b");
        classification_set_file=fopen(filepath_to_classification_file.c_str(),"w");

        std::vector<int> class_count;
        class_count.resize(10,0);

        for(int i=0;i<(int)vector_training_labels.size();i++)
        {
            //if vector_training_labels[i] ==10 this means that we dont have training
            //data for this junction

            if(vector_training_labels[i]!=10)
            {
                //vector to store the label and the feature values
                std::vector<float> label_and_features;

                //first entry in labels_and_features should be the label
                label_and_features.push_back((float)vector_training_labels[i]);

                //we want so know how many data is from wich class
                class_count[vector_training_labels[i]]=class_count[vector_training_labels[i]]+1;

                std::cout<<"junction "<<i<<" with label "<<vector_training_labels[i]<<" has features: ";

                //j is the index the junction feature array
                for(int j=(i)*(int)nr_of_features;j<(i+1)*(int)nr_of_features;j++ )
                {
                    float feature_value=junction_features_array[j];
                    std::cout<<feature_value<<" ";
                    label_and_features.push_back(feature_value);
                }

                std::cout<<std::endl;

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

        //Now we must look know how many data has been written to the file
        int nr_of_new_data=0;

        if(class_count[0]!=class_count[1])
        {
            std::cout<<"unbalanced data in "<<get_filename(filepath_to_junction_training_file)<<std::endl;
            //exit(-1);
        }

        std::cout<<"Nr of training data class 0: "<<class_count[0]<<std::endl;
        std::cout<<"Nr of training data class 1: "<<class_count[1]<<std::endl;

        float all_data=(float)(class_count[0]+class_count[1]);

        //creates a new classifaction set info file
        FILE * info;
        info=fopen(filepath_to_classification_info.c_str(),"w+b");

        fwrite(&all_data,sizeof(float),1,info);
        fwrite(&nr_of_features,sizeof(float),1,info);

        fclose(info);

        //and we can delete the junction feature array
        delete junction_features_array;
    }
}
