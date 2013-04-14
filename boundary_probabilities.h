/*!  \file boundary_probabilities.h
 * \brief Boundary probability calculations.
 *
 * This header file contains functions for extracting and saving boundary probabilities.
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

#include <vector>
#include <cstring>
#include <fstream>

#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/impex.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/random_forest.hxx>

#include <vigra/hdf5impex.hxx>
#include <vigra/random_forest_hdf5_impex.hxx>

#include "cgp_structure.hxx"
#include "marray.hxx"
#include <math.h>

/*! \fn train_boundary_rf(std::string filepath_to_classification_set_file,std::string dest_filepath_to_random_forest_file,std::string param_file_name)
 * \brief Train a boundary random forest.
 * \param filepath_to_classification_set_file Filepath to the classification set file
 * \param dest_filepath_to_random_forest_file Filepath to the target random forest file
 * \param param_file_name Filename of the parameter file
 */
void train_boundary_rf(std::string filepath_to_classification_set_file,std::string dest_filepath_to_random_forest_file,std::string param_file_name)
{
    ParameterFile paramFile;

    if( !paramFile.load(param_file_name) )
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }

	Parameter<int> nr_of_trees;
	nr_of_trees.assign("", "nr_of_trees", 256);
    nr_of_trees.load(paramFile,"config");

    int nr_of_boundary_features=0;

    //some filepointers to get the training data
    FILE * info;
    std::string filepath_to_classification_info=filepath_to_classification_set_file;
    filepath_to_classification_info.append(".info.bin");
    info=fopen(filepath_to_classification_info.c_str(),"rb");

    FILE *fp;
    filepath_to_classification_set_file.append(".bin");
    fp=fopen(filepath_to_classification_set_file.c_str(),"rb");

    if(info!=NULL)
    {
        //read in the the nr of training data (should be the first float in the file)
        float float_nr_of_training_data=0;
        int         nr_of_training_data=0;
        int read_in_size=fread(&float_nr_of_training_data,sizeof(float),1,info);
        if(read_in_size!=1)
        {
            std::cout<<"train_boundary_rf(..) was not able to read in the nr of training data"<<std::endl;
            exit(-1);
        }
        else  //classification set has been found,AND we know the nr of training data (as a float)
        {
            //read in the nr of boundary features (should be the second float in the file)            
            float float_nr_of_boundary_features=0;    
            int read_in_boundary_f=fread(&float_nr_of_boundary_features,sizeof(float),1,info);

            if(read_in_boundary_f!=1)
            {
                std::cout<<"train_boundary_rf(..) was not able to read in the nr of boundary features"<<std::endl;
                exit(-1);
            }

            nr_of_training_data=(int)float_nr_of_training_data;
            nr_of_boundary_features=(int)float_nr_of_boundary_features;

            std::cout<<"Nr of training data "<<nr_of_training_data<<std::endl;
            std::cout<<"Nr of boundary features "<<nr_of_boundary_features<<std::endl;

            //now we can create an array to store all the data of the training file
            //except the first two float (nr of training data, nr of boundary features),that has already been readed in
            float * data_cset_array=new float[((nr_of_boundary_features+1)*nr_of_training_data)];
            float * feature_array=new float[((nr_of_boundary_features)*nr_of_training_data)];
            int * label_array=new int[(nr_of_training_data)];

            std::cout<<"read in the labels and features..."<<std::endl;

            int read_in_cset=fread(data_cset_array,sizeof(float),((nr_of_boundary_features+1)*nr_of_training_data),fp);
            if(read_in_cset!=((nr_of_boundary_features+1)*nr_of_training_data))
            {
                std::cout<<"train_boundary_rf(..) was not able to read in all the floats of "<<filepath_to_classification_set_file<<std::endl;
                exit(-1);
            }

            std::cout<<"...done"<<std::endl;

            //we can close the filepointer
            fclose(fp);
            fclose(info);

            //now me must fill the feature_array and the label array
            int ii=0;
            int label_index=0;
            int feature_index=0;
            while(ii<((nr_of_boundary_features+1)*nr_of_training_data))
            {
               //first entry is always the class-label
               if(data_cset_array[ii]>(nr_of_classes-1) || data_cset_array[ii]<0)
               {
                   std::cout<<"error in train_boundary_rf(..)"<<std::endl;
                   std::cout<<"label "<<data_cset_array[ii]<<" is not valid. Nr of classes: "<<nr_of_classes<<std::endl;
                   exit(-1);
               }
               label_array[label_index]=(int)data_cset_array[ii];
               label_index++;
               ii++;

               //after the label there are nr_of_boundary features floats
               int jj=0;
               while(jj<nr_of_boundary_features)
               {
                   feature_array[feature_index]=data_cset_array[ii];
                   feature_index++;
                   jj++;
                   ii++;
               }
            }

            //now we can delete the data_cset_array
            delete data_cset_array;

            //now we have to store the data in a vigra::multiarray!
            vigra::MultiArray<2, float> features(vigra::MultiArrayShape<2>::type(nr_of_training_data, nr_of_boundary_features),feature_array);
            for(int x=0;x<nr_of_training_data;x++)
            {
                for(int y=0;y<nr_of_boundary_features;y++)
                {
                features(x,y)=feature_array[nr_of_boundary_features*x+y];
                }
            }

            vigra::MultiArray<2, int> labels(vigra::MultiArrayShape<2>::type(nr_of_training_data, 1),label_array);


            vigra::RandomForest<> rf(vigra::RandomForestOptions().tree_count(nr_of_trees));
            std::cout<<"Train Random Forest... "<<std::endl;
 
            //vigra::rf::visitors::VariableImportanceVisitor var_imp;

            rf.learn(features,labels);
            //rf.learn(features,labels, create_visitor(var_imp));
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

            /*
            //EXPORT an other "info" file,in this info file is written how many
            //features has been used to train the random forest
            FILE *info_rf;
            dest_filepath_to_random_forest_info_file.append(".info.bin");
            info_rf =fopen(dest_filepath_to_random_forest_info_file.c_str(),"wb");
            fwrite(&nr_of_features,sizeof(int),1,info_rf);
            fclose(info_rf);
            */

            //delete the arrays
            delete feature_array;
            delete label_array;
        }
    }
    else // fp==NULL
    {
        std::cout<<"train_boundary_rf(..) was not able to open "<<filepath_to_classification_set_file<<std::endl;
        exit(-1);
    }
}

/*! \fn save_boundary_probabilities(vigra::MultiArray<2,float> const & probability,
                                 std::string filepath_to_ws_region_image,
                                 std::string path_to_output_folder,
                                 vigra::BasicImage<unsigned int> & ws_region_image,
                                 marray::Marray<unsigned int> & one_boundings,
                                 marray::Marray<unsigned int> & two_boundings,
                                 std::vector< std::vector<point> > & arcs,
                                 std::vector<point> & junctions,
                                 int & dim_x,
                                 int & dim_y,
                                 std::string param_file_name
                                 )
 * \brief Write boundary probabilities into a file.
 * \param probability 2-dimensional marray representing the probabilities
 * \param filepath_to_ws_region_image Filepath to the watershed segmentation image
 * \param path_to_output_folder Path to the output folder
 * \param ws_region_image Watershed segmentation image
 * \param one_boundings Marray representing one boundings
 * \param two_boundings Marray representing two boundings
 * \param arcs Vector containing the segmentation arcs
 * \param junctions Vector containing the segmentation junctions
 * \param dim_x Image width
 * \param dim_y Image height
 * \param param_file_name Filename of the parameter file
 */ 
void save_boundary_probabilities(vigra::MultiArray<2,float> const & probability,
                                 std::string filepath_to_ws_region_image,
                                 std::string path_to_output_folder,
                                 vigra::BasicImage<unsigned int> & ws_region_image,
                                 marray::Marray<unsigned int> & one_boundings,
                                 marray::Marray<unsigned int> & two_boundings,
                                 std::vector< std::vector<point> > & arcs,
                                 std::vector<point> & junctions,
                                 int & dim_x,
                                 int & dim_y,
                                 std::string param_file_name
                                 )
{
    //IMPORT RESULTS FROM HDF5 file
    filepath_to_ws_region_image.resize(filepath_to_ws_region_image.size()-4);
    filepath_to_ws_region_image.append(".h5");
 
    load_cgp_data_structure(filepath_to_ws_region_image,
                            ws_region_image,
                            one_boundings,
                            two_boundings,
                            arcs,
                            junctions,
                            dim_x,
                            dim_y,
                            true);

    vigra::FImage * probability_image = new vigra::FImage[nr_of_classes];

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            for (int i=0; i<nr_of_classes; i++)
            {
                probability_image[i].resize(dim_x,dim_y);
                probability_image[i](x,y)=255;
            }
        }
    }
    
    //load validation training data
    int last_read_arc=-1;
    std::vector<int> vector_training_labels;
    vector_training_labels.resize((int)arcs.size(),10);

    std::string validation_file="boundary-labels/validation/";
    validation_file.append(get_filename(filepath_to_ws_region_image));
    validation_file.resize(validation_file.size()-3);
    validation_file.append(".dat");

    //string is read from temp file to check whether file is empty
    std::string teststring;

    std::ifstream training_file(validation_file.c_str());
    std::ifstream temp_training_file(validation_file.c_str());
    temp_training_file>>teststring;

    bool validation=false;
    if(training_file && teststring.size()!=0) validation=true;

    std::vector<int> class_count;
    std::vector<float> class_error;

    if(validation)
    {
        class_count.resize(nr_of_classes,0);
        class_error.resize(nr_of_classes,0.0f);

        //temp file is used to get one line out of the training file
        std::vector<int> temp;
        temp.resize(2);
        while(!training_file.eof())
        {
            //the vector of the training data is
            training_file>>temp[0];
            training_file>>temp[1];

            if(last_read_arc!=temp[0])
            {
                vector_training_labels[temp[0]]=temp[1];
                //HERE WE SET LAST READ TO THE LAST READED /WRITTEN POINTS TO AVOID THE LAST LINE
                //IS READED TWICE
                last_read_arc=temp[0];
                //we want to know how many data is from wich class
                class_count[temp[1]]++;
            }
        }
        training_file.close();
    }

    for(int a=0;a<(int)arcs.size();a++)
    {
        std::vector<point>  this_arc;
        this_arc=arcs[a];

        //now we loop over the points in this arc
        for(int p=0;p<(int)this_arc.size();p++)
        {
            int x=this_arc[p].x;
            int y=this_arc[p].y;
            for (int i=0; i<nr_of_classes; i++)
            {
                probability_image[i](x,y)=255*(1-probability(a,i));
            }
        }

        if(validation) for (int i=0; i<nr_of_classes; i++)
        {
            if (vector_training_labels[a]==i) class_error[i]=class_error[i]+((1-probability(a,i))/class_count[i]);
        }
    }

    if(validation)
    {
        std::string filepath_log_file="/export/home/tobinder/src/CIS/error_log.txt";
        std::ofstream log_file(filepath_log_file.c_str(), std::ios_base::out | std::ios_base::app);

        log_file <<param_file_name.c_str();

        for (int i=0; i<nr_of_classes; i++)
        {
            std::cout<<"Error Class "<<i<<": "<<class_error[i]<<std::endl;
            log_file <<" "<<class_error[i];
        }

        log_file << "\n";
        log_file.close();
    }

    std::string * filepath_prob = new std::string[nr_of_classes];

    for (int i=0; i<nr_of_classes; i++)
    {
        filepath_prob[i]=path_to_output_folder;

        std::stringstream s;
        s << i <<"/";
        filepath_prob[i].append(s.str());
        filepath_prob[i].append(param_file_name.c_str());
        filepath_prob[i].append(get_filename(filepath_to_ws_region_image));
        filepath_prob[i].resize(filepath_prob[i].size()-6);
        filepath_prob[i].append("jpg");

        exportImage(srcImageRange(probability_image[i]), vigra::ImageExportInfo(filepath_prob[i].c_str()));
    }
}

/*! \fn extract_boundary_probabilities(std::string filepath_to_feature_file,std::string path_to_ws_image,std::string filepath_to_random_forest_file,std::string path_to_output_folder,std::string path_to_gm_output_folder,std::string param_file)
 * \brief Extract boundary probabilities.
 * \param filepath_to_feature_file Filepath to the feature file
 * \param path_to_ws_image Filepath to the watershed segmentation image
 * \param filepath_to_random_forest_file Filepath to the random forest file
 * \param path_to_output_folder Path to the output folder
 * \param param_file Parameter file
 */
//fn is the filepath to the image file
void extract_boundary_probabilities(std::string filepath_to_feature_file,std::string path_to_ws_image,std::string filepath_to_random_forest_file,std::string path_to_output_folder,std::string path_to_gm_output_folder,std::string param_file)
{
    // some important parameters should not need to be changed in the source code
    // so let's load them from a parameter file
    ParameterFile paramFile;
    
    if( !paramFile.load(param_file) )
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    int nr_of_features=0;
    int nr_of_arcs=0;

    //First we open the random_forest file
    filepath_to_random_forest_file.append("_rf.hdf5");
    std::ifstream check_file(filepath_to_random_forest_file.c_str());
    if(!check_file)
    {
        std::cout<<"extract_boundary_probabilities(..)  Error: Random Forest File "<<filepath_to_random_forest_file<<" is NOT existend"<<std::endl;
        exit(-1);
    }
    check_file.close();

    vigra::RandomForest<> rf;
    rf_import_HDF5(rf,filepath_to_random_forest_file.c_str());

    //NOW FOR UNKNOWN IMAGES/FEATURES
    FILE *fp_unknown_features;
    fp_unknown_features =fopen(filepath_to_feature_file.c_str(),"rb");

    if(fp_unknown_features==NULL)
    {
        std::cout<<"extract_boundary_probabilities(..)  Error: Feature File "<<filepath_to_feature_file<<" is NOT existend"<<std::endl;
        exit(-1);
    }

    float float_nr_of_arcs=0;
    float float_nr_of_features=0;

    int read_in_arc_size=fread(&float_nr_of_arcs,sizeof(float),1,fp_unknown_features);
    if(read_in_arc_size!=1)
    {
        std::cout<<" extract_boundary_probabilities(..) error, could not read in arc size"<<std::endl;
        exit(-1);
    }

    int read_in_nr_of_features=fread(&float_nr_of_features,sizeof(float),1,fp_unknown_features);
    if(read_in_nr_of_features!=1)
    {
        std::cout<<"extract_boundary_probabilities(..) error, could not read in boundary feature size"<<std::endl;
        exit(-1);
    }

    //NOW WE KNOW HOW BIG THE unknow_feature_array has to be
    nr_of_features=(int)float_nr_of_features;
    nr_of_arcs=(int)float_nr_of_arcs;
    std::cout<<"nr of features:"<<nr_of_features<<std::endl;
    float * unknown_features_array = new float[nr_of_features*nr_of_arcs];

    int read_in=fread(unknown_features_array,sizeof(float),nr_of_features*nr_of_arcs,fp_unknown_features);
    if(read_in!=nr_of_features*nr_of_arcs)
    {
        std::cout<<"extract_boundary_probabilities() error, in unknown features array, could not load all the floats"<<std::endl;
        exit(-1);
    }

	for(size_t jj=0;jj<nr_of_arcs*nr_of_features;jj++)
	{
		if (std::isnan(unknown_features_array[jj]))
		{
			unknown_features_array[jj] = 0;
		}
	}

    fclose(fp_unknown_features);

    vigra::MultiArray<2, float> unknown_features(vigra::MultiArrayShape<2>::type(nr_of_arcs,nr_of_features),unknown_features_array);
    for(int x=0;x<nr_of_arcs;x++)
    {
        for(int y=0;y<nr_of_features;y++)
        {
            unknown_features(x,y)=unknown_features_array[nr_of_features*x+y];
        }
    }

    vigra::MultiArray<2, double> unknown_probability(vigra::MultiArrayShape<2>::type(nr_of_arcs,nr_of_classes));  //nr_of_classes is a global variable

    std::cout<<"Probability prediction...."<<std::endl;
    rf.predictProbabilities(unknown_features,unknown_probability);
    std::cout<<"....done"<<std::endl;

    std::string filepath_to_ws_image=path_to_ws_image;
    std::string feature_file_name_parameter=(get_filename(filepath_to_feature_file));
    //remove parameter name from feature file name
    char feature_file_name[80]="";
    feature_file_name_parameter.copy(feature_file_name,feature_file_name_parameter.size()-param_file_name.size(),param_file_name.size());
    filepath_to_ws_image.append(feature_file_name);
    //remove the ".bin"
    filepath_to_ws_image.resize(filepath_to_ws_image.size()-4);
    filepath_to_ws_image.append(".bmp");

    vigra::BasicImage<unsigned int> ws_region_image;
    marray::Marray<unsigned int> one_boundings;
    marray::Marray<unsigned int> two_boundings;
    std::vector< std::vector<point> > arcs;
    std::vector<point> junctions;
    int dim_x;
    int dim_y;

    delete unknown_features_array;

    save_boundary_probabilities(unknown_probability,
                                filepath_to_ws_image,
                                path_to_output_folder,
                                ws_region_image,
                                one_boundings,
                                two_boundings,
                                arcs,
                                junctions,
                                dim_x,
                                dim_y,
                                param_file_name
                                );

	Parameter<int> low_x;
	low_x.assign("", "low_x", 0);
    low_x.load(paramFile,"config");

	Parameter<int> low_y;
	low_y.assign("", "low_y", 0);
    low_y.load(paramFile,"config");

	Parameter<int> high_x;
	high_x.assign("", "high_x", dim_x-1);
    high_x.load(paramFile,"config");

	Parameter<int> high_y;
	high_y.assign("", "high_y", dim_y-1);
    high_y.load(paramFile,"config");

    vigra::MultiArray<2, double> boundary_probability(vigra::MultiArrayShape<2>::type(nr_of_arcs,2));
    for (int arcindex=0; arcindex<nr_of_arcs; arcindex++)
    {
        if ((arcs[arcindex][0].x > low_x  && arcs[arcindex][arcs[arcindex].size()-1].x > low_x )  &&
            (arcs[arcindex][0].x < high_x && arcs[arcindex][arcs[arcindex].size()-1].x < high_x)  &&
            (arcs[arcindex][0].y > low_y  && arcs[arcindex][arcs[arcindex].size()-1].y > low_y )  &&
            (arcs[arcindex][0].y < high_y && arcs[arcindex][arcs[arcindex].size()-1].y < high_y))
        {//both arc endpoint are within selected region
            boundary_probability(arcindex,1)=unknown_probability(arcindex,1);
            boundary_probability(arcindex,0)=unknown_probability(arcindex,0);
            for (int i=2;i<nr_of_classes;i++) boundary_probability(arcindex,0)+=unknown_probability(arcindex,i);
        }
        else
        {
            boundary_probability(arcindex,1)=0.0f;  
            boundary_probability(arcindex,0)=1.0f;
        }
    }

    /*std::cout<<"compute angles..."<<std::endl;
    std::vector<std::vector<float> > angles;
    compute_angles(filepath_to_ws_image, boundary_probability, one_boundings, arcs, junctions, angles, 9);
    std::cout<<"...done"<<std::endl;

    std::string filepath_to_gm_output=path_to_gm_output_folder;
    filepath_to_gm_output.append(get_filename(filepath_to_ws_image));
    filepath_to_gm_output.resize(filepath_to_gm_output.size()-4);

    std::cout<<"start graphical model predictions..."<<std::endl;
    buildGraphicalModel<float>(filepath_to_gm_output,
                               ws_region_image,
                               one_boundings,
                               two_boundings,
                               boundary_probability,
                               angles,
                               arcs,
                               dim_x,
                               dim_y,
                               0.75f,
                               20.0f
                               );

    std::cout<<"...done"<<std::endl;
*/
}
