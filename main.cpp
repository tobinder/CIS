/*! \file main.cpp
 * \brief CIS's main file.
 */
/*! \mainpage CIS Documentation
 * \tableofcontents
 *
 * \section section_licence Licence
 *
 *  Copyright (c) 2013 Tobias Binder.
 *  
 *  This software was developed at the University of Heidelberg by
 *  Tobias Binder, Bjoern Andres, Thorsten Beier and Arthur Kuehlwein.
 *  Enquiries shall be directed to tobias.binder@iwr.uni-heidelberg.de.
 * 
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met:
 * 
 *  - Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright notice, 
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  - All advertising materials mentioning features or use of this software must 
 *    display the following acknowledgement: ``This product includes the CIS
 *    package developed by Tobias Binder and others''.
 *  - The name of the author must not be used to endorse or promote products 
 *    derived from this software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
 *  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 *  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
 *  EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
 *  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
 *  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * \section synopsis Synopsis
 * The project Colour Image Segmentation (\b CIS) goes back to the project SubgrainBoundaries by Kai Ueltzhöffer, which enabled segmentation and classification of microstructure images of polar ice.  With the aim of creating a versatile segmentation framework, the project SubgrainBoundaries was reconstructed and extended to the project \b CIS by Thorsten Beier and Björn Andres.
 * Since March 2010, further developments have been made by Tobias Binder, which are twofold: conserving or even extending the "universality" of the segmentation algorithm on the one hand, developing an optional adaption towards microstructure images of ice on the other hand.
 * The choice of features is done in the parameter file, only the number of edge classes is a hardcoded global constant. Extensive tests of \b CIS have been conducted for black and white images only, but processing of coloured images is implemented as well.
 * \section cargs Command line arguments
 * The \b CIS's <tt>main.cpp</tt> contains the following command line arguments. Note that there are command line arguments with a suffix \b -gui. These are specifically designed to be used with a GUI and thus not meant to be explicitly called from a terminal.
 *
 * \subsection cargs_resize -resize
 * Doubles the size of an image and executes a B-spline interpolation. This function is yet to be used for ice data.
 * \subsection cargs_pixellabels -import-pixel-labels
 * Loads the training data that has been saved in the image - not implemented. Not needed for ice data.
 * \subsection cargs_pixellabels2 -pixel-labels
 * Loads existing training data and enables adding of training points via means of a graphical user interface afterwards.\n Usage: <tt>./cis -pixel-labels images/*.bmp images/as/ pixel-labels/ </tt>
 * \subsection cargs_pixelfeatures -pixel-features
 * Compute and save pixel features defined in the parameter file for an image.\n Usage: <tt>./cis -pixel-features images/*.bmp images/as/ pixel-features/ </tt>
 * \subsection cargs_pixelclassset -pixel-classification-set
 * Compute pixel features on the labeled training data.\n Usage: <tt>./cis -pixel-classification-set images/*.bmp images/as/ pixel-labels/ pixel-classification/name </tt>
 * \subsection cargs_pixelrftrain -pixel-rf-train
 * Train the random forest with the classification set.\n Usage: <tt>./cis -pixel-rf-train pixel-classification/name pixel-classification/random-forests/name </tt>
 * \subsection cargs_pixelrfpredict -pixel-rf-predict
 * Create a probability map from given pixel features and random forest.\n Usage: <tt> ./cis -pixel-rf-predict images/*.bmp pixel-features/ pixel-classification/random-forests/name pixel-classification</tt>
 * \subsection cargs_wsregions -ws-regions
 * Create a watershed segmentation from a given probability map and CGP-structure calculation.\n Usage: <tt>./cis -ws-regions pixel-classification/probability-maps/*.bmp pixel-classification/watershed-segmentation/ </tt>
 * \subsection cargs_wsman -ws-manual
 * Create a watershed segmentation with an externally generated probability map.\n Usage: <tt>./cis -ws-manual binary-maps/*.bmp pixel-classification/watershed-segmentation/ </tt>
 * \subsection cargs_ws -ws
 * Create a watershed segmentation from an externally preprocessed probability map.\n Usage: <tt>./cis -ws pixel-classification/probability-maps/*.bmp image/image_name.bmp pixel-classification/watershed-segmentation/ </tt>
 * \subsection cargs_cgpstruct -cgp-structure 
 * Compute the CGP boundary structure and save this combined with the CPG-structure calculation. \b Not \b necessary, actually included in \b -ws-regions.\n Usage: <tt>./cis -cgp-structure pixel-classification/watershed-segmentation/*.h5 pixel-classification/watershed-segmentation/ </tt>
 * \subsection cargs_boundarylabels -boundary-labels
 * Manual labeling of boundaries in different classes.\n Usage: <tt>./cis -boundary-labels images/*.bmp pixel-classification/watershed-segmentation/ boundary-labels/ </tt>
 * \subsection cargs_boundaryfeatures -boundary-features
 * Compute and save boundary features for an image.\n Usage: <tt>./cis -boundary-features images/*.bmp pixel-classification/watershed-segmentation/ pixel-classification/probability-maps/ pixel-features/ boundary-features/ </tt>
 * \subsection cargs_boundaryclassset -boundary-classification-set
 * Create a boundary training set.\n Usage: <tt>./cis -boundary-classification-set boundary-labels/*.bmp.dat boundary-features/ boundary-classification/trainingsets/name </tt>
 * \subsection cargs_boundaryrftrain -boundary-rf-train
 * Train the boundary random forest with the boundary training set.\n Usage: <tt>./cis -boundary-rf-train boundary-classification/trainingsets/name boundary-classification/random-forests/name </tt>
 * \subsection cargs_boundaryrfpred -boundary-rf-predict
 * Boundary random forest prediction.\n Usage: <tt>./cis -boundary-rf-predict boundary-features/*.bmp.bin pixel-classification/watershed-segmentation/ boundary-classification </tt>
 * \subsection cargs_junctions -junctions
 * Compute junction features for a training image. Create a classification set and train the junction random forest.\n Usage: <tt>./cis -junctions /images/image.bmp pixel-classification/watershed-segmentation/ boundary-labels/image.bmp.junctions.dat </tt>
 * \subsection cargs_bubbles -bubbles
 * Create a bubble image from two existing preprocessed images (*.as.bmp and *.de.bmp).\n Usage: <tt>./cis -bubbles images/*.bmp images/as/ </tt>
 * 
 * In addition to the above command line arguments, there are composite arguments, which use filepaths defined in the parameter file and combine some of the function calls.
 * \subsection cargs_pixelcomplim -pixel-complete-image
 * All functions up until the computation of the boundary features are executed as needed for ice data evaluation. The suffix "image" meaning, that multiple images can be processed in batch mode.\n Usage: <tt>./cis -pixel-complete-image images/*.bmp </tt>
 * \subsection cargs_pixelcomplpa -pixel-complete-param
 * The suffix "param" means, that multiple parameter files can be processed for an image in order to evaluate the influence of different parameters. The intermediate and end results are named after the parameter file used.\n Usage: <tt>./cis -pixel-complete-param images/image.bmp parameterfiles/*.txt </tt>
 * \subsection cargs_boundarycomplsuf -pixel-complete-suffix and -boundary-complete-suffix
 * The "suffix" function can be used when there exists a large amount of images in different folders and the intermediate and end results are to be saved in different folders as well. It appends the suffix on all paths. That way, the paths do not have to be adjusted in the parameter file.\n Usage: <tt>./cis -pixel-complete-suffix suffix/ images/*.bmp </tt> or <tt>./cis -boundary-complete-suffix suffix/ images/*.bmp </tt>
 * \subsection cargs_boundarycomplim -boundary-complete-image
 * Analogous to the pixel composite function \b -pixel-complete-image, there exists a boundary composite function. The suffix "image" again indicating image batch mode.\n Usage: <tt>./cis -boundary-complete-image images/*.bmp </tt>
 * \subsection cargs_boundarycomplpa -boundary-complete-param
 * Boundary composite function with parameter file in batch mode, analogous to \b -pixel-complete-param.\n Usage: <tt>./cis -boundary-complete-param images/image.bmp parameterfiles/*.txt </tt>
 * 
 * \section includes Libraries
 * This project uses the following (C++) libraries:
 * - <a href="http://www.cimg.sourceforge.net/">CIMG</a> enables easy use, processing and display of images. 
 * - <a href="http://hci.iwr.uni-heidelberg.de/vigra/">VIGRA</a> is a superb collection of algorithms developed in the HCI. It includes algorithms for image analysis and classification with random forests.
 * - <a href="https://gorgonzola.iwr.uni-heidelberg.de/svn/CGP">CGP</a> enables calculating an image's geometry based on its segmentation. It has been specifically developed for big volume data, so small modifications are
 * necessary to use it for two-dimensional data. It has been developed in the HCI. Note that SVN rights are required for web-access.
 * - <a href="http://www.boost.org/">Boost</a> is a free library consisting of a multitude of portable sublibraries, which serve a number of tasks. The library <B>Boost</B>.MPI enables parallelisation of GCP-computations.
 * - <a href="http://www.hdfgroup.org/">HDF5</a> is a data format used for storing scientific data efficiently and flexible. It is used for random forests at the moment, but further data can be integrated.
 */
#include <cgp/cgp_config.hxx>

#if HAVE_MPI
#include "mpi.h"
#endif

#include <vigra/colorconversions.hxx>
#include <vigra/convolution.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/impex.hxx>
#include <vigra/edgedetection.hxx>

#include <iostream>
#include <string>

#include "bubbles.h"

namespace std
{
    unsigned int abs(unsigned int u)
    {
	    return u;
    }
}

const int nr_of_classes=3;
//Static class attribute 
std::string ParameterFile::filepath = "";

#include "ParameteredObject.hxx"

#include "path_functions.h"
#include "pixel_training.h"
#include "feature_extractor.h"
#include "feature_extractor_gray.h"
#include "pixel_probabilities.h"
#include "compute_watershed_regions.h"

#include "boundary_training.h"
#include "boundary_features.h"
#include "boundary_features_gray.h"
#include "boundary_probabilities.h"
#include "junction_features.h"

#include "pixel_complete.h"
#include "boundary_complete.h"

#include "ElleSim.h"

SplitStream sout(std::cout);

int main(int argc, char *argv[])
{
    if (argc<1)
    {
        std::cout<<"you need more arguments for the CIS tool"<<std::endl;
        return 0;
    }
    else
    {
        std::string test=argv[argc-1];
        if (test == "-test")
        {
            return 0;
        }

        int size=1;
        int rank=0;

        #if HAVE_MPI
            MPI_Init(&argc, &argv);
            MPI_Comm_size(MPI_COMM_WORLD, &size);  
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        #endif //HAVE_MPI

        if (rank==0) std::cout<<std::endl<<"CIS -Color Image Segmentation"<<std::endl<<std::endl;

        std::string command_line_option=argv[1];

        /*
        OPTION -resize  (used to double the image size an to do b-spline interpolation)
        */
        if(command_line_option=="-resize")
        {  
            std::cout<<"Option -resize"<<std::endl;

            std::string dest_path=argv[argc-1];
            std::cout<<"Path to save resized images: "<<dest_path<<std::endl;
            int i=2;
            while(i<(argc-1))
            {
                std::string path_and_filename_of_image=argv[i];
                std::string filename_of_image=get_filename(path_and_filename_of_image);
                std::string dest_filepath=dest_path;

                filename_of_image.resize(filename_of_image.size()-4);
                filename_of_image.append(".bmp");
                dest_filepath.append(filename_of_image);

                vigra::ImageImportInfo info(argv[i]);

                if(info.isGrayscale())
                {
                    vigra::BImage in(info.width(),info.height());
                    importImage(info, destImage(in));
                    vigra::BImage out(0.2*info.width(),0.2*info.height());
                    resizeImageSplineInterpolation(srcImageRange(in),destImageRange(out));

                    exportImage(srcImageRange(out), vigra::ImageExportInfo( dest_filepath.c_str() ));
                    i++;
                }
                else
                {
                    vigra::BRGBImage in(info.width(),info.height());
                    importImage(info, destImage(in));
                    vigra::BRGBImage out(2*info.width(),2*info.height());
                    resizeImageSplineInterpolation(srcImageRange(in),destImageRange(out));

                    exportImage(srcImageRange(out), vigra::ImageExportInfo( dest_filepath.c_str() ));

                    i++;
                }
                std::cout<<std::endl;
            }
        }

        /*
        OPTION -import-pixel-labels  (used to load trainingsdata saved in an image)
        */
        else if(command_line_option=="-import-pixel-labels")
        {
            std::cout<<"Option -import-pixel-labels"<<std::endl;

            std::string dest_path=argv[argc-1];
            std::cout<<"Path to save pixel-labels: "<<dest_path<<std::endl;
            int i=2;
            while(i<(argc-1))
            {
                std::string path_and_filename_of_image=argv[i];
                std::string path_of_image=get_path(path_and_filename_of_image);
                std::string filename_of_image=get_filename(path_and_filename_of_image);

                std::cout<<"path of the image: "<<path_of_image<<std::endl;
                std::cout<<"filename of the image : "<<filename_of_image<<std::endl;
                import_and_interpolate_handmade_segmentatations(argv[i],dest_path);
                i++;
                std::cout<<i<<"/"<<argc-1<<std::endl;
            }
        }

        /*
        OPTION -pixel-labels  (used to get trainingsdata for an image
                               it writes x_coord, y_coord, label in a *.dat file)  
        */
        else if(command_line_option=="-pixel-labels")
        {
            std::cout<<"Option -pixel-labels"<<std::endl;

            std::string dest_path=argv[argc-1];
            std::cout<<"Path to save pixel-labels: "<<dest_path<<std::endl;
            std::string as_path=argv[argc-2];   
            int i=2;
            while(i<(argc-2))
            {
                std::string path_and_filename_of_image=argv[i];
                std::string path_of_image=get_path(path_and_filename_of_image);
                std::string filename_of_image=get_filename(path_and_filename_of_image);

                std::string as_image_filepath=as_path;
                as_image_filepath.append(filename_of_image);
                as_image_filepath.append(".as.bmp");

                std::cout<<"path of the image: "<<path_of_image<<std::endl;
                std::cout<<"filename of the image: "<<filename_of_image<<std::endl;
                std::cout<<"filepath anisotropic smoothed: "<<as_image_filepath<<std::endl;

                Pixel_training_image pixel_train_image(argv[i],as_image_filepath,dest_path);
                pixel_train_image.file_to_training_data();//load existing trainingsdata
                //pixel_train_image.do_analysis();
                pixel_train_image.do_training();
                i++;
                std::cout<<std::endl;
            }
        }

        /*
        OPTION -pixel-labels-gui
        */
        else if(command_line_option=="-pixel-labels-gui")
        {
            std::string segmentation_path=argv[argc-3];
            std::string prob_map_path=argv[argc-4];
            std::string dest_path=argv[argc-5];
            std::string as_path=argv[argc-6];

            std::string folder=argv[argc-2];
            if (folder=="no") folder="";
            prob_map_path.append(folder.c_str());
            segmentation_path.append(folder.c_str());
            as_path.append(folder.c_str());

            std::cout<<"Path to save pixel-labels: "<<dest_path<<std::endl;
            if (prob_map_path!="no") std::cout<<"Path to probability maps: "<<prob_map_path<<std::endl;
            if (segmentation_path!="no") std::cout<<"Path to segmentation: "<<segmentation_path<<std::endl;
            int i=2;
            while(i<(argc-6))
            {
                std::string path_and_filename_of_image=argv[i];
                std::string path_of_image=get_path(path_and_filename_of_image);
                std::string filename_of_image=get_filename(path_and_filename_of_image);

                std::string as_image_filepath=as_path;
                as_image_filepath.append(filename_of_image);
                as_image_filepath.append(".as.bmp");

                std::cout<<"path of the image: "<<path_of_image<<std::endl;
                std::cout<<"filename of the image: "<<filename_of_image<<std::endl;
                std::cout<<"filepath anisotropic smoothed: "<<as_image_filepath<<std::endl;

                std::string filepath_prob_map=prob_map_path;
                if (prob_map_path!="no") filepath_prob_map.append(filename_of_image);
                else filepath_prob_map="";
                std::string filepath_segmentation=segmentation_path;
                if (prob_map_path!="no") filepath_segmentation.append(filename_of_image);
                else filepath_segmentation="";

                Pixel_training_image pixel_train_image(argv[i],as_image_filepath,dest_path,argv[argc-1],
                    filepath_prob_map,filepath_segmentation);
                pixel_train_image.file_to_training_data();//load existing trainingsdata               
                pixel_train_image.do_training();
                i++;
                std::cout<<std::endl;
            }
        }
	  
        /*
        OPTION -pixel-features  (used to compute the features of an image)
        */
        else if(command_line_option=="-pixel-features")
        {
            std::cout<<"Option -pixel-features"<<std::endl;
            std::string as_path=argv[argc-2];
            std::string dest_path=argv[argc-1];
            std::cout<<"Path to save pixel-features: "<<dest_path<<std::endl;
            int i=2;
            while(i<(argc-2))
            {
                std::string path_and_filename_of_image=argv[i];
                std::string path_of_image=get_path(path_and_filename_of_image);
                std::string filename_of_image=get_filename(path_and_filename_of_image);

                std::string as_image_filepath=as_path;
                as_image_filepath.append(filename_of_image);
                as_image_filepath.append(".as.bmp");

                std::cout<<"path of the image: "<<path_of_image<<std::endl;
                std::cout<<"filename of the image: "<<filename_of_image<<std::endl;
                std::cout<<"filepath anisotropic smoothed: "<<as_image_filepath<<std::endl;

                //Check if grayscale or color
                vigra::ImageImportInfo info(argv[i]);
                if(info.isGrayscale())
                {
                    FeatureExtractorGray pixel_features(argv[i],as_image_filepath,dest_path);
                    pixel_features.extract(true,"parameters.txt");
                }
                else
                {
                    FeatureExtractor pixel_features(argv[i],as_image_filepath,dest_path);
                    pixel_features.extract(true,"parameters.txt");
                }
                i++;
                std::cout<<std::endl;
            }
        }

        /*
        OPTION -pixel-classification-set  (used to compute the features on the labled training data
                                            it appends the data to one big *.0.bin and *.1.bin file)
        */
        else if(command_line_option=="-pixel-classification-set")
        {
            if(argc>=4)
            {
                std::cout<<"Option -pixel-classification-set"<<std::endl;
                std::string dest_filepath=argv[argc-1];
                std::cout<<"Filepath to save pixel-classification-set: "<<dest_filepath<<std::endl;
                std::string path_to_labels=argv[argc-2];
                std::string as_path=argv[argc-3];

                int i=2;
                while(i<(argc-3))
                {
                    std::string path_and_filename_of_image=argv[i];
                    std::string path_of_image=get_path(path_and_filename_of_image);
                    std::string filename_of_image=get_filename(path_and_filename_of_image);

                    std::string as_image_filepath=as_path;
                    as_image_filepath.append(filename_of_image);
                    as_image_filepath.append(".as.bmp");

                    std::cout<<"path of the image: "<<path_of_image<<std::endl;
                    std::cout<<"filename of the image: "<<filename_of_image<<std::endl;
                    std::cout<<"filepath anisotropic smoothed: "<<as_image_filepath<<std::endl;

                    //Check if grayscale or color
                    vigra::ImageImportInfo info(argv[i]);
                    if(info.isGrayscale())
                    {
                        FeatureExtractorGray pixel_features(argv[i],as_image_filepath,"");
                        int nr_of_features = pixel_features.extract(false,"parameters.txt");
                        pixel_features.extract_training_data(path_to_labels,dest_filepath,nr_of_features);
                    }
                    else
                    {
                        FeatureExtractor pixel_features(argv[i],as_image_filepath,"");
                        int nr_of_features = pixel_features.extract(false,"parameters.txt");
                        pixel_features.extract_training_data(path_to_labels,dest_filepath,nr_of_features);
                    }
                    i++;
                    std::cout<<std::endl;
                }
            }
            else
            {
                std::cout<<"you need more arguments for CIS -pixel-classification-set"<<std::endl;
            }
        }

        /*
        OPTION -pixel-classification-gui
        */
        else if(command_line_option=="-pixel-classification-gui")
        {
            std::string as_path=argv[argc-6];
            std::string path_to_labels=argv[argc-5];
            std::string dest_filepath=argv[argc-4];
            dest_filepath.append(argv[argc-3]);
            std::cout<<"Filepath to save pixel-classification-set: "<<dest_filepath<<std::endl;

            bool overwrite;
            if (atoi(argv[argc-1])==1) overwrite = true;
            else overwrite = false;

            int i=2;
            while(i<(argc-6))
            {
                std::string path_and_filename_of_image=argv[i];
                std::string path_of_image=get_path(path_and_filename_of_image);
                std::string filename_of_image=get_filename(path_and_filename_of_image);

                std::string as_image_filepath=as_path;
                as_image_filepath.append(filename_of_image);
                as_image_filepath.append(".as.bmp");

                std::cout<<"path of the image: "<<path_of_image<<std::endl;
                std::cout<<"filename of the image: "<<filename_of_image<<std::endl;
                std::cout<<"filepath anisotropic smoothed: "<<as_image_filepath<<std::endl;

                //Check if grayscale or color
                vigra::ImageImportInfo info(argv[i]);
                if(info.isGrayscale())
                {
                    FeatureExtractorGray pixel_features(argv[i],as_image_filepath,"");
                    int nr_of_features = pixel_features.extract(false,argv[argc-2]);
                    pixel_features.extract_training_data(path_to_labels,dest_filepath,nr_of_features,overwrite);
                }
                else
                {
                    FeatureExtractor pixel_features(argv[i],as_image_filepath,"");
                    int nr_of_features = pixel_features.extract(false,argv[argc-2]);
                    pixel_features.extract_training_data(path_to_labels,dest_filepath,nr_of_features);
                }
                i++;
                std::cout<<std::endl;
            }
        }

        /*
        OPTION -pixel-rf-train  (used to train the random-forest with the classification-set)
        */
        else if(command_line_option=="-pixel-rf-train")
        {
            std::cout<<"Option -rf-train"<<std::endl;

            std::string filepath_to_classification_set_file=argv[argc-2];
            std::cout<<"Filepath to classification set file: "<<filepath_to_classification_set_file<<std::endl;
            std::string dest_filepath_to_random_forest_file=argv[argc-1];
            train_pixel_rf(filepath_to_classification_set_file,dest_filepath_to_random_forest_file,"parameters.txt");
        }

        /*
        OPTION -pixel-rf-train-gui
        */
        else if(command_line_option=="-pixel-rf-train-gui")
        {
            std::string filepath_to_classification_set_file=argv[argc-4];
            filepath_to_classification_set_file.append(argv[argc-2]);
            std::cout<<"Filepath to classification set file: "<<filepath_to_classification_set_file<<std::endl;
            std::string dest_filepath_to_random_forest_file=argv[argc-3];
            dest_filepath_to_random_forest_file.append(argv[argc-2]);
            train_pixel_rf(filepath_to_classification_set_file,dest_filepath_to_random_forest_file,argv[argc-1]);
        }

        /*
        OPTION -pixel-rf-predict  (used to create a probability-map
                                   from given pixel-features and random-forest)
        */
        else if(command_line_option=="-pixel-rf-predict")
        {
            std::cout<<"Option -rf-predict"<<std::endl;

            std::string path_to_feature_files=argv[argc-3];
            std::cout<<"path to feature files: "<<path_to_feature_files<<std::endl;

            std::string filepath_to_random_forest_file=argv[argc-2];
            std::cout<<"Filepath to random forest file: "<<filepath_to_random_forest_file<<std::endl;

            std::string dest_path_probability_maps=argv[argc-1];
            std::cout<<"Dest. path to probability maps : "<<dest_path_probability_maps<<std::endl;

            int i=2;
            while(i<argc-3)
            {
                std::cout<<"NR: "<<i-1<<"/"<<argc-5<<std::endl;
                extract_pixel_probabilities(argv[i],path_to_feature_files,filepath_to_random_forest_file,dest_path_probability_maps,
                    "parameters.txt");
                i++;
            }

        }

        /*
        OPTION -ws-regions  (used to create watershed-segmentation
                             which is saved als *.bmp and *.bmp.h5 file)
        */
        else if(command_line_option=="-ws-regions")
        {
            std::string dest_image_path=argv[argc-1];
            int i=2;

            int threshold=26;
//            std::cout<<"threshold: "<<std::endl;
//            std::cin>>threshold;
//            std::cout<<"threshold="<<threshold<<std::endl;
 
            double scale=1;
//            std::cout<<"scale: "<<std::endl;
//            std::cin>>scale;
//            std::cout<<"scale="<<scale<<std::endl;
 
            int equalTolerance=1;
//            std::cout<<"equalTolerance: "<<std::endl;
//            std::cin>>equalTolerance;
//            std::cout<<"equalTolerance="<<equalTolerance<<std::endl;

            std::cout<<".0.bmp images will be ignored"<<std::endl;

            while(i<argc-1)
            {
                std::cout<<"nr: "<<i-1<<" of " <<argc-3<<std::endl;
                std::string source_image_path=argv[i];
                char check_ending[6]="";
                char compare[7]="";                
                strcpy(compare,".0.bmp");
                source_image_path.copy(check_ending,6,source_image_path.size()-6);
                if (strcmp(check_ending, compare)==0) std::cout<<source_image_path<<" ignored"<<std::endl;

                std::string source_image_ext = source_image_path;
                source_image_ext.append(".eed_ms.bmp");
                FILE * eed_ext;
                eed_ext = fopen(source_image_ext.c_str(), "r");

                if(eed_ext != NULL)
                {
                    compute_watershed_regions(source_image_ext,dest_image_path,threshold,scale,equalTolerance,size,rank);
                    i++;
                }
                else
                {
                    compute_watershed_regions(source_image_path,dest_image_path,threshold,scale,equalTolerance,size,rank);
                    i++;
                }
            }
 
        }

        /*
        OPTION -ws  (used to create watershed-segmentation with extern probability-map)
        */
        else if(command_line_option=="-ws")
        {
            std::string source_image_filepath=argv[argc-2];
            std::string dest_image_path=argv[argc-1];
            int i=2;
            
            int threshold;
            std::cout<<"threshold: "<<std::endl;
            std::cin>>threshold;
            std::cout<<"threshold="<<threshold<<std::endl;

            int equalTolerance;
            std::cout<<"equalTolerance: "<<std::endl;
            std::cin>>equalTolerance;
            std::cout<<"equalTolerance="<<equalTolerance<<std::endl;

            while(i<argc-2)
            {
                std::cout<<"nr: "<<i-1<<" of " <<argc-4<<std::endl;
                std::string preprocessed_probmap_filepath=argv[i];
                compute_ws_regions(source_image_filepath,preprocessed_probmap_filepath,dest_image_path,threshold,equalTolerance);
                i++;
            }
        }

        /*
        OPTION -ws-manual  (used to create watershed-segmentationfrom manually generated binary map)
        */
        else if(command_line_option=="-ws-manual")
        {
            std::string dest_image_path=argv[argc-2];

            std::string folder=argv[argc-1];
            if (folder=="no") folder="";
            dest_image_path.append(folder.c_str());

            int equalTolerance=1;

            int i=2;
            while(i < argc-2)
            {
                std::cout << "nr: " << i-1 << " of " << argc-4 << std::endl;
                std::string source_image_path = argv[i];
                
                //Check if image is grayscale or color
                vigra::ImageImportInfo info(source_image_path.c_str());
                if(info.isGrayscale())
                {
                    compute_watershed_regions_binary(source_image_path,dest_image_path,equalTolerance,size,rank);
                }
                else
                {
                    compute_watershed_regions_binary_color(source_image_path,dest_image_path,equalTolerance,size,rank);                    
                }
                i++;
            }
        }

        /*
        OPTION -cgp-structure  (used to compute the cgp boundary-structure)
        */
        else if(command_line_option=="-cgp-structure")
        {
            #if HAVE_MPI
                if(size < 2)
                {
                    std::cout<<"Notice: CIS has been compiled with MPI support"<<std::endl;
                    std::cout<<"but only one process has been started."<<std::endl;
                    std::cout<<"Try to run CIS like this: mpirun -np 16 ./cis -cgp-structure"<<std::endl;
                }
                if(argc>4 && rank==0 && size>1)
                {
                    std::cout<<"Option -cgp-structure with MPI support handles only single files,"<<std::endl;
                    std::cout<<"use a shell script instead or start only one process."<<std::endl;
                    exit(-1);
                }
            #endif //HAVE_MPI

            if (rank==0) std::cout<<"Option -cgp-structure"<<std::endl;

            std::string dest_path=argv[argc-1];
            if (rank==0) std::cout<<"Dest. path to cgp structure file : "<<dest_path<<std::endl;

            int i=2;
            while(i<argc-1)
            {
                if (rank==0) std::cout<<"Image nr: "<<i-1<<" of " <<argc-3<<std::endl;
                std::string ws_image_filepath=argv[i];

                std::string dest_cgp_path=dest_path;
                dest_cgp_path.append(get_filename(ws_image_filepath));
                dest_cgp_path.resize(dest_cgp_path.size()-3);
                //ws_image_filepath.append(".vec");

                compute_and_save_cgp_data_structure(dest_cgp_path,size,rank);
                i++;
            }
        }

        /*
        OPTION -pixel-complete-param  (used when labeling for at least one image is given to do all steps
                                       until ws-regions use standard path, less arguments
                                       BATCH mode for parameter files)
        */
        else if(command_line_option=="-pixel-complete-param")
        {
            std::cout<<"Option -pixel-complete-param"<<std::endl;

            // some important parameters should not need to be changed in the source code
            // so let's load them from a parameter file
            ParameterFile paramFile;

            int i=3;
            while(i<argc)
            {
                bool default_rf_trained=false;
    
                std::cout<<"Parameter file "<<i-2<<" of " <<argc-3<<": "<<argv[i]<<std::endl;

                if( !paramFile.load(argv[i]) )
                {
                    std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                    return 0;
                }
                do_pixel_complete(paramFile, argv[2], argv[i], &default_rf_trained,size,rank);
                i++;
            }
        }

        /*
        OPTION -pixel-complete-image  (used when labeling for at least one image is given to do all steps
                                       until ws-regions use standard path, less arguments
                                       BATCH mode for images)
        */
        else if(command_line_option=="-pixel-complete-image")
        {
            std::cout<<"Option -pixel-complete-image"<<std::endl;

            // some important parameters should not need to be changed in the source code
            // so let's load them from a parameter file
            ParameterFile paramFile;

            if( !paramFile.load("parameters.txt") )//default parameter file for image batch mode
            {
                std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                return 0;
            }

            //bool default_rf_trained=false;
            bool default_rf_trained=true;

            int i=2;
            while(i<argc)
            {
                std::cout<<"Image nr: "<<i-1<<" of " <<argc-2<<": "<<argv[i]<<std::endl;
                do_pixel_complete(paramFile, argv[i], "", &default_rf_trained,size,rank);
                i++;
            }
        }

        /*
        OPTION -pixel-complete-gui  (segmentation of single images or batches via gui)
        */
        else if(command_line_option=="-pixel-complete-gui")
        {
            // some important parameters should not need to be changed in the source code
            // so let's load them from a parameter file
            ParameterFile paramFile;

            if( !paramFile.load(argv[argc-3]) )//parameter file defined in gui
            {
                std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                return 0;
            }

            std::string path_as_image=argv[2];
            std::string path_pixel_feature=argv[3];
            std::string path_prob_map=argv[4];
            std::string path_watershed=argv[5];
            std::string path_boundary_feature=argv[6];
            std::string path_thumbs=argv[argc-1];

            std::string folder=argv[argc-2];
            if (folder=="no") folder="";
            path_as_image.append(folder.c_str());

            if(path_pixel_feature.compare(path_pixel_feature.size()-5, 5, "FALSE") == 0)
            {
                path_pixel_feature.resize(path_pixel_feature.size()-5);
                path_pixel_feature.append(folder.c_str());
                path_pixel_feature.append("FALSE");
            }
            else if(path_pixel_feature.compare(path_pixel_feature.size()-4, 4, "TRUE") == 0)
            {
                path_pixel_feature.resize(path_pixel_feature.size()-4);
                path_pixel_feature.append(folder.c_str());
                path_pixel_feature.append("TRUE");
            }
            else path_pixel_feature.append(folder.c_str());
            path_prob_map.append(folder.c_str());
            path_watershed.append(folder.c_str());
            path_boundary_feature.append(folder.c_str());

            Parameter<std::string> p_thumbs;
            p_thumbs.assign("", "path_thumbs", "no");
            if (path_thumbs!="no") path_thumbs.append(folder.c_str());
            p_thumbs=path_thumbs;
            p_thumbs.save(paramFile,"config");

            int i=7;
            while(i<argc-3)
            {
                std::cout<<"Image nr: "<<i-6<<" of " <<argc-10<<": "<<argv[i]<<std::endl;
                do_pixel_complete(paramFile, argv[i], argv[argc-3], path_as_image, path_pixel_feature, path_prob_map, path_watershed,
                                  path_boundary_feature, size, rank);
                i++;
            }
        }

        /*
        OPTION -pixel-complete-suffix  (similar to pixel-complete-image, using suffix for all folders specified in parameter file)
        */
        else if(command_line_option=="-pixel-complete-suffix")
        {
            std::cout<<"Option -pixel-complete-suffix"<<std::endl;

            // some important parameters should not need to be changed in the source code
            // so let's load them from a parameter file
            ParameterFile paramFile;

            if( !paramFile.load("parameters.txt") )//default parameter file for image batch mode
            {
                std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                return 0;
            }

            std::string suffix=argv[2];
            std::string buffer;

            Parameter<std::string> fp_as_image;
            fp_as_image.assign("", "filepath_as_image", "as/");
            fp_as_image.load(paramFile,"config");
            buffer=fp_as_image;
            buffer.append(suffix.c_str());
            fp_as_image=buffer;
            fp_as_image.save(paramFile,"config");

            Parameter<std::string> p_pixel_feature;
            p_pixel_feature.assign("", "path_pixel_feature", "pixel-features/");
            p_pixel_feature.load(paramFile,"config");
            buffer=p_pixel_feature;
            buffer.append(suffix.c_str());
            p_pixel_feature=buffer;
            p_pixel_feature.save(paramFile,"config");

            Parameter<std::string> path_prob_map;
            path_prob_map.assign("", "path_prob_map", "pixel-classification/probability-maps/");
            path_prob_map.load(paramFile,"config");
            buffer=path_prob_map;
            buffer.append(suffix.c_str());
            path_prob_map=buffer;
            path_prob_map.save(paramFile,"config");

            Parameter<std::string> p_watershed;
            p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
            p_watershed.load(paramFile,"config");
            buffer=p_watershed;
            buffer.append(suffix.c_str());
            p_watershed=buffer;
            p_watershed.save(paramFile,"config");

            Parameter<std::string> p_boundary_feature;
            p_boundary_feature.assign("", "path_boundary_feature", "boundary-features/");
            p_boundary_feature.load(paramFile,"config");
            buffer=p_boundary_feature;
            buffer.append(suffix.c_str());
            p_boundary_feature=buffer;
            p_boundary_feature.save(paramFile,"config");

            bool default_rf_trained=true;

            int i=3;
            while(i<argc)
            {
                std::cout<<"Image nr: "<<i-2<<" of " <<argc-3<<": "<<argv[i]<<std::endl;
                do_pixel_complete(paramFile, argv[i], "", &default_rf_trained,size,rank);
                i++;
            }
        }

        /*
        OPTION -boundary-features  (used to compute the boundary-features of an image)
        */
        else if(command_line_option=="-boundary-features")
        {
            std::cout<<"Option -boundary-features"<<std::endl;

            // some important parameters should not need to be changed in the source code
            // so let's load them from a parameter file
            ParameterFile paramFile;

            if( !paramFile.load("parameters.txt") )//default parameter file for image batch mode
            {
                std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                return 0;
            }

            //std::string source_path_image         =argv[argc-4];
            std::string source_ws_image_path=argv[argc-4];
            std::string source_path_pixel_features=argv[argc-2];
            std::string source_path_prob_map      =argv[argc-3];
            std::string dest_path_boundary_features=argv[argc-1];

            Parameter<std::string> p_thumbs;
            p_thumbs.assign("", "path_thumbs", "no");
            p_thumbs.load(paramFile,"config");
            std::string path_thumbs=p_thumbs;

            std::string suffix=dest_path_boundary_features;
            suffix.resize(suffix.size()-1);
            suffix=get_filename(suffix);
            if (atoi(suffix.c_str())>0)
            {
                suffix.append("/");
                std::cout<<"Found suffix: "<<suffix<<std::endl;
            }
            else suffix="";
            
            if (path_thumbs!="no")
            {
                path_thumbs.append(suffix.c_str());
            }

            int i=2;
            while(i<argc-4)
            {
                //Check if grayscale or color
                vigra::ImageImportInfo info(argv[i]);

                if(info.isGrayscale())
                {
                    boundary_features_gray(argv[i],source_ws_image_path,source_path_prob_map,source_path_pixel_features,
                        dest_path_boundary_features,"parameters.txt",path_thumbs);
                }
                else
                {
                    boundary_features(argv[i],source_ws_image_path,source_path_prob_map,source_path_pixel_features,
                        dest_path_boundary_features);
                }
                i++;
            }
        }

        /*
        OPTION -boundary-features-gui
        */
        else if(command_line_option=="-boundary-features-gui")
        {
            std::string source_ws_image_path=argv[argc-7];
            std::string source_path_pixel_features=argv[argc-5];
            std::string source_path_prob_map      =argv[argc-6];
            std::string dest_path_boundary_features=argv[argc-4];
            std::string path_thumbs=argv[argc-1];

            std::string folder=argv[argc-2];
            if (folder=="no") folder="";
            source_ws_image_path.append(folder.c_str());
            source_path_pixel_features.append(folder.c_str());
            source_path_prob_map.append(folder.c_str());
            dest_path_boundary_features.append(folder.c_str());
            if (path_thumbs!= "no") path_thumbs.append(folder.c_str());

            int i=2;
            while(i<argc-7)
            {
                //Check if grayscale or color
                vigra::ImageImportInfo info(argv[i]);
                if(info.isGrayscale())
                {
                    boundary_features_gray(argv[i],source_ws_image_path,source_path_prob_map,source_path_pixel_features,
                        dest_path_boundary_features,argv[argc-3],path_thumbs);
                }
                else
                {
                    boundary_features(argv[i],source_ws_image_path,source_path_prob_map,source_path_pixel_features,
                        dest_path_boundary_features);
                }
                i++;
            }
        }

        /*
        OPTION -boundary-labels  (used to label boundaries in different classes)
        */
        else if(command_line_option=="-boundary-labels")
        {
            std::cout<<"Option -boundary-labels"<<std::endl;

            //some strings for the filepath and path of the porb.map and the watershed image and the output
            std::string path_to_watershed_images=argv[argc-2];
            std::string path_to_training_data_output=argv[argc-1];

            std::cout<<"path ws: "<<path_to_watershed_images<<std::endl;
            std::cout<<"path out: "<<path_to_training_data_output<<std::endl;

            int i=2;
            while(i<argc-2)
            {
                std::string filename_of_the_image=get_filename(argv[i]);
                std::cout<<"filename of the image: "<<filename_of_the_image<<std::endl;

                std::string watershed_image_filepath=path_to_watershed_images;
                watershed_image_filepath.append(filename_of_the_image);
                watershed_image_filepath.append(".vec");

                //NOW WE CAN START THE TRAINING
                Boundary_training_image boundary_train_image(argv[i],watershed_image_filepath,path_to_training_data_output);
                boundary_train_image.file_to_training_data();
                boundary_train_image.do_training();
                i++;
            }
        }

        /*
        OPTION -boundary-labels-gui
        */
        else if(command_line_option=="-boundary-labels-gui")
        {
            //some strings for the filepath and path of the porb.map and the watershed image and the output
            std::string path_to_watershed_images=argv[argc-3];
            std::string path_to_training_data_output=argv[argc-2];

            std::cout<<"path ws: "<<path_to_watershed_images<<std::endl;
            std::cout<<"path out: "<<path_to_training_data_output<<std::endl;

            int i=2;
            while(i<argc-3)
            {
                std::string filename_of_the_image=get_filename(argv[i]);
                std::cout<<"filename of the image: "<<filename_of_the_image<<std::endl;

                std::string watershed_image_filepath=path_to_watershed_images;
                watershed_image_filepath.append(filename_of_the_image);
                watershed_image_filepath.append(".vec");

                //NOW WE CAN START THE TRAINING
                Boundary_training_image boundary_train_image(argv[i],watershed_image_filepath,path_to_training_data_output,argv[argc-1]);
                boundary_train_image.file_to_training_data();
                boundary_train_image.do_training();
                i++;
            }
        }

        /*
        OPTION -boundary-classification-set  (used to create a boundary trainingset)
        */
        else if(command_line_option=="-boundary-classification-set")
        {
            std::cout<<"-boundary-classification-set"<<std::endl;

            //needs 3 pathes ,
                //one filepath to the boundary feature file (if batch-modi this are more than one filepathes
                //one path to the folder where the "boundary training files" are stored
                //one path to the folder where the classification file should be stored

            //filepath to boundary training file == argv[i] (i>=2)
            std::string path_to_boundary_features=argv[argc-2];
            std::string filepath_to_classification_file=argv[argc-1];

            int i=2;
            while(i<argc-2)
            {
                std::cout<<"append: "<<argv[i]<<std::endl;

                create_boundary_classification_set(argv[i],path_to_boundary_features,filepath_to_classification_file);
                i++;
            }
        }

        /*
        OPTION -boundary-classification-gui
        */
        else if(command_line_option=="-boundary-classification-gui")
        {
            std::string path_to_boundary_features=argv[argc-4];
            std::string filepath_to_classification_file=argv[argc-3];
            filepath_to_classification_file.append(argv[argc-2]);

            bool overwrite;
            if (atoi(argv[argc-1])==1) overwrite = true;
            else overwrite = false;

            int i=2;
            while(i<argc-5)
            {
                std::string filepath_to_boundary_labels=argv[argc-5];
                std::string filename=get_filename(argv[i]);
                filepath_to_boundary_labels.append(filename.c_str());
                filepath_to_boundary_labels.append(".dat");

                create_boundary_classification_set(filepath_to_boundary_labels,path_to_boundary_features,filepath_to_classification_file,
                    overwrite);
                i++;
            }
        }

        /*
        OPTION -boundary-rf-train  (used to train the boundary-random-forest with the boundary trainingset)
        */
        else if(command_line_option=="-boundary-rf-train")
        {
            std::cout<<"-boundary-rf-train"<<std::endl;

            std::string filepath_to_classification_set_file=argv[argc-2];
            std::cout<<"Filepath to classification set file: "<<filepath_to_classification_set_file<<std::endl;
            std::string dest_filepath_to_random_forest_file=argv[argc-1];
            train_boundary_rf(filepath_to_classification_set_file,dest_filepath_to_random_forest_file,"parameters.txt");
        }

        /*
        OPTION -boundary-rf-train-gui
        */
        else if(command_line_option=="-boundary-rf-train-gui")
        {
            std::string filepath_to_classification_set_file=argv[argc-4];
            filepath_to_classification_set_file.append(argv[argc-2]);
            std::cout<<"Filepath to classification set file: "<<filepath_to_classification_set_file<<std::endl;
            std::string dest_filepath_to_random_forest_file=argv[argc-3];
            dest_filepath_to_random_forest_file.append(argv[argc-2]);
            train_boundary_rf(filepath_to_classification_set_file,dest_filepath_to_random_forest_file,argv[argc-1]);
        }

        /*
        OPTION -boundary-rf-predict  (used to ...)
        */
        else if(command_line_option=="-boundary-rf-predict")
        {
            std::cout<<"-boundary-rf-predict"<<std::endl;

            std::string path_to_ws_image=argv[argc-3];
            std::string filepath_to_random_forest_file=argv[argc-2];
            std::string path_to_output_folder=argv[argc-1];

            std::string path_to_gm=path_to_output_folder;

            int i=2;
            while(i<argc-3)
            {
                extract_boundary_probabilities(argv[i],path_to_ws_image,filepath_to_random_forest_file,path_to_output_folder,path_to_gm,
                    "parameters.txt");
                i++;
            }

        }

        /*
        OPTION -boundary-complete-param  (used when labeling for at least one image is given to do all steps
                                       until ws-regions use standard path, less arguments
                                       BATCH mode for parameter files)
        */
        else if(command_line_option=="-boundary-complete-param")
        {
            std::cout<<"Option -boundary-complete-param"<<std::endl;

            // some important parameters should not need to be changed in the source code
            // so let's load them from a parameter file
            ParameterFile paramFile;

            int i=3;
            while(i<argc)
            {
                bool default_rf_trained=false;

                std::cout<<"Parameter file "<<i-2<<" of " <<argc-3<<": "<<argv[i]<<std::endl;

                if( !paramFile.load(argv[i]) )
                {
                    std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                    return 0;
                }
                do_boundary_complete(paramFile, argv[2], argv[i], &default_rf_trained);
                i++;
            }
        }

        /*
        OPTION -boundary-complete-image  (used when labeling for at least one image is given to do all steps
                                       until ws-regions use standard path, less arguments
                                       BATCH mode for images)
        */
        else if(command_line_option=="-boundary-complete-image")
        {
            std::cout<<"Option -boundary-complete-image"<<std::endl;

            // some important parameters should not need to be changed in the source code
            // so let's load them from a parameter file
            ParameterFile paramFile;

            if( !paramFile.load("parameters.txt") )//default parameter file for image batch mode
            {
                std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                return 0;
            }

            bool default_rf_trained=false;

            int i=2;
            while(i<argc)
            {
                std::cout<<"Image nr: "<<i-1<<" of " <<argc-2<<": "<<argv[i]<<std::endl;
                do_boundary_complete(paramFile, argv[i], "", &default_rf_trained);
                i++;
            }
        }

        /*
        OPTION -junctions (used to compute the junctions-features of an image with labels, create classification set and train random forest)
        */
        else if(command_line_option=="-junctions")
        {
            std::cout<<"Option -junction"<<std::endl;
            std::string source_ws_image_path=argv[argc-4];
            std::string filepath_to_junction_training_file=argv[argc-3];
            std::string filepath_to_classification_file=argv[argc-2];
            std::string dest_filepath_to_random_forest_file=argv[argc-1];
            junction_features_and_classification_set(argv[argc-5],source_ws_image_path,filepath_to_junction_training_file,
                filepath_to_classification_file);
            train_boundary_rf(filepath_to_classification_file,dest_filepath_to_random_forest_file,"parameters.txt");
        }        

        /*
         * OPTION -bubbles (used to mark bubbles within a picture)
         */
        else if(command_line_option == "-bubbles")
        {
            std::cout << "Option -bubbles" << std::endl;
            std::string source_proc_image_path = argv[argc-1];

            bool firn=false;
            std::cout<<"Are bubbles closed off? Press 'n' for firn or 'y' for ice: ";
            std::string in;
            std::cin >> in;
            if(in=="n") firn=true;

            int i=2;
            while(i<argc-1)
            {
                std::string source_image_path = argv[i];
           
                std::cout << "Source Image: " << source_image_path << std::endl;
                std::cout << "Preprocessed Image: " << source_proc_image_path << std::endl;
                extract_bubbles(source_image_path,source_proc_image_path,firn);
                i++;
            }
        }

        /*
         * OPTION -bubbles-gui (gui version)
         */
        else if(command_line_option == "-bubbles-gui")
        {
            std::cout << "Option -bubbles" << std::endl;
            std::string source_proc_image_path = argv[argc-3];

            bool firn=false;
            if (atoi(argv[argc-1])==1) firn = true;

            int i=2;
            while(i<argc-3)
            {
                std::string source_image_path = argv[i];
           
                std::cout << "Source Image: " << source_image_path << std::endl;
                std::cout << "Preprocessed Image: " << source_proc_image_path << std::endl;
                extract_bubbles(source_image_path,source_proc_image_path,firn,argv[argc-2]);
                i++;
            }
        }

        /*
         * OPTION -elle-export
         */
        else if(command_line_option == "-elle-export")
        {
            std::cout << "Option -elle-export" << std::endl;
            std::string elleDataSetFilepath = argv[2];
            std::string elleExportFilepath = argv[3];
            std::string elleImageFilename = elleDataSetFilepath;
            elleImageFilename = get_filename(elleImageFilename);
            elleImageFilename = elleImageFilename.substr(0, elleImageFilename.length() - 5);
            elleImageFilename.append(".bmp");

            elleExportFilepath.append(elleImageFilename);

            int dim_x = atoi(argv[4]);
            int dim_y = atoi(argv[5]);

            ElleDataSet elleDS(elleDataSetFilepath);
            elleDS.exportToImage(dim_x, dim_y, elleExportFilepath);
        }

        else
        {
             std::cout<<"-WRONG OPTION!!!!"<<std::endl;
        }    
    }

    return 0;
}
