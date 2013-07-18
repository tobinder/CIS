/*! \file compute_watershed_regions.h
 * \brief Calculation of watershed regions.
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

#include <vigra/watersheds.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/localminmax.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/flatmorphology.hxx>
#include <vigra/pixelneighborhood.hxx>

#include <vigra/orientedtensorfilters.hxx>

#include <math.h>
#include <cgp/cgp_hdf5.hxx>
#include "cgp_structure.hxx"
#include "marray.hxx"

#include "CImg.h"
#include <omp.h>

//Used for handling colors
#include <vigra/rgbvalue.hxx>

#include "ParameterFile.hxx"

using namespace cimg_library;
#undef min
#undef max

template <class T> struct EqualWithToleranceFunctor
{
    EqualWithToleranceFunctor(T tolerance)
    : t(tolerance)
    {}

    bool operator()(T l, T r) const
    {
        return vigra::abs(l-r) <= t;
    }

    T t;
};

/*! \fn compute_ws_regions(std::string source_image_filepath,std::string preprocessed_probmap_filepath,std::string dest_image_path,int limit,int equalTolerance)
 * \brief Watershed segmentation with given preprocessed image.
 * \param source_image_filepath Filepath to the source image
 * \param preprocessed_probmap_filepath Filepath to the preprocessed probability map image
 * \param dest_image_path Path to the destination image
 * \param limit Image value threshold
 * \param equalTolerance Equal tolerance functor
 */
//Exports the watershed segmentation to a hdf5 file
void export_ws_hdf5(std::string &dest_path, unsigned int** ws_image, int dim_x, int dim_y)
{
    //EXPORT SEGMENTATION TO A HDF5 file
    dest_path.append(".h5");

    std::cout<<"Exporting watershed image to file: "<<std::endl;
    std::cout<<dest_path<<std::endl;
    
    //creating the file
    hid_t file_save=H5Fcreate(dest_path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    //dataspace for ws_image
    hsize_t dims[3];
    dims[0] = dim_y; 
    dims[1] = dim_x; 
    dims[2] = 1;
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

    //dataset for ws_image    
    hid_t dataset_id = H5Dcreate1(file_save, "/ws_image", H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT);

    //dataspace for one row
    hsize_t row[2];
    row[0] = 1;
    row[1] = dim_x;  
    row[2] = 1;
    hid_t mdataspace_id = H5Screate_simple(3, row, NULL);

    //loop over rows in ws_image
    for(int i= 0; i<dim_y; i++)
    {
        // select file hyperslab 
        hsize_t start[3];// start of hyperslab
        hsize_t count[3];// block count
        
        count[0]  = 1; 
        count[1]  = dim_x;
        count[2]  = 1; 

        start[0]  = i; 
        start[1]  = 0;
        start[2]  = 0;
        
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
        
        unsigned int *row_values = ws_image[i];
        H5Dwrite(dataset_id, H5T_NATIVE_UINT, mdataspace_id, dataspace_id, H5P_DEFAULT, row_values);

        delete row_values;
    }

    {
        //Create group Parameters
        hid_t group_id = H5Gcreate1(file_save, "/Parameters", 0);

        std::cout << "Loading parameter file from: " << ParameterFile::filepath;
        ParameterFile paramFile;
        paramFile.load(ParameterFile::filepath);
        
        //Write the parameters into the group
        Parameter<std::string>::writeParamToHDF5("filepath_default_pixel_rf", "", group_id, &paramFile);
       
        Parameter<float>::writeParamToHDF5("gray_threshold", (float)0.51, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("light_threshold", (float)131, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("feature_threshold", (float)0.01, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("no_boundary_probability", (float)0.951, group_id, &paramFile);

        Parameter<int>::writeParamToHDF5("threshold", 26, group_id, &paramFile);
        Parameter<int>::writeParamToHDF5("scale", 1, group_id, &paramFile);
        Parameter<int>::writeParamToHDF5("equal_tolerance", 1, group_id, &paramFile);
        
        Parameter<int>::writeParamToHDF5("blocksize", 100, group_id, &paramFile);
        Parameter<int>::writeParamToHDF5("grayvalue_difference", 300, group_id, &paramFile);
        
        Parameter<std::string>::writeParamToHDF5("grayvalue", "false", group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("gradient_magnitude1", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("gradient_magnitude2", (float)1.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("gradient_magnitude3", (float)3.0, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor1a", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor1b", (float)0.6, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure1", (float)0.0, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor2a", (float)1.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor2b", (float)0.9, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure2", (float)0.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor3a", (float)2.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor3b", (float)1.2, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure3", (float)0.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_matrix1", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian1", (float)0.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_matrix2", (float)1.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian2", (float)0.0, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_matrix3", (float)2.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian3", (float)0.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_boundary_strength1", (float)0.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_boundary_strength2", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_boundary_strength3", (float)1.5, group_id, &paramFile);
        
        Parameter<std::string>::writeParamToHDF5("grayvalue_as", "false", group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("gradient_magnitude_as1", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("gradient_magnitude_as2", (float)1.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("gradient_magnitude_as3", (float)3.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor_as1a", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor_as1b", (float)0.6, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_as1", (float)1.0, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor_as2a", (float)1.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor_as2b", (float)0.9, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_as2", (float)1.0, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor_as3a", (float)2.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_tensor_as3b", (float)1.2, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_structure_as3", (float)1.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_matrix_as1", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_as1", (float)1.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_matrix_as2", (float)1.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_as2", (float)1.0, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_matrix_as3", (float)3.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("eigenvalues_hessian_as3", (float)1.0, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("boundary_strength_as1", (float)0.5, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("boundary_strength_as2", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("boundary_strength_as3", (float)1.5, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("difference_of_gaussians1a", (float)6.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("difference_of_gaussians1b", (float)0.8, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("difference_of_gaussians2a", (float)3.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("difference_of_gaussians2b", (float)0.4, group_id, &paramFile);
        
        Parameter<float>::writeParamToHDF5("laplace_of_gaussian1", (float)1.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("laplace_of_gaussian2", (float)2.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("laplace_of_gaussian3", (float)3.0, group_id, &paramFile);

        Parameter<float>::writeParamToHDF5("abs_difference_of_gaussians1a", (float)6.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("abs_difference_of_gaussians1b", (float)0.8, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("abs_difference_of_gaussians2a", (float)3.0, group_id, &paramFile);
        Parameter<float>::writeParamToHDF5("abs_difference_of_gaussians2b", (float)0.4, group_id, &paramFile);

        //Close the group
        H5Gclose(group_id);

    }

    H5Sclose(mdataspace_id);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    delete ws_image;

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;
}

template <class InImage, class OutImage, class WSLabels>
void watershedSegmentation_regions(InImage & in, OutImage & out, WSLabels & labels_ws, int equal_tolerance)
{
    int w = in.width();
    int h = in.height();

    //In this image we save the watershed-segmentation
    vigra::IImage labels_img(w,h);
    labels_img = 0;

    vigra::EightNeighborhood::NeighborCode neighbor;

    //TODO can this be done more efficient?
    for(int y=0;y<h;y++)
    {
        for(int x=0;x<w;x++)
        {
       		if (labels_ws(x,y)==255) labels_img(x,y)=1;
        }
    }

    // find the local minima of the gradient magnitude
    // (might be larger than one pixel)
    vigra::extendedLocalMinima(vigra::srcImageRange(in),vigra::destImage(labels_img),1,neighbor,EqualWithToleranceFunctor<unsigned char>(equal_tolerance));

    // the water will start to rise from those pixels
	// exportImage(srcImageRange(labels_img), vigra::ImageExportInfo("localminima.bmp"));

    // label the minima just found
    int max_region_label = vigra::labelImageWithBackground(vigra::srcImageRange(labels_img), vigra::destImage(labels_img), false, 0);

    // create a statistics functor for region growing
    vigra::ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<float> > gradstat(max_region_label);

    // perform region growing, starting from the minima of the gradient magnitude;
    // as the feature (first input) image contains the gradient magnitude,
    // this calculates the catchment basin of each minimum
    vigra::seededRegionGrowing(vigra::srcImageRange(in), vigra::srcImage(labels_img), vigra::destImage(labels_img), gradstat, vigra::CompleteGrow);

    vigra::copyImage(srcImageRange(labels_img), destImage(out));
}

/*! \fn compute_ws_regions(std::string source_image_filepath,std::string preprocessed_probmap_filepath,std::string dest_image_path,int limit,int equalTolerance)
 * \brief Watershed segmentation with given preprocessed image.
 * \param source_image_filepath Filepath to the source image
 * \param preprocessed_probmap_filepath Filepath to the preprocessed probability map image
 * \param dest_image_path Path to the destination image
 * \param limit Image value threshold
 * \param equalTolerance Equal tolerance functor
 */
//Watershed-segmentation with given preprocessed image
void compute_ws_regions(std::string source_image_filepath,std::string preprocessed_probmap_filepath,std::string dest_image_path,int limit,int equalTolerance)
{
    std::string filename_of_image=get_filename(preprocessed_probmap_filepath);
    //std::string filename_of_image=get_filename(source_image_path);
    std::string dest_path=dest_image_path;
    std::string filepath_preprocessed_image=preprocessed_probmap_filepath;
    //std::string filepath_preprocessed_image=preprocessed_image_path;

    dest_path.append(filename_of_image);
    //filepath_preprocessed_image.append(filename_of_image);
    //filepath_preprocessed_image.append(".as.bmp");

    //import the preprocessed image
    vigra::ImageImportInfo info(filepath_preprocessed_image.c_str());
    int dim_x=info.width();
    int dim_y=info.height();

    vigra::FImage preprocessed_image;
    preprocessed_image.resize(dim_x,dim_y);
    importImage(info, destImage(preprocessed_image));

	transformImage(srcImageRange(preprocessed_image), destImage(preprocessed_image), vigra::linearIntensityTransform(10,0));

 	for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
			
       		if(preprocessed_image(x,y)<=limit)
			{
				preprocessed_image(x,y)=0;
			}
        }
    }

    vigra::BasicImage<unsigned int> ws_region_image(dim_x,dim_y);

    filepath_preprocessed_image.resize(filepath_preprocessed_image.size()-3);
    filepath_preprocessed_image.append("0.bmp");    

    //Check if ws-labels for this image exist
    FILE *fp_test;
    fp_test=fopen(filepath_preprocessed_image.c_str(),"r");
    //In this image we have stored ws_label
    vigra::IImage labels_img(dim_x,dim_y);

    if(fp_test!=NULL) //labels for this picture
    {
        vigra::ImageImportInfo info_ws(filepath_preprocessed_image.c_str());
        importImage(info_ws, destImage(labels_img));
    }
    else
    {
        labels_img = 0;
    }

    //do the watershed segmentation
    watershedSegmentation_regions(preprocessed_image,ws_region_image,labels_img,equalTolerance);

    //export the ws region image als std::vector
    std::vector<unsigned int>  labeling_vector_x;
    labeling_vector_x.resize(dim_x);

    std::vector<std::vector<unsigned int > > labeling_vector;
    labeling_vector.resize(dim_y,labeling_vector_x);

    //just for visualisation of the regions
    vigra::FImage regions(dim_x,dim_y);
    int x_old=0;
    int y_old[dim_x];
    for (int x=0;x<dim_x;x++) y_old[x]=0;

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           labeling_vector[y][x]=ws_region_image(x,y);

           if (ws_region_image(x,y)!=x_old || ws_region_image(x,y)!=y_old[x] || x==(dim_x-1) || y==(dim_y-1)) regions(x,y)=0;
           else regions(x,y)=1;
           x_old=ws_region_image(x,y);
           y_old[x]=ws_region_image(x,y);
        }
    }

    //export visualisation
    exportImage(srcImageRange(regions), vigra::ImageExportInfo(dest_path.c_str()));

    //an array to store ws region image
    unsigned int** ws_image=new unsigned int* [dim_y];
    for (int i=0; i<dim_y; i++) ws_image[i]=new unsigned int[dim_x];

    for(int y=0;y<(int)labeling_vector.size();y++)
    {
        for(int x=0;x<(int)labeling_vector[0].size();x++)
        {
            ws_image[y][x]=labeling_vector[y][x];
        }
    }

    //EXPORT SEGMENTATION TO A HDF5 file
    /*dest_path.append(".h5");

    std::cout<<"Exporting watershed image to file: "<<std::endl;
    std::cout<<dest_path<<std::endl;
    
    //creating the file
    hid_t file_save=H5Fcreate(dest_path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    //dataspace for ws_image
    hsize_t dims[3];
    dims[0] = dim_y; 
    dims[1] = dim_x; 
    dims[2] = 1;
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

    //dataset for ws_image    
    hid_t dataset_id = H5Dcreate1(file_save, "/ws_image", H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT);

    //dataspace for one row
    hsize_t row[2];
    row[0] = 1;
    row[1] = dim_x;  
    row[2] = 1;
    hid_t mdataspace_id = H5Screate_simple(3, row, NULL);

    //loop over rows in ws_image
    for(int i= 0; i<dim_y; i++)
    {
        // select file hyperslab 
        hsize_t start[3];// start of hyperslab
        hsize_t count[3];// block count
        
        count[0]  = 1; 
        count[1]  = dim_x;
        count[2]  = 1; 

        start[0]  = i; 
        start[1]  = 0;
        start[2]  = 0;
        
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
        
        unsigned int *row_values = ws_image[i];
        H5Dwrite(dataset_id, H5T_NATIVE_UINT, mdataspace_id, dataspace_id, H5P_DEFAULT, row_values);

        delete row_values;
    }

    H5Sclose(mdataspace_id);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    delete ws_image;

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;*/
    export_ws_hdf5(dest_path, ws_image, dim_x, dim_y);
}

/*! \fn compute_watershed_regions(std::string source_image_path,std::string dest_image_path,int limit,double scale,int equalTolerance, int size, int rank)
 * \brief Watershed segmentation with calculation of a preprocessed image.
 * \param source_image_filepath Filepath to the source image
 * \param dest_image_path Path to the destination image
 * \param limit Image value threshold
 * \param scale Scale
 * \param equalTolerance Equal tolerance functor
 * \param size Size
 * \param rank Rank
 */
//Watershed-segmentation with calculation of a preprocessed image
void compute_watershed_regions(std::string source_image_path,std::string dest_image_path,int limit,double scale,int equalTolerance, int size, int rank)
{
    std::cout<<"Compute Watershed Segmentation on image: "<<get_filename(source_image_path)<<std::endl;

    vigra::ImageImportInfo info(source_image_path.c_str());
    int dim_x = info.width();
    int dim_y = info.height();

    vigra::BImage in(dim_x, dim_y);
    importImage(info, destImage(in));

    //Adjust image filename for "eed_ms.bmp" suffix
    if(source_image_path.size()>11)
        if(source_image_path.compare(source_image_path.length()-11, 11, ".eed_ms.bmp") == 0)
        {
            source_image_path.erase(source_image_path.length()-11, 11);
        }

    std::cout<<"Start Preprocessing..."<<std::endl;
    int blocksize = std::min(2000,std::min(dim_x,dim_y));//has to be larger than all structures in image!
    std::cout<<"Image size: ("<<dim_x<<","<<dim_y<<")"<<std::endl;
    std::cout<<"Block size: ("<<blocksize<<","<<blocksize<<")"<<std::endl;

    vigra::BImage preprocessing_float_image(dim_x, dim_y);
    vigra::BImage preprocessing_float_image_s(dim_x, dim_y);

    #pragma omp parallel
    {
        #pragma omp for
        for(int y=0;y<(int)(dim_y/blocksize);y++)
        {
            for(int x=0;x<(int)(dim_x/blocksize);x++)
	        {
                vigra::BImage block;
                if ((x+2)*blocksize<=dim_x && (y+2)*blocksize<=dim_y)
                {
                    block.resize(blocksize, blocksize);
                }
                else if ((x+2)*blocksize<=dim_x)
                {
                    block.resize(blocksize, dim_y-(y*blocksize));
                }
                else if ((y+2)*blocksize<=dim_y)
                {
                    block.resize(dim_x-(x*blocksize), blocksize);
                }
                else
                {
                    block.resize(dim_x-(x*blocksize), dim_y-(y*blocksize));
                }

                //fill the block
                std::cout<<"Block "<<y*(int)(dim_x/blocksize)+x+1<<"/"<<(int)(dim_x/blocksize)*(int)(dim_y/blocksize)<<", Offset: ("<<
                    x*blocksize<<","<<y*blocksize<<"), Block size: ("<<block.width()<<","<<block.height()<<")"<<std::endl;
                for(int y_block=0; y_block<block.height(); y_block++)
                {
                    for(int x_block=0; x_block<block.width(); x_block++)
                    {
                        block(x_block,y_block)=in(x*blocksize+x_block,y*blocksize+y_block);
                    }
                }

                vigra::BImage reversed(block.width(), block.height());
                vigra::BImage smoothed(block.width(), block.height());

	            CImg<float>  src(block.width(),block.height());

	            transformImage(srcImageRange(block), destImage(reversed), vigra::linearIntensityTransform(-1, -255));

                for(int y_block=0;y_block<block.height();y_block++)
                {
                    for(int x_block=0;x_block<block.width();x_block++)
                    {
                   		src(x_block,y_block)=(float)reversed(x_block,y_block);
                    }
                }

                CImgList<double> images(src,src);

                for (unsigned int iter = 0; iter<200; ++iter)
                {
                    // Compute PDE velocity field.
                    CImg_3x3(I,double);
                    CImg<> veloc(src);
                    cimg_forV(src,k) cimg_for3x3(images[1],x,y,0,k,I)
                    {
                        const double
                        ix = (Inc - Ipc)/2,
                        iy = (Icn - Icp)/2,
                        ng = (double)std::sqrt(1e-10f + ix*ix + iy*iy),
                        ixx = Inc + Ipc - 2*Icc,
                        iyy = Icn + Icp - 2*Icc,
                        ixy = 0.25f*(Inn + Ipp - Ipn - Inp),
                        iee = (ix*ix*iyy + iy*iy*ixx - 2*ix*iy*ixy)/(ng*ng);
                        veloc(x,y,k) = iee/(0.1f+ng);
                    }

                    // Find adaptive time step and update current image.
                    double m = 0, M = veloc.maxmin(m);
                    veloc*=40.0f/cimg::max(cimg::abs(m),cimg::abs(M));
                    images[1]+=veloc;
                }
	
                images[1].normalize(0,255);

                //write back
                for(int y_block=0;y_block<block.height();y_block++)
                {
                    for(int x_block=0;x_block<block.width();x_block++)
                    {
			            preprocessing_float_image(x*blocksize+x_block,y*blocksize+y_block)=images[1](x_block,y_block);
                    }
                }
            }
        }
    }//end of parallelisation

    vigra::gaussianSmoothing(srcImageRange(preprocessing_float_image),destImage(preprocessing_float_image_s),scale);
 	for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
       		if(preprocessing_float_image_s(x,y)<=limit)
            {
	            preprocessing_float_image_s(x,y)=0;
            }
        }
    }

	//exportImage(srcImageRange(preprocessing_float_image_s), vigra::ImageExportInfo("prob_map_preprocessed.bmp"));
    std::cout<<"...done"<<std::endl;

    vigra::BasicImage<unsigned int> ws_region_image(dim_x,dim_y);

    std::string filename_of_image=get_filename(source_image_path);
    source_image_path.resize(source_image_path.size()-3);
    source_image_path.append("0.bmp");    

    //Check if ws-labels for this image exist
    FILE *fp_test;
    fp_test=fopen(source_image_path.c_str(),"r");
    //In this image we have stored ws_label
    vigra::IImage labels_img(dim_x,dim_y);

    if(fp_test!=NULL) //labels for this picture
    {
        std::cout<<"Watershed growing with additional given areas..."<<std::endl;
        vigra::ImageImportInfo info_ws(source_image_path.c_str());
        importImage(info_ws, destImage(labels_img));
    }
    else
    {
        std::cout<<"Watershed growing..."<<std::endl;
        labels_img = 0;
    }
 
    watershedSegmentation_regions(preprocessing_float_image_s,ws_region_image,labels_img,equalTolerance);
    std::string dest_path=dest_image_path;
    dest_path.append(filename_of_image);

    std::cout<<"...done"<<std::endl;

    //int nr_areas=0;

    //export the ws region image als std::vector
    std::vector<unsigned int>  labeling_vector_x;
    labeling_vector_x.resize(dim_x);

    std::vector<std::vector<unsigned int > > labeling_vector;
    labeling_vector.resize(dim_y,labeling_vector_x);

    //visualisation of the regions
    vigra::FImage regions(dim_x,dim_y);
    int x_old=0;
    int y_old[dim_x];

    for (int x=0;x<dim_x;x++) y_old[x]=0;
 
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           labeling_vector[y][x]=ws_region_image(x,y);

           if (ws_region_image(x,y)!=x_old || ws_region_image(x,y)!=y_old[x] || x==(dim_x-1) || y==(dim_y-1)) regions(x,y)=0;
           else regions(x,y)=1;
           x_old=ws_region_image(x,y);
           y_old[x]=ws_region_image(x,y);

           //if (ws_region_image(x,y)>nr_areas) nr_areas=ws_region_image(x,y);
        }
    }

    //std::cout<<"Areas: "<<nr_areas<<std::endl;

    //export the nr of areas
    //std::ofstream areas_file(dest_path.append(".area.dat").c_str(), std::ios::out | std::ios::app);
    //areas_file <<nr_areas<< "\n";
    //areas_file.close();
    //dest_path.resize(dest_path.size()-9);

    //export visualisation
    exportImage(srcImageRange(regions), vigra::ImageExportInfo(dest_path.c_str()));

    //an array to store ws region image
    unsigned int** ws_image=new unsigned int* [dim_y];
    for (int i=0; i<dim_y; i++) ws_image[i]=new unsigned int[dim_x];

    for(int y=0;y<(int)labeling_vector.size();y++)
    {
        for(int x=0;x<(int)labeling_vector[0].size();x++)
        {
            ws_image[y][x]=labeling_vector[y][x];
        }
    }

    //EXPORT SEGMENTATION TO A HDF5 file
    /*dest_path.append(".h5");

    std::cout<<"Exporting watershed image to file: "<<std::endl;
    std::cout<<dest_path<<std::endl;
    
    //creating the file
    hid_t file_save=H5Fcreate(dest_path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    //dataspace for ws_image
    hsize_t dims[3];
    dims[0] = dim_y; 
    dims[1] = dim_x; 
    dims[2] = 1;
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

    //dataset for ws_image    
    hid_t dataset_id = H5Dcreate1(file_save, "/ws_image", H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT);

    //dataspace for one row
    hsize_t row[3];
    row[0] = 1;
    row[1] = dim_x;  
    row[2] = 1;
    hid_t mdataspace_id = H5Screate_simple(3, row, NULL);

    //loop over rows in ws_image
    for(int i= 0; i<dim_y; i++)
    {
        // select file hyperslab 
        hsize_t start[3];// start of hyperslab
        hsize_t count[3];// block count
        
        count[0]  = 1; 
        count[1]  = dim_x;
        count[2]  = 1; 

        start[0]  = i; 
        start[1]  = 0;
        start[2]  = 0;
        
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
        
        unsigned int *row_values = ws_image[i];
        H5Dwrite(dataset_id, H5T_NATIVE_UINT, mdataspace_id, dataspace_id, H5P_DEFAULT, row_values);

        delete row_values;
    }

    H5Sclose(mdataspace_id);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    delete ws_image;

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;*/

    export_ws_hdf5(dest_path, ws_image, dim_x, dim_y);

    std::string dest_cgp_path=dest_path;
    dest_cgp_path.resize(dest_cgp_path.size()-3);

    compute_and_save_cgp_data_structure(dest_cgp_path,size,rank);
}

/*! \fn compute_watershed_regions_binary(std::string source_image_path,std::string dest_image_path,int equalTolerance, int size, int rank)
 * \brief Watershed segmentation from a manually generated binary map.
 * \param source_image_filepath Filepath to the source image
 * \param dest_image_path Path to the destination image
 * \param equalTolerance Equal tolerance functor
 * \param size Size
 * \param rank Rank
 */
//Watershed-segmentation from manually generated binary map
void compute_watershed_regions_binary(std::string source_image_path,std::string dest_image_path,int equalTolerance, int size, int rank, vigra::IImage* converted_source_image = NULL, bool converted_source_image_available = false)
{
    std::cout<<"Compute Watershed Segmentation from manually generated image: "<<get_filename(source_image_path)<<std::endl;
    vigra::ImageImportInfo info(source_image_path.c_str());
    int dim_x = info.width();
    int dim_y = info.height();

    vigra::IImage labels_img(dim_x, dim_y);
    if(converted_source_image_available == false)
    {
        importImage(info, destImage(labels_img));
    }
    else
    {
        labels_img = *converted_source_image;
        delete converted_source_image;
    }

    vigra::BImage in(dim_x, dim_y);
    transformImage(srcImageRange(labels_img), destImage(in), vigra::linearIntensityTransform(-1, -255));

    vigra::BasicImage<unsigned int> ws_region_image(dim_x,dim_y);

    std::string filename_of_image=get_filename(source_image_path);
 
    watershedSegmentation_regions(in,ws_region_image,labels_img,equalTolerance);
    std::string dest_path=dest_image_path;
    dest_path.append(filename_of_image);

    std::cout<<"...done"<<std::endl;

    //export the ws region image als std::vector
    std::vector<unsigned int>  labeling_vector_x;
    labeling_vector_x.resize(dim_x);

    std::vector<std::vector<unsigned int > > labeling_vector;
    labeling_vector.resize(dim_y,labeling_vector_x);

    //visualisation of the regions
    vigra::FImage regions(dim_x,dim_y);
    int x_old=0;
    int y_old[dim_x];

    for (int x=0;x<dim_x;x++) y_old[x]=0;
 
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           labeling_vector[y][x]=ws_region_image(x,y);

           if (ws_region_image(x,y)!=x_old || ws_region_image(x,y)!=y_old[x] || x==(dim_x-1) || y==(dim_y-1)) regions(x,y)=0;
           else regions(x,y)=1;
           x_old=ws_region_image(x,y);
           y_old[x]=ws_region_image(x,y);
        }
    }

    //export visualisation
    exportImage(srcImageRange(regions), vigra::ImageExportInfo(dest_path.c_str()));

    //an array to store ws region image
    unsigned int** ws_image=new unsigned int* [dim_y];
    for (int i=0; i<dim_y; i++) ws_image[i]=new unsigned int[dim_x];

    for(int y=0;y<(int)labeling_vector.size();y++)
    {
        for(int x=0;x<(int)labeling_vector[0].size();x++)
        {
            ws_image[y][x]=labeling_vector[y][x];
        }
    }

    //EXPORT SEGMENTATION TO A HDF5 file
    /*dest_path.append(".h5");

    std::cout<<"Exporting watershed image to file: "<<std::endl;
    std::cout<<dest_path<<std::endl;
    
    //creating the file
    hid_t file_save=H5Fcreate(dest_path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    //dataspace for ws_image
    hsize_t dims[3];
    dims[0] = dim_y; 
    dims[1] = dim_x; 
    dims[2] = 1;
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

    //dataset for ws_image    
    hid_t dataset_id = H5Dcreate1(file_save, "/ws_image", H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT);

    //dataspace for one row
    hsize_t row[3];
    row[0] = 1;
    row[1] = dim_x;  
    row[2] = 1;
    hid_t mdataspace_id = H5Screate_simple(3, row, NULL);

    //loop over rows in ws_image
    for(int i= 0; i<dim_y; i++)
    {
        // select file hyperslab 
        hsize_t start[3];// start of hyperslab
        hsize_t count[3];// block count
        
        count[0]  = 1; 
        count[1]  = dim_x;
        count[2]  = 1; 

        start[0]  = i; 
        start[1]  = 0;
        start[2]  = 0;
        
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
        
        unsigned int *row_values = ws_image[i];
        H5Dwrite(dataset_id, H5T_NATIVE_UINT, mdataspace_id, dataspace_id, H5P_DEFAULT, row_values);

        delete row_values;
    }

    H5Sclose(mdataspace_id);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    delete ws_image;

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;*/

    export_ws_hdf5(dest_path, ws_image, dim_x, dim_y);

    std::string dest_cgp_path=dest_path;
    dest_cgp_path.resize(dest_cgp_path.size()-3);

    compute_and_save_cgp_data_structure(dest_cgp_path,size,rank);
}

void compute_watershed_regions_binary_color(std::string source_image_path, std::string dest_image_path,int equalTolerance, int size, int rank)
{
    vigra::ImageImportInfo info(source_image_path.c_str());
    int dim_x = info.width();
    int dim_y = info.height();

    vigra::BRGBImage source_image(dim_x, dim_y);
    importImage(info, destImage(source_image));

    //Might be implemented later; check if the image has a border drawn with the segmentation color. If it has one, use the first labeling image, if not, use the second one. Problem: images with boundaries that do not exceed the image border.
    bool hasBorder = true;

    //Image containing the grayscale version of the segmentation
    vigra::IImage * gScaleImage = new vigra::IImage(dim_x, dim_y);
    //vigra::IImage gScaleImage(dim_x, dim_y);

    int totalPixelNumber = dim_x*dim_y;

    //Vector containing all available colors. Must not contain more than two
    std::vector<vigra::RGBValue<int, 0, 1, 2> > source_image_colors;
    std::vector<int> source_image_colors_pixels;
    //source_image_colors_pixels.push_back(0);
    //A rather crafty and resource intensive approach. Can be simplified if we assume that the image is garuanteed to contain exactly two colors.
    bool doNotInsert = false;

    for(int x = 0; x < dim_x; x++)
    {
        for(int y = 0; y < dim_y; y++)
        {
            //Make the grayscale image white
            (*gScaleImage)(x,y) = 255;
            //Get the color value for the pixel and check if it is already in the color vector. Count the instances of this color
            doNotInsert = false;
            vigra::RGBValue<int, 0, 1, 2> xyColor = source_image(x,y);
            for(int i = 0; i < source_image_colors.size(); i++)
            {
                if(xyColor.red() == source_image_colors[i].red() && xyColor.green() == source_image_colors[i].green() && xyColor.blue() == source_image_colors[i].blue())
                {
                    source_image_colors_pixels[i]++;
                    doNotInsert = true;
                    break;
                }
            }
            if(doNotInsert == false)
            {
                source_image_colors.push_back(xyColor);
                source_image_colors_pixels.push_back(1);
            }
        }
    }

    //Abort if the image does not contain two colors
    if(source_image_colors.size() != 2) {
        std::cout << "Error: image does not contain two colors!" << std::endl;
        return;
    }

    std::vector<float> source_image_colors_perc;
    //Print out the image's colors and their percentage within the image
    std::cout << "Image colors (R, G, B):" << "\n";
    for(int i = 0; i < source_image_colors.size(); i++)
    {
        float iPerc = (float)source_image_colors_pixels[i]/(float)totalPixelNumber;
        source_image_colors_perc.push_back(iPerc);
        std::cout << "(" << source_image_colors[i].red() << ", " << source_image_colors[i].green() << ", "<< source_image_colors[i].blue() << "), Percentage = " << iPerc << "\n";
    }
    std::cout << std::endl;
    
    vigra::RGBValue<int, 0, 1, 2> segmentation_color;

    //Take the color with the lowest percentage. This *should* be the color used for the segmentation
    if(source_image_colors_perc[0] <= source_image_colors_perc[1])
    {
        segmentation_color = source_image_colors[0];
    }
    else
    {
        segmentation_color = source_image_colors[1];
    }
    std::cout << "Using (" << segmentation_color.red() << ", " << segmentation_color.green() << ", "<< segmentation_color.blue() << ") as segmentation color" << std::endl;

    //Paint the grayscale image
    for(int x = 0; x < dim_x; x++)
    {
        for(int y = 0; y < dim_y; y++)
        {
            if(source_image(x,y) == segmentation_color)
            {
                (*gScaleImage)(x,y) = 0;
            }
        }
    }

    //Find out which color is marking the grains and which is marking the boundaries. Grains will be disconnected, so we will use this fact as a criterion.
    //Use Vigra's label image feature to label disjoint areas within the image.
    vigra::IImage label_image0(dim_x,dim_y);
    vigra::labelImageWithBackground(vigra::srcImageRange(source_image), vigra::destImage(label_image0), true, source_image_colors[0]);

    vigra::IImage label_image1(dim_x,dim_y);
    vigra::labelImageWithBackground(vigra::srcImageRange(source_image), vigra::destImage(label_image1), true, source_image_colors[1]);

    /*vigra::IImage label_image2(dim_x,dim_y);
    vigra::labelImageWithBackground(vigra::srcImageRange(source_image), vigra::destImage(label_image2), false, source_image_colors[0]);

    vigra::IImage label_image3(dim_x,dim_y);
    vigra::labelImageWithBackground(vigra::srcImageRange(source_image), vigra::destImage(label_image3), false, source_image_colors[1]);*/

    //export the label image
    /*exportImage(srcImageRange(label_image0), vigra::ImageExportInfo("/home/akuehlwe/test_data/pixel-classification/fft_label0.bmp"));
    exportImage(srcImageRange(label_image1), vigra::ImageExportInfo("/home/akuehlwe/test_data/pixel-classification/fft_label1.bmp"));*/
    /*exportImage(srcImageRange(label_image2), vigra::ImageExportInfo("/home/akuehlwe/test_data/pixel-classification/fft_label2.bmp"));
    exportImage(srcImageRange(label_image3), vigra::ImageExportInfo("/home/akuehlwe/test_data/pixel-classification/fft_label3.bmp"));*/
    std::string source_image_path_suffix = source_image_path;
    source_image_path_suffix.substr(0, source_image_path_suffix.length() - 4);
    source_image_path_suffix.append("_color.bmp");
    exportImage(srcImageRange(source_image), vigra::ImageExportInfo(source_image_path_suffix.c_str()));
    exportImage(srcImageRange(*gScaleImage), vigra::ImageExportInfo(source_image_path.c_str()));

    compute_watershed_regions_binary(source_image_path, dest_image_path, equalTolerance, size, rank, gScaleImage, true);
}
