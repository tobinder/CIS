/*! \file feature_extractor_gray.h
 * \brief Feature extraction for grayscale images.
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
#include <sstream>
#include <cstring>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vigra/tensorutilities.hxx>
#include <vigra/boundarytensor.hxx>

class FeatureExtractorGray
{
private:

    std::string filename; /**< filename of the image to extract features on*/
    std::string an_filename;
    std::string path_to_save_feature_file;

    vigra::FImage image; /**< file/image to extract feature on in luv color space*/
    vigra::FImage as_image;
    vigra::FImage gradient_magnitude_image; /**< file/image to extract feature on in luv color space*/
    vigra::FImage temp_image; /**< file/image to  work on,can be overwritten */
    vigra::FImage temp_image2; /**< file/image to  work on,can be overwritten */
    vigra::FImage empty_image;

    std::vector<vigra::FImage> vector_of_images; //ALL FEATURE IMAGES ARE SAFED INTO THIS VECTor

    bool feature_to_image; /**< if feature_to_image is true the images of the extracted features are saved to the filepath*/
    int dim_x; /**<X-Dimension of the image to extract features on*/
    int dim_y; /**<Y-Dimension of the image to extract features on*/

    void export_images(std::string name,double scale1,int channel);
    void export_images(std::string name,double scale1,double scale2,int channel);
    void export_images(std::string name,double scale1,double scale2,double scale3,int channel);

    //FEATURE FUNKTIONS --implemented
    void laplace_of_gaussian(double scale);

    void eigenvalues_structure_tensor(double inner_scale,double outer_scale,int choice);
    void eigenvalues_structure_tensor_as(double inner_scale,double outer_scale,int choice);

    void eigenvalues_hessian_matrix(double scale,int choice);
    void eigenvalues_hessian_matrix_as(double scale,int choice);

    void difference_of_gaussians(double scale1,double scale2);
    void abs_difference_of_gaussians(double scale1,double scale2);

    void gradient_magnitude(double scale);
    void gradient_magnitude_as(double scale);

    void standart_deviation_image(int size); //STANDART ABWEICHUNG IN EINER SIZExSIZE NACHBARSCHAFT //GAUß SMOOTHING VORHER
    void canny_edge_detector(double scale1,double scale2,double scale3);

    void boundary_strength(double scale);
    void boundary_strength_as(double scale);

    vigra::FImage decision_feature();

public:

    FeatureExtractorGray(std::string fn,std::string as_fn,std::string dest_path);
    void get_gradient_magnitude(double scale);
    int extract(bool to_file,std::string param_file_name);
    void extract_training_data(std::string path_labels,std::string filepath_classification_set,int& nr_of_pixel_features,bool overwrite);

};

FeatureExtractorGray::FeatureExtractorGray(std::string fn,std::string as_fn,std::string dest_path)
{
    //SET THE FILENAME
    filename=fn;
    an_filename=as_fn;
    path_to_save_feature_file=dest_path;
    std::cout<<"extract pixel-features on: "<<filename<<std::endl;

    //ITS USED TO SWTCH PICTURES ON OR OFF
    feature_to_image=false;

    //NOW LOAD IMAGES
    vigra::ImageImportInfo info(filename.c_str());
    vigra::ImageImportInfo as_info(an_filename.c_str());

    dim_x=info.width();
    dim_y=info.height();

    //RESIZE IMAGE AND GRADIENT IMAGE TO THE SIZE
    as_image.resize(dim_x,dim_y);
    image.resize(dim_x,dim_y);
    temp_image.resize(dim_x,dim_y);
    temp_image2.resize(dim_x,dim_y);
    empty_image.resize(dim_x,dim_y);
    gradient_magnitude_image.resize(dim_x,dim_y);

    importImage(info, destImage(image));

    //temp image
    vigra::FImage as_temp(as_info.width(),as_info.height());
    importImage(as_info,destImage(as_temp));
    resizeImageSplineInterpolation(srcImageRange(as_temp),destImageRange(as_image));
}


void FeatureExtractorGray::get_gradient_magnitude(double scale)
{
    gaussianGradientMagnitude(srcImageRange(image), destImage(gradient_magnitude_image), scale);
}

void FeatureExtractorGray::export_images(std::string name,double scale1 ,int channel)
{
    std::stringstream ss0;
    ss0<<filename<<"_"<<name<<"_"<<scale1<<"B(0)_"<<".png";
    if(channel==1) exportImage(srcImageRange(temp_image) , vigra::ImageExportInfo(  ss0.str().c_str() ));
    if(channel==2) exportImage(srcImageRange(temp_image2) , vigra::ImageExportInfo(  ss0.str().c_str() ));
}

void FeatureExtractorGray::export_images(std::string name,double scale1,double scale2,int channel)
{
    std::stringstream ss0;
    ss0<<filename<<"_"<<name<<"_"<<scale1<<"_"<<scale2<<"B(0)_"<<".png";
    if(channel==1) exportImage(srcImageRange(temp_image) , vigra::ImageExportInfo(  ss0.str().c_str() ));
    if(channel==2) exportImage(srcImageRange(temp_image2) , vigra::ImageExportInfo(  ss0.str().c_str() ));
}

void FeatureExtractorGray::export_images(std::string name,double scale1,double scale2,double scale3,int channel)
{
    std::stringstream ss0;
    ss0<<filename<<"_"<<name<<"_"<<scale1<<"_"<<scale2<<"_"<<scale3<<"_"<<"B(0)_"<<".png";
    if(channel==1) exportImage(srcImageRange(temp_image) , vigra::ImageExportInfo(  ss0.str().c_str() ));
    if(channel==2) exportImage(srcImageRange(temp_image2) , vigra::ImageExportInfo(  ss0.str().c_str() ));
}

void FeatureExtractorGray::laplace_of_gaussian(double scale)
{
    //RUN laplaceOFGaussian
    vigra::laplacianOfGaussian(srcImageRange(image), destImage(temp_image), scale);

    vector_of_images.push_back(temp_image);

    if(feature_to_image==true)
    {
        export_images("LOG",scale,1);
    }
}

void FeatureExtractorGray::eigenvalues_structure_tensor(double inner_scale,double outer_scale,int choice)
{
    //SCORE IMAGES
    vigra::FImage stxx(dim_x,dim_y), stxy(dim_x,dim_y), styy(dim_x,dim_y);
    structureTensor(srcImageRange(image),destImage(stxx), destImage(stxy), destImage(styy), inner_scale,outer_scale);

    #pragma omp parallel
    {
        //visit all pixels,and compute the eigenvalues
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
	        {
                float lambda1;
                float lambda2;
                float d=((stxx(x,y)-styy(x,y))/2);
                lambda1= (  ((stxx(x,y)+styy(x,y))/2) + sqrt(d*d+ stxy(x,y)*stxy(x,y))  ) ;
                lambda2= (  ((stxx(x,y)+styy(x,y))/2) - sqrt(d*d+ stxy(x,y)*stxy(x,y))  ) ;

                temp_image(x,y)=lambda1;
                temp_image2(x,y)=lambda2;
            }
        }
    }//end of parallelisation

	if (choice!=2) vector_of_images.push_back(temp_image);
	if (choice!=1) vector_of_images.push_back(temp_image2);

    if(feature_to_image==true)
    {
        export_images("EW_struct_tensor lambda1",inner_scale,outer_scale,1);
        export_images("EW_struct_tensor lambda2",inner_scale,outer_scale,2);
    }
}

void FeatureExtractorGray::eigenvalues_structure_tensor_as(double inner_scale,double outer_scale,int choice)
{
    //SCORE IMAGES
    vigra::FImage stxx(dim_x,dim_y), stxy(dim_x,dim_y), styy(dim_x,dim_y);
    structureTensor(srcImageRange(as_image),destImage(stxx), destImage(stxy), destImage(styy), inner_scale,outer_scale);

    #pragma omp parallel
    {
        //visit all pixels,and compute the eigenvalues
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
	        {
                float lambda1;
                float lambda2;
                float d=((stxx(x,y)-styy(x,y))/2.0f);
                lambda1= (  ((stxx(x,y)+styy(x,y))/2.0f) + sqrt(d*d + stxy(x,y)*stxy(x,y))  ) ;
                lambda2= (  ((stxx(x,y)+styy(x,y))/2.0f) - sqrt(d*d + stxy(x,y)*stxy(x,y))  ) ;

                temp_image(x,y)=lambda1;
                temp_image2(x,y)=lambda2;
            }
        }
    }//end of parallelisation

	if (choice!=2) vector_of_images.push_back(temp_image);
	if (choice!=1) vector_of_images.push_back(temp_image2);

    if(feature_to_image==true)
    {
        export_images("AS EW_struct_tensor lambda1",inner_scale,outer_scale,1);
        export_images("AS EW_struct_tensor lambda2",inner_scale,outer_scale,2);
    }
}

void FeatureExtractorGray::eigenvalues_hessian_matrix(double scale,int choice)
{
    //SCORE IMAGES
	vigra::FImage stxx(dim_x,dim_y), stxy(dim_x,dim_y), styy(dim_x,dim_y);
	hessianMatrixOfGaussian(srcImageRange(image),destImage(stxx), destImage(stxy), destImage(styy), scale);

    #pragma omp parallel
    {
        //visit all pixels
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
	        {
                float lambda1;
                float lambda2;
                float d=((stxx(x,y)-styy(x,y))/2);
                lambda1= (  ((stxx(x,y)+styy(x,y))/2) + sqrt(d*d+ stxy(x,y)*stxy(x,y))  ) ;
                lambda2= (  ((stxx(x,y)+styy(x,y))/2) - sqrt(d*d+ stxy(x,y)*stxy(x,y))  ) ;

                temp_image(x,y)=lambda1;
                temp_image2(x,y)=lambda2;
            }
        }
    }//end of parallelisation

    if (choice!=2) vector_of_images.push_back(temp_image);
    if (choice!=1) vector_of_images.push_back(temp_image2);

    if(feature_to_image==true)
    {
        export_images("EW_hessian_matrix_l1",scale,1);
        export_images("EW_hessian_matrix_l2",scale,2);
    }
}

void FeatureExtractorGray::eigenvalues_hessian_matrix_as(double scale, int choice)
{
    //SCORE IMAGES
    vigra::FImage stxx(dim_x,dim_y), stxy(dim_x,dim_y), styy(dim_x,dim_y);
    hessianMatrixOfGaussian(srcImageRange(as_image),destImage(stxx), destImage(stxy), destImage(styy), scale);

    #pragma omp parallel
    {
        //visit all pixels
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
	        {
                float lambda1;
                float lambda2;
                float d=((stxx(x,y)-styy(x,y))/2.0f);
                lambda1= (  ((stxx(x,y)+styy(x,y))/2.0f) + sqrt(d*d+ stxy(x,y)*stxy(x,y))  ) ;
                lambda2= (  ((stxx(x,y)+styy(x,y))/2.0f) - sqrt(d*d+ stxy(x,y)*stxy(x,y))  ) ;

                if(lambda1>=lambda2)
                {
                temp_image(x,y)=lambda1;
                temp_image2(x,y)=lambda2;
                }
                else
                {
                temp_image(x,y)=lambda2;
                temp_image2(x,y)=lambda1;
                }
            }
        }
    }//end of parallelisation
 
    if (choice!=2) vector_of_images.push_back(temp_image);
    if (choice!=1) vector_of_images.push_back(temp_image2);

    if(feature_to_image==true)
    {
        export_images("AS_EW_hessian_matrix_l1",scale,1);
        export_images("AS_EW_hessian_matrix_l2",scale,2);
    }
}

void FeatureExtractorGray::difference_of_gaussians(double scale1,double scale2)
{
    vigra_precondition(scale1>scale2, "scale1 musst be bigger than scale 2");

    //SCORE IMAGES
    vigra::FImage a(dim_x,dim_y), b(dim_x,dim_y) ;

    gaussianSmoothing(srcImageRange(image), destImage(a), scale1);
    gaussianSmoothing(srcImageRange(image), destImage(b), scale2);

    #pragma omp parallel
    {
        //visit all pixels
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
	        {
                temp_image(x,y)=a(x,y)-b(x,y);
            }
        }
    }//end of parallelisation

    vector_of_images.push_back(temp_image);

    if(feature_to_image==true)
    {
        export_images("DOG",scale1,scale2,1);
    }
}

void FeatureExtractorGray::abs_difference_of_gaussians(double scale1,double scale2)
{
    vigra_precondition(scale1>scale2, "scale1 musst be bigger than scale 2");

    //SCORE IMAGES
    vigra::FImage a(dim_x,dim_y), b(dim_x,dim_y) ;

    gaussianSmoothing(srcImageRange(image), destImage(a), scale1);
    gaussianSmoothing(srcImageRange(image), destImage(b), scale2);

    #pragma omp parallel
    {
        //visit all pixels
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
	        {
                temp_image(x,y)=-abs(a(x,y)-b(x,y));
            }
        }
    }//end of parallelisation

    vector_of_images.push_back(temp_image);

    if(feature_to_image==true)
    {
        export_images("ABS(DOG)",scale1,scale2,1);
    }

}

void FeatureExtractorGray::gradient_magnitude(double scale)
{
    //gradient_magnitude at scale scale..
    gaussianGradientMagnitude(srcImageRange(image), destImage(temp_image), scale);

    vector_of_images.push_back(temp_image);

    if(feature_to_image==true)
    {
        export_images("GRDM ",scale,1);
    }
}

void FeatureExtractorGray::gradient_magnitude_as(double scale)
{
    //gradient_magnitude at scale scale..
    gaussianGradientMagnitude(srcImageRange(as_image), destImage(temp_image), scale);

    vector_of_images.push_back(temp_image);

    if(feature_to_image==true)
    {
        export_images("AS GRDM ",scale,1);
    }

}

void FeatureExtractorGray::standart_deviation_image(int size) //SIZE MUSS UNGERADE SEIN!!!
{
    #pragma omp parallel
    {
        //COMPUTE MEAN IMAGE
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
		    {
			    float sum0=0;
			    int y_offset=-1*((size-1)/2);
			    while(y_offset<=(size-1)/2)
			    {
				    int x_offset=-1*((size-1)/2);
				    while(x_offset<=(size-1)/2)
				    {
					    //NOW WE HAVE TO MAKE SAVE VARIABLES
					    int x_acces=x+x_offset;
					    int y_acces=y+y_offset;

					    if(x_acces<0){x_acces=-1*x_acces;}
					    if(y_acces<0){y_acces=-1*y_acces;}

					    if(x_acces>=dim_x){x_acces=dim_x-1 - ( x_acces-(dim_x-1) ) ; }
					    if(y_acces>=dim_y){y_acces=dim_y-1 - ( y_acces-(dim_y-1) ) ; }

					    sum0=sum0+image(x_acces,y_acces);

					    x_offset++;
				    }
				    y_offset++;
			    }
			
			    //SUM IS THE AVERAGE OF THE PIXELS NEIGBOURHOOD
			    sum0=sum0/(size*size);

			    float sd0=0;

			    y_offset=-1*((size-1)/2);
			    while(y_offset<=(size-1)/2)
			    {
				    int x_offset=-1*((size-1)/2);
				    while(x_offset<=(size-1)/2)
				    {
					    //NOW WE HAVE TO MAKE SAVE VARIABLES
					    int x_acces=x+x_offset;
					    int y_acces=y+y_offset;

					    if(x_acces<0){x_acces=-1*x_acces;}
					    if(y_acces<0){y_acces=-1*y_acces;}
	
					    if(x_acces>=dim_x){x_acces=dim_x-1 - ( x_acces-(dim_x-1) ) ; }
					    if(y_acces>=dim_y){y_acces=dim_y-1 - ( y_acces-(dim_y-1) ) ; }

					    sd0=sd0+  ((sum0-image(x_acces,y_acces))*(sum0-image(x_acces,y_acces)));

					    x_offset++;
				    }
				    y_offset++;
			    }
			
			    temp_image(x,y)=sqrt(sd0/(size*size-1));

             }
        }
    }//end of parallelisation

    vector_of_images.push_back(temp_image);

    if(feature_to_image==true)
    {
        export_images("SDV ",size,1);
    }
}

void FeatureExtractorGray::canny_edge_detector(double scale1,double scale2,double scale3)
{
    gaussianSmoothing(srcImageRange(image),destImage(temp_image2),scale3);

    temp_image=empty_image;

    cannyEdgeImage(srcImageRange(temp_image2),destImage(temp_image),scale1,scale2,1);

    gaussianSmoothing(srcImageRange(temp_image),destImage(temp_image),scale3);

    vector_of_images.push_back(temp_image);

    if(feature_to_image==true)
    {
        export_images("canny",scale1,scale2,scale3,1);
    }
}

void FeatureExtractorGray::boundary_strength(double scale)
{
    //Images to store the boundary tensor for each color channel
    vigra::FVector3Image bt_0(dim_x,dim_y);

    //perform boundaryTensor1
    vigra::boundaryTensor1(vigra::srcImageRange(image), vigra::destImage(bt_0), scale);

    //Images to store the trace
    vigra::tensorTrace(vigra::srcImageRange(bt_0), vigra::destImage(temp_image));

    vector_of_images.push_back(temp_image);
    if(feature_to_image==true)
    {
        export_images("boundary_tensor",scale,1);
    }
}

void FeatureExtractorGray::boundary_strength_as(double scale)
{
    //Images to store the boundary tensor for each color channel
    vigra::FVector3Image bt_0(dim_x,dim_y);

    //perform boundaryTensor1
    vigra::boundaryTensor1(vigra::srcImageRange(as_image), vigra::destImage(bt_0), scale);

    //Images to store the trace
    vigra::tensorTrace(vigra::srcImageRange(bt_0), vigra::destImage(temp_image));

    vector_of_images.push_back(temp_image);
    if(feature_to_image==true)
    {
        export_images("AS boundary_tensor",scale,1);
    }
}

vigra::FImage FeatureExtractorGray::decision_feature()
{
    vigra::FImage feature_image(dim_x,dim_y);

    //SCORE IMAGES
    vigra::FImage stxx(dim_x,dim_y), stxy(dim_x,dim_y), styy(dim_x,dim_y);
    hessianMatrixOfGaussian(srcImageRange(as_image),destImage(stxx), destImage(stxy), destImage(styy), 1.5);

    #pragma omp parallel
    {
        //visit all pixels
        #pragma omp for
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
	        {
                float lambda;
                float d=((stxx(x,y)-styy(x,y))/2);
                lambda= ((stxx(x,y)+styy(x,y))/2) - sqrt(d*d+ stxy(x,y)*stxy(x,y));
                feature_image(x,y)=fabs(lambda);
            }
        }
    }//end of parallelisation
 
    return feature_image;
}

//if bool to_file is true a feature file is saved
int FeatureExtractorGray::extract(bool to_file,std::string param_file_name)
{
    // some important parameters should not need to be changed in the source code
    // so let's load them from a parameter file
    ParameterFile paramFile;

    if( !paramFile.load(param_file_name) )
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        return 0;
    }

	Parameter<bool> grayvalue;
	grayvalue.assign("", "grayvalue", false);
    grayvalue.load(paramFile,"config");

//    std::cout << "Parameter grayvalue: " << grayvalue << std::endl;

	Parameter<float> gradient_magnitude1;
	gradient_magnitude1.assign("", "gradient_magnitude1", 1.0f);
    gradient_magnitude1.load(paramFile,"config");

//    std::cout << "Parameter gradient_magnitude1: " << gradient_magnitude1 << std::endl;

	Parameter<float> gradient_magnitude2;
	gradient_magnitude2.assign("", "gradient_magnitude2", 1.5f);
    gradient_magnitude2.load(paramFile,"config");

//    std::cout << "Parameter gradient_magnitude2: " << gradient_magnitude2 << std::endl;

	Parameter<float> gradient_magnitude3;
	gradient_magnitude3.assign("", "gradient_magnitude3", 3.0f);
    gradient_magnitude3.load(paramFile,"config");

//    std::cout << "Parameter gradient_magnitude3: " << gradient_magnitude3 << std::endl;



	Parameter<float> eigenvalues_structure_tensor1a;
	eigenvalues_structure_tensor1a.assign("", "eigenvalues_structure_tensor1a", 1.0f);
    eigenvalues_structure_tensor1a.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor1a: " << eigenvalues_structure_tensor1a << std::endl;

	Parameter<float> eigenvalues_structure_tensor1b;
	eigenvalues_structure_tensor1b.assign("", "eigenvalues_structure_tensor1b", 0.6f);
    eigenvalues_structure_tensor1b.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor1b: " << eigenvalues_structure_tensor1b << std::endl;

	Parameter<float> eigenvalues_structure_tensor2a;
	eigenvalues_structure_tensor2a.assign("", "eigenvalues_structure_tensor2a", 1.5f);
    eigenvalues_structure_tensor2a.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor2a: " << eigenvalues_structure_tensor2a << std::endl;

	Parameter<float> eigenvalues_structure_tensor2b;
	eigenvalues_structure_tensor2b.assign("", "eigenvalues_structure_tensor2b", 0.9f);
    eigenvalues_structure_tensor2b.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor2b: " << eigenvalues_structure_tensor2b << std::endl;

	Parameter<float> eigenvalues_structure_tensor3a;
	eigenvalues_structure_tensor3a.assign("", "eigenvalues_structure_tensor3a", 2.0f);
    eigenvalues_structure_tensor3a.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor3a: " << eigenvalues_structure_tensor3a << std::endl;

	Parameter<float> eigenvalues_structure_tensor3b;
	eigenvalues_structure_tensor3b.assign("", "eigenvalues_structure_tensor3b", 1.2f);
    eigenvalues_structure_tensor3b.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor3b: " << eigenvalues_structure_tensor3b << std::endl;



	Parameter<float> eigenvalues_hessian_matrix1;
	eigenvalues_hessian_matrix1.assign("", "eigenvalues_hessian_matrix1", 1.0f);
    eigenvalues_hessian_matrix1.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_hessian_matrix1: " << eigenvalues_hessian_matrix1 << std::endl;

	Parameter<float> eigenvalues_hessian_matrix2;
	eigenvalues_hessian_matrix2.assign("", "eigenvalues_hessian_matrix2", 1.5f);
    eigenvalues_hessian_matrix2.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_hessian_matrix2: " << eigenvalues_hessian_matrix2 << std::endl;

	Parameter<float> eigenvalues_hessian_matrix3;
	eigenvalues_hessian_matrix3.assign("", "eigenvalues_hessian_matrix3", 2.0f);
    eigenvalues_hessian_matrix3.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_hessian_matrix3: " << eigenvalues_hessian_matrix3 << std::endl;



	Parameter<float> boundary_strength1;
	boundary_strength1.assign("", "boundary_strength1", 0.5f);
    boundary_strength1.load(paramFile,"config");

//    std::cout << "Parameter boundary_strength1: " << boundary_strength1 << std::endl;

	Parameter<float> boundary_strength2;
	boundary_strength2.assign("", "boundary_strength2", 1.0f);
    boundary_strength2.load(paramFile,"config");

//    std::cout << "Parameter boundary_strength2: " << boundary_strength2 << std::endl;

	Parameter<float> boundary_strength3;
	boundary_strength3.assign("", "boundary_strength3", 1.5f);
    boundary_strength3.load(paramFile,"config");

//    std::cout << "Parameter boundary_strength3: " << boundary_strength3 << std::endl;


	Parameter<bool> grayvalue_as;
	grayvalue_as.assign("", "grayvalue_as", false);
    grayvalue_as.load(paramFile,"config");

//    std::cout << "Parameter grayvalue_as: " << grayvalue_as << std::endl;

	Parameter<float> gradient_magnitude_as1;
	gradient_magnitude_as1.assign("", "gradient_magnitude_as1", 1.0f);
    gradient_magnitude_as1.load(paramFile,"config");

//    std::cout << "Parameter gradient_magnitude_as1: " << gradient_magnitude_as1 << std::endl;

	Parameter<float> gradient_magnitude_as2;
	gradient_magnitude_as2.assign("", "gradient_magnitude_as2", 1.5f);
    gradient_magnitude_as2.load(paramFile,"config");

//    std::cout << "Parameter gradient_magnitude_as2: " << gradient_magnitude_as2 << std::endl;

	Parameter<float> gradient_magnitude_as3;
	gradient_magnitude_as3.assign("", "gradient_magnitude_as3", 3.0f);
    gradient_magnitude_as3.load(paramFile,"config");

//    std::cout << "Parameter gradient_magnitude_as3: " << gradient_magnitude_as3 << std::endl;



	Parameter<float> eigenvalues_structure_tensor_as1a;
	eigenvalues_structure_tensor_as1a.assign("", "eigenvalues_structure_tensor_as1a", 1.0f);
    eigenvalues_structure_tensor_as1a.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor_as1a: " << eigenvalues_structure_tensor_as1a << std::endl;

	Parameter<float> eigenvalues_structure_tensor_as1b;
	eigenvalues_structure_tensor_as1b.assign("", "eigenvalues_structure_tensor_as1b", 0.6f);
    eigenvalues_structure_tensor_as1b.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor_as1b: " << eigenvalues_structure_tensor_as1b << std::endl;

	Parameter<float> eigenvalues_structure_tensor_as2a;
	eigenvalues_structure_tensor_as2a.assign("", "eigenvalues_structure_tensor_as2a", 1.5f);
    eigenvalues_structure_tensor_as2a.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor_as2a: " << eigenvalues_structure_tensor_as2a << std::endl;

	Parameter<float> eigenvalues_structure_tensor_as2b;
	eigenvalues_structure_tensor_as2b.assign("", "eigenvalues_structure_tensor_as2b", 0.9f);
    eigenvalues_structure_tensor_as2b.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor_as2b: " << eigenvalues_structure_tensor_as2b << std::endl;

	Parameter<float> eigenvalues_structure_tensor_as3a;
	eigenvalues_structure_tensor_as3a.assign("", "eigenvalues_structure_tensor_as3a", 2.0f);
    eigenvalues_structure_tensor_as3a.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor_as3a: " << eigenvalues_structure_tensor_as3a << std::endl;

	Parameter<float> eigenvalues_structure_tensor_as3b;
	eigenvalues_structure_tensor_as3b.assign("", "eigenvalues_structure_tensor_as3b", 1.2f);
    eigenvalues_structure_tensor_as3b.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_structure_tensor_as3b: " << eigenvalues_structure_tensor_as3b << std::endl;



	Parameter<float> eigenvalues_hessian_matrix_as1;
	eigenvalues_hessian_matrix_as1.assign("", "eigenvalues_hessian_matrix_as1", 1.0f);
    eigenvalues_hessian_matrix_as1.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_hessian_matrix_as1: " << eigenvalues_hessian_matrix_as1 << std::endl;

	Parameter<float> eigenvalues_hessian_matrix_as2;
	eigenvalues_hessian_matrix_as2.assign("", "eigenvalues_hessian_matrix_as2", 1.5f);
    eigenvalues_hessian_matrix_as2.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_hessian_matrix_as2: " << eigenvalues_hessian_matrix_as2 << std::endl;

	Parameter<float> eigenvalues_hessian_matrix_as3;
	eigenvalues_hessian_matrix_as3.assign("", "eigenvalues_hessian_matrix_as3", 2.0f);
    eigenvalues_hessian_matrix_as3.load(paramFile,"config");

//    std::cout << "Parameter eigenvalues_hessian_matrix_as3: " << eigenvalues_hessian_matrix_as3 << std::endl;



	Parameter<float> boundary_strength_as1;
	boundary_strength_as1.assign("", "boundary_strength_as1", 0.5f);
    boundary_strength_as1.load(paramFile,"config");

//    std::cout << "Parameter boundary_strength_as1: " << boundary_strength_as1 << std::endl;

	Parameter<float> boundary_strength_as2;
	boundary_strength_as2.assign("", "boundary_strength_as2", 1.0f);
    boundary_strength_as2.load(paramFile,"config");

//    std::cout << "Parameter boundary_strength_as2: " << boundary_strength_as2 << std::endl;

	Parameter<float> boundary_strength_as3;
	boundary_strength_as3.assign("", "boundary_strength_as3", 1.5f);
    boundary_strength_as3.load(paramFile,"config");

//    std::cout << "Parameter boundary_strength_as3: " << boundary_strength_as3 << std::endl;



    Parameter<float> difference_of_gaussians1a;
	difference_of_gaussians1a.assign("", "difference_of_gaussians1a", 6.0f);
    difference_of_gaussians1a.load(paramFile,"config");

//    std::cout << "Parameter difference_of_gaussians1a: " << difference_of_gaussians1a << std::endl;

    Parameter<float> difference_of_gaussians1b;
	difference_of_gaussians1b.assign("", "difference_of_gaussians1b", 0.8f);
    difference_of_gaussians1b.load(paramFile,"config");

//    std::cout << "Parameter difference_of_gaussians1b: " << difference_of_gaussians1b << std::endl;

    Parameter<float> difference_of_gaussians2a;
	difference_of_gaussians2a.assign("", "difference_of_gaussians2a", 3.0f);
    difference_of_gaussians2a.load(paramFile,"config");

//    std::cout << "Parameter difference_of_gaussians2a: " << difference_of_gaussians2a << std::endl;

    Parameter<float> difference_of_gaussians2b;
	difference_of_gaussians2b.assign("", "difference_of_gaussians2b", 0.4f);
    difference_of_gaussians2b.load(paramFile,"config");

//    std::cout << "Parameter difference_of_gaussians2b: " << difference_of_gaussians2b << std::endl;



	Parameter<float> laplace_of_gaussian1;
	laplace_of_gaussian1.assign("", "laplace_of_gaussian1", 1.0f);
    laplace_of_gaussian1.load(paramFile,"config");

//    std::cout << "Parameter laplace_of_gaussian1: " << laplace_of_gaussian1 << std::endl;

	Parameter<float> laplace_of_gaussian2;
	laplace_of_gaussian2.assign("", "laplace_of_gaussian2", 2.0f);
    laplace_of_gaussian2.load(paramFile,"config");

//    std::cout << "Parameter laplace_of_gaussian2: " << laplace_of_gaussian2 << std::endl;

	Parameter<float> laplace_of_gaussian3;
	laplace_of_gaussian3.assign("", "laplace_of_gaussian3", 3.0f);
    laplace_of_gaussian3.load(paramFile,"config");

//    std::cout << "Parameter laplace_of_gaussian3: " << laplace_of_gaussian3 << std::endl;



    Parameter<float> abs_difference_of_gaussians1a;
	abs_difference_of_gaussians1a.assign("", "abs_difference_of_gaussians1a", 6.0f);
    abs_difference_of_gaussians1a.load(paramFile,"config");

//    std::cout << "Parameter abs_difference_of_gaussians1a: " << abs_difference_of_gaussians1a << std::endl;

    Parameter<float> abs_difference_of_gaussians1b;
	abs_difference_of_gaussians1b.assign("", "abs_difference_of_gaussians1b", 0.8f);
    abs_difference_of_gaussians1b.load(paramFile,"config");

//    std::cout << "Parameter abs_difference_of_gaussians1b: " << abs_difference_of_gaussians1b << std::endl;

    Parameter<float> abs_difference_of_gaussians2a;
	abs_difference_of_gaussians2a.assign("", "abs_difference_of_gaussians2a", 3.0f);
    abs_difference_of_gaussians2a.load(paramFile,"config");

//    std::cout << "Parameter abs_difference_of_gaussians2a: " << abs_difference_of_gaussians2a << std::endl;

    Parameter<float> abs_difference_of_gaussians2b;
	abs_difference_of_gaussians2b.assign("", "abs_difference_of_gaussians2b", 0.4f);
    abs_difference_of_gaussians2b.load(paramFile,"config");

//    std::cout << "Parameter abs_difference_of_gaussians2b: " << abs_difference_of_gaussians2b << std::endl;

    //threshold -1 means no threshold set
//	Parameter<int> gray_threshold;
//	gray_threshold.assign("", "gray_threshold", -1);
//  gray_threshold.load(paramFile,"config");

    Parameter<int> eigenvalues_hessian1;
	eigenvalues_hessian1.assign("", "eigenvalues_hessian1", 0);
    eigenvalues_hessian1.load(paramFile,"config");

    Parameter<int> eigenvalues_hessian2;
	eigenvalues_hessian2.assign("", "eigenvalues_hessian2", 0);
    eigenvalues_hessian2.load(paramFile,"config");

    Parameter<int> eigenvalues_hessian3;
	eigenvalues_hessian3.assign("", "eigenvalues_hessian3", 0);
    eigenvalues_hessian3.load(paramFile,"config");


    Parameter<int> eigenvalues_structure1;
	eigenvalues_structure1.assign("", "eigenvalues_structure1", 0);
    eigenvalues_structure1.load(paramFile,"config");

    Parameter<int> eigenvalues_structure2;
	eigenvalues_structure2.assign("", "eigenvalues_structure2", 0);
    eigenvalues_structure2.load(paramFile,"config");

    Parameter<int> eigenvalues_structure3;
	eigenvalues_structure3.assign("", "eigenvalues_structure3", 0);
    eigenvalues_structure3.load(paramFile,"config");


    Parameter<int> eigenvalues_hessian_as1;
	eigenvalues_hessian_as1.assign("", "eigenvalues_hessian_as1", 0);
    eigenvalues_hessian_as1.load(paramFile,"config");

    Parameter<int> eigenvalues_hessian_as2;
	eigenvalues_hessian_as2.assign("", "eigenvalues_hessian_as2", 0);
    eigenvalues_hessian_as2.load(paramFile,"config");

    Parameter<int> eigenvalues_hessian_as3;
	eigenvalues_hessian_as3.assign("", "eigenvalues_hessian_as3", 0);
    eigenvalues_hessian_as3.load(paramFile,"config");


    Parameter<int> eigenvalues_structure_as1;
	eigenvalues_structure_as1.assign("", "eigenvalues_structure_as1", 0);
    eigenvalues_structure_as1.load(paramFile,"config");

    Parameter<int> eigenvalues_structure_as2;
	eigenvalues_structure_as2.assign("", "eigenvalues_structure_as2", 0);
    eigenvalues_structure_as2.load(paramFile,"config");

    Parameter<int> eigenvalues_structure_as3;
	eigenvalues_structure_as3.assign("", "eigenvalues_structure_as3", 0);
    eigenvalues_structure_as3.load(paramFile,"config");

    timeval start, end;
    gettimeofday(&start, 0);

    if (grayvalue!=false)
    {
        vector_of_images.push_back(image);
        std::cout<<"Compute feature grayvalue"<<std::endl;
    }

    if (gradient_magnitude1!=0.0f)
    {
        gradient_magnitude(gradient_magnitude1);
        std::cout<<"Compute feature gradient_magnitude"<<std::endl; 
    }
    if (gradient_magnitude2!=0.0f) gradient_magnitude(gradient_magnitude2);
    if (gradient_magnitude3!=0.0f) gradient_magnitude(gradient_magnitude3);


    if (eigenvalues_structure_tensor1a!=0.0f) 
    {
        eigenvalues_structure_tensor(eigenvalues_structure_tensor1a,eigenvalues_structure_tensor1b,eigenvalues_structure1);
        std::cout<<"Compute feature eigenvalues_structure_tensor"<<std::endl; 
    }
    if (eigenvalues_structure_tensor2a!=0.0f) eigenvalues_structure_tensor(eigenvalues_structure_tensor2a,eigenvalues_structure_tensor2b,
                                                                           eigenvalues_structure2);
    if (eigenvalues_structure_tensor3a!=0.0f) eigenvalues_structure_tensor(eigenvalues_structure_tensor3a,eigenvalues_structure_tensor3b,
                                                                           eigenvalues_structure3);


    if (eigenvalues_hessian_matrix1!=0.0f) 
    {
        eigenvalues_hessian_matrix(eigenvalues_hessian_matrix1,eigenvalues_hessian1);
        std::cout<<"Compute feature eigenvalues_hessian_matrix"<<std::endl; 
    }
    if (eigenvalues_hessian_matrix2!=0.0f) eigenvalues_hessian_matrix(eigenvalues_hessian_matrix2,eigenvalues_hessian2);
    if (eigenvalues_hessian_matrix3!=0.0f) eigenvalues_hessian_matrix(eigenvalues_hessian_matrix3,eigenvalues_hessian3);


    if (boundary_strength1!=0.0f) 
    {
        boundary_strength(boundary_strength1);
        std::cout<<"Compute feature boundary_strength"<<std::endl; 
    }
    if (boundary_strength2!=0.0f) boundary_strength(boundary_strength2);
    if (boundary_strength3!=0.0f) boundary_strength(boundary_strength3);


    if (grayvalue_as!=false)
    {
        vector_of_images.push_back(as_image);
        std::cout<<"Compute feature grayvalue_as"<<std::endl;
    }

    if (gradient_magnitude_as1!=0.0f) 
    {
        gradient_magnitude_as(gradient_magnitude_as1);
        std::cout<<"Compute feature gradient_magnitude_as"<<std::endl; 
    }
    if (gradient_magnitude_as2!=0.0f) gradient_magnitude_as(gradient_magnitude_as2);
    if (gradient_magnitude_as3!=0.0f) gradient_magnitude_as(gradient_magnitude_as3);

    if (eigenvalues_structure_tensor_as1a!=0.0f) 
    {
        eigenvalues_structure_tensor_as(eigenvalues_structure_tensor_as1a, eigenvalues_structure_tensor_as1b,eigenvalues_structure_as1);
        std::cout<<"Compute feature eigenvalues_structure_tensor_as"<<std::endl; 
    }
    if (eigenvalues_structure_tensor_as2a!=0.0f) eigenvalues_structure_tensor_as(eigenvalues_structure_tensor_as2a, eigenvalues_structure_tensor_as2b,
                                                                                 eigenvalues_structure_as2);
    if (eigenvalues_structure_tensor_as3a!=0.0f) eigenvalues_structure_tensor_as(eigenvalues_structure_tensor_as3a, eigenvalues_structure_tensor_as3b,
                                                                                 eigenvalues_structure_as3);

    if (eigenvalues_hessian_matrix_as1!=0.0f) 
    {
        eigenvalues_hessian_matrix_as(eigenvalues_hessian_matrix_as1,eigenvalues_hessian_as1);
        std::cout<<"Compute feature eigenvalues_hessian_matrix_as"<<std::endl; 
    }
    if (eigenvalues_hessian_matrix_as2!=0.0f) eigenvalues_hessian_matrix_as(eigenvalues_hessian_matrix_as2,eigenvalues_hessian_as2);
    if (eigenvalues_hessian_matrix_as3!=0.0f) eigenvalues_hessian_matrix_as(eigenvalues_hessian_matrix_as3,eigenvalues_hessian_as3);

    if (boundary_strength_as1!=0.0f)
    {
        boundary_strength_as(boundary_strength_as1);
        std::cout<<"Compute feature boundary_strength_as"<<std::endl; 
    }
    if (boundary_strength_as2!=0.0f) boundary_strength_as(boundary_strength_as2);
    if (boundary_strength_as3!=0.0f) boundary_strength_as(boundary_strength_as3);


    if (difference_of_gaussians1a!=0.0f) 
    {
        difference_of_gaussians(difference_of_gaussians1a,difference_of_gaussians1b);
        std::cout<<"Compute feature difference_of_gaussians"<<std::endl; 
    }
    if (difference_of_gaussians2a!=0.0f) difference_of_gaussians(difference_of_gaussians2a,difference_of_gaussians2b);

    if (laplace_of_gaussian1!=0.0f) 
    {
        laplace_of_gaussian(laplace_of_gaussian1);
        std::cout<<"Compute feature laplace_of_gaussian"<<std::endl; 
    }
    if (laplace_of_gaussian2!=0.0f) laplace_of_gaussian(laplace_of_gaussian2);
    if (laplace_of_gaussian3!=0.0f) laplace_of_gaussian(laplace_of_gaussian3);

    if (abs_difference_of_gaussians1a!=0.0f)
    {
        abs_difference_of_gaussians(abs_difference_of_gaussians1a,abs_difference_of_gaussians1b);
        std::cout<<"Compute feature abs_difference_of_gaussians"<<std::endl; 
    }
    if (abs_difference_of_gaussians2a!=0.0f) abs_difference_of_gaussians(abs_difference_of_gaussians2a,abs_difference_of_gaussians2b);

    vigra::FImage feature=decision_feature();

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

    std::cout<<"Nr of features:"<<vector_of_images.size()<<std::endl;

    if(to_file==true)
    {
        FILE *fp;
        std::string filename_and_path=filename;  //filename is path+filename
        std::string only_filename=get_filename(filename_and_path);
        std::string path_and_filename=path_to_save_feature_file;
        path_and_filename.append(only_filename);
        fp =fopen(path_and_filename.append(".bin").c_str(),"wb");
        float f0;  //float for images;

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                for(int v=0;v<(int)vector_of_images.size();v++)
                {
                    f0=vector_of_images[v](x,y);
                    fwrite(&f0,sizeof(float),1,fp);
                }
            }
        }

        fclose(fp);

        /*
        FILE *fp_decision_feature;
        path_and_filename.resize(path_and_filename.size()-4);
        fp_decision_feature =fopen(path_and_filename.append(".decision.bin").c_str(),"wb");

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                f0=feature(x,y);
                fwrite(&f0,sizeof(float),1,fp_decision_feature);
            }
        }

        fclose(fp_decision_feature);
        */
    }

    return vector_of_images.size();
}

void FeatureExtractorGray::extract_training_data(std::string path_labels,std::string filepath_classification_set,int& nr_of_pixel_features,
                                                 bool overwrite = false)
{
    std::cout<<"writing training data to: "<<filepath_classification_set<<std::endl;
    int counter0=0;
    int counter1=0;
    int c0=0;
    int c1=0;
    int nr_of_written_floats0=0;
    int nr_of_written_floats1=0;

    //here is some dirty stuff,the last line is readed twice(eof cant look in future),so we define some variables,what the last wirte was
    //and we check if the last read hast already been written ,and i use -1 as start value because there cant be a regular x (and y) with -1
    int last_read_x=-1;
    int last_read_y=-1;

    //FIRST WE HAVE A LOOK IF THERE IS A FILENAME "filename.dat"
    std::string path_and_filename_of_image=filename;
    std::string only_filename=get_filename(path_and_filename_of_image);
    std::string path_and_filename_of_label=path_labels;
    only_filename.append(".dat");
    path_and_filename_of_label.append(only_filename);

    //TRAINING FILE IS THE FILE WE TRY TO SEARCH,IF IT IS EXISTEND WE GOT TRAINING DATA ON THIS IMAGE
    std::ifstream training_file(path_and_filename_of_label.c_str());
    std::ifstream temp_training_file(path_and_filename_of_label.c_str());

    //string is just for testing stuff
    std::string teststring;
    temp_training_file>>teststring;

    //info file described below
    FILE *info;
    std::string filepath_classification_set_info=filepath_classification_set;
    filepath_classification_set_info.append(".info.bin");

    std::string filepath_classification_set_0=filepath_classification_set;
    filepath_classification_set_0.append(".0.bin");

    std::string filepath_classification_set_1=filepath_classification_set;
    filepath_classification_set_1.append(".1.bin");

    if(!training_file)
    {
        std::cout<<"no training data"<<std::endl;
    }
    else  //FILE IS EXISTEND
    {
        std::cout<<"training data is existent"<<std::endl;
        if(teststring.size()!=0)
        {
            //FILE IS NOT EMPTY

            //info file is used to check if the nr of features has changed
            info =fopen(filepath_classification_set_info.c_str(),"rb");

            if (info==NULL)
            {
                std::cout<<"No classification set"<<std::endl;                
            }
            else
            {
                std::cout<<"Classification set is existent"<<std::endl;                
                fread(&c0,sizeof(int),1,info);
                fread(&c1,sizeof(int),1,info);
                int nr_of_old_features;
                fread(&nr_of_old_features,sizeof(int),1,info);
                fclose(info);

                if (nr_of_old_features!=nr_of_pixel_features)
                {
                    std::cout<<"Nr of old features: "<<nr_of_old_features<<std::endl;
                    std::cout<<"Old classification files contain information with different nr of features"<<std::endl;
                    std::cout<<"Old classification files are deleted"<<std::endl;
                    remove(filepath_classification_set_0.c_str());
                    remove(filepath_classification_set_1.c_str());
                    c0=0;
                    c1=0;
                }
                else if (overwrite)
                {
                    std::cout<<"Old classification files are deleted"<<std::endl;
                    remove(filepath_classification_set_0.c_str());
                    remove(filepath_classification_set_1.c_str());
                    c0=0;
                    c1=0;
                }
            }

            //temp file is used to get one line out of the training file
            std::vector<int> temp;
            temp.resize(3);
            while(!training_file.eof())
            {
                training_file>>temp[0];
                training_file>>temp[1];
                training_file>>temp[2];

                if(temp[2]==0)
                {
                    //NOW WE CREATE A BINARY FILE training0.bin
                    FILE *fp0;
                    //with "a+b" we write /append at the end of the file
                    // If the file is nonexistend it will be created

                    fp0 =fopen(filepath_classification_set_0.c_str(),"a+b");

                    float f0;
                    for(int v=0;v<(int)vector_of_images.size();v++)
                    {
                        if(last_read_x!=temp[0] || last_read_y!=temp[1])
                        {
                            f0=vector_of_images[v](temp[0],temp[1]);
                            fwrite(&f0,sizeof(float),1,fp0);
                            nr_of_written_floats0=nr_of_written_floats0+1;
                            counter0++;
                        }
                    }
                    fclose(fp0);
                }
                    
                if(temp[2]==1)
                {
                    //NOW WE CREATE A BINARY FILE training1.bin
                    FILE *fp1;
                    //with "ab+" we write /append at the end of the file
                    // If the file is nonexistend it will be created

                    fp1 =fopen(filepath_classification_set_1.c_str(),"a+b");

                    float f0;
                    for(int v=0;v<(int)vector_of_images.size();v++)
                    {
                        if(last_read_x!=temp[0] || last_read_y!=temp[1])
                        {
                            f0=vector_of_images[v](temp[0],temp[1]);
                            fwrite(&f0,sizeof(float),1,fp1);
                            nr_of_written_floats1=nr_of_written_floats1+1;
                            counter1++;
                        }
                    }
                    fclose(fp1);
                }
                    
                if(temp[2]!=0 && temp[2] !=1)
                {
                    std::cout<<"something is wrong here,training data label must be 0 or 1"<<std::endl;
                }

                //HERE WE SET LAST READ TO THE LAST READED /WRITTEN POINTS TO AVOID THE LAST LINE
                //IS READED TWICE
                last_read_x=temp[0];
                last_read_y=temp[1];

            }
            training_file.close();


            // now we create a litte "info" file with some infos about the classification set
            // How many c0 ,c1 and how many features have been used
            info =fopen(filepath_classification_set_info.c_str(),"wb");

            std::cout<<"Nr of OLD training data with label=0:"<<c0<<std::endl;
            std::cout<<"Nr of NEW training data with label=0:"<<counter0/(nr_of_pixel_features)<<std::endl;

            c0=c0+(counter0/(nr_of_pixel_features));
            c1=c1+(counter1/(nr_of_pixel_features));
            std::cout<<"Nr of ALL training data with label=0:"<<c0<<std::endl;
            std::cout<<"Nr of ALL training data with label=1:"<<c1<<std::endl;

            int nr_of_features=(int)vector_of_images.size();

            fwrite(&c0,sizeof(int),1,info);
            fwrite(&c1,sizeof(int),1,info);
            fwrite(&nr_of_features,sizeof(int),1,info);

            fclose(info);
        }

        if(teststring.size()==0)
        {
            std::cout<<"training data is existend but empty"<<std::endl;
        }

        if(counter0!=counter1)
        {
            std::cout<<"WARNING NO EQUAl NUMBER in"<<filename<<std::endl;
        }

        if(nr_of_written_floats0!=nr_of_written_floats1)
        {
            std::cout<<"WARNING NO EQUAl NUMBER in"<<filename<<std::endl;
        }

    }

}
