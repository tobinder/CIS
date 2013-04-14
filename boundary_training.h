/*! \file boundary_training.h
 * \brief Boundary training.
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

#include "marray.hxx"
#include "boundary_data_structure.h"

/*
training data as file:
arc nr    mode
0         0         =>arcs[0] is labled as class 0 (no-boundary)
44        1         =>arcs[44] is a boundary class 1
*/

class Boundary_training_image
{
    private:

    marray::Marray<unsigned int> one_boundings;
    marray::Marray<unsigned int> two_boundings;

    std::vector< std::vector<point> > arcs;
    std::vector<point> junctions;

    cimg_library::CImg<unsigned char> image_to_train_on;

    //TODO image aus trainingsdaten erstellen
    cimg_library::CImg<unsigned char> menu;

    //vector_of_training_data is used to represent the training data during the training
    std::vector<int>    vector_training_data_boundaries;
    std::vector<int>    vector_training_data_junctions;

    std::string         filename_of_image; //filename of the actual image

    std::string         path_to_save_label_file;
    std::string         filepath_param_file;
    int                 posx;
    int                 posy;

    //Private methods:
    cimg_library::CImg<unsigned char> select_image_section(cimg_library::CImg<unsigned char> image,int display_x,int display_y);

    public:

    Boundary_training_image(std::string filepath_img,
                            std::string filepath_to_ws_region_image,
                            std::string filepath_out,
                            std::string param_file="parameters.txt");

    void training_data_to_file(void);
    void file_to_training_data(void);
    void do_training(void);

    int get_nr_of_boundaries(void);
    int get_nr_of_class_boundaries(int);

    int get_nr_of_junctions(void);
    int get_nr_of_class_junctions(int);
};

Boundary_training_image::Boundary_training_image(std::string filepath_img,
                                                std::string filepath_to_ws_region_image,
                                                std::string filepath_out,
                                                std::string param_file)
{
    //IMPORT RESULTS FROM HDF5 file

    vigra::BasicImage<unsigned int> ws_region_image;
    int dim_x, dim_y;

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

    //DO SOME VIGRA IMAGE IMPORT AND CONVERT IT BACK TO A CIMG IMAGE
    std::cout<<"start image importing..."<<std::endl;

    vigra::ImageImportInfo info_img(filepath_img.c_str());

    std::cout<<"dim_X="<<dim_x<<std::endl;
    std::cout<<"dim_Y="<<dim_y<<std::endl;

    //Check if grayscale or color
    if(info_img.isGrayscale())
    {
        vigra::BImage temp_img(dim_x,dim_y);
        cimg_library::CImg<unsigned char> temp_img_cimg(dim_x,dim_y,1,3);
        importImage(info_img, destImage(temp_img));

        //vigra=>CImg
        for(int x=0;x<dim_x;x++)
        {
            for(int y=0;y<dim_y;y++)
            {
                temp_img_cimg(x,y,0)=temp_img(x,y);
                temp_img_cimg(x,y,1)=temp_img(x,y);
                temp_img_cimg(x,y,2)=temp_img(x,y);
            }
        }

	    image_to_train_on=temp_img_cimg;
    }
    else
    {
        vigra::BRGBImage temp_img(dim_x,dim_y);
        cimg_library::CImg<unsigned char> temp_img_cimg(dim_x,dim_y,1,3);
        importImage(info_img, destImage(temp_img));

        //vigra=>CImg
        for(int x=0;x<dim_x;x++)
        {
            for(int y=0;y<dim_y;y++)
            {
                temp_img_cimg(x,y,0)=temp_img(x,y)[0];
                temp_img_cimg(x,y,1)=temp_img(x,y)[1];
                temp_img_cimg(x,y,2)=temp_img(x,y)[2];
            }
        }

	    image_to_train_on=temp_img_cimg;
    }

    // Initialization of menu
    //temp CImg file is used to avoid an error
	cimg_library::CImg<unsigned char> temp_menu(450,20+nr_of_classes*10,1,3,0);
	menu=temp_menu;

    //Set Filename
	filename_of_image=filepath_img;
	path_to_save_label_file=filepath_out;

    //Set parameters file
    filepath_param_file=param_file;

    //we have to set up the vectorof the training data;
    vector_training_data_boundaries.clear();
    vector_training_data_boundaries.resize( arcs.size(), 10 );  // 10 means that it has not been decided if it is a boundary or not (to initialize)
    vector_training_data_junctions.clear();
    vector_training_data_junctions.resize( junctions.size(), 10 );  // 10 means that it has not been decided if it is a boundary or not (to initialize)

    std::cout<<"...done"<<std::endl;
}

void Boundary_training_image::training_data_to_file(void)
{
    std::string filename_and_path=filename_of_image;
    std::string filename=get_filename(filename_and_path);
    std::string path_and_filename_to_save=path_to_save_label_file;

    //BOUNDARIES
    path_and_filename_to_save.append(filename);
    path_and_filename_to_save.append(".dat");

    std::ofstream training_file(path_and_filename_to_save.c_str());

    std::cout<<"boundary training data to file: "<<std::endl;
    std::cout<<"size: "<<vector_training_data_boundaries.size()<<std::endl;
    int i=0;
    while(i<(int)vector_training_data_boundaries.size())
    {
        if(vector_training_data_boundaries[i]!=10)
        {
        training_file<<i<<" "<<vector_training_data_boundaries[i]<< "\n";
        }
        i++;
    }
    training_file.close();

    //JUNCTIONS
    path_and_filename_to_save.resize(path_and_filename_to_save.size()-4);
    path_and_filename_to_save.append(".junctions.dat");

    std::ofstream training_file_junctions(path_and_filename_to_save.c_str());

    std::cout<<"junction training data to file: "<<std::endl;
    std::cout<<"size: "<<vector_training_data_junctions.size()<<std::endl;
    i=0;
    while(i<(int)vector_training_data_junctions.size())
    {
        if(vector_training_data_junctions[i]!=10)
        {
        training_file_junctions<<i<<" "<<vector_training_data_junctions[i]<< "\n";
        }
        i++;
    }
    training_file_junctions.close();
}

void Boundary_training_image::file_to_training_data(void)
{
    std::string filename_and_path=filename_of_image;
    std::string filename=get_filename(filename_and_path);
    std::string path_and_filename_to_save=path_to_save_label_file;

    //BOUNDARIES
    path_and_filename_to_save.append(filename);
    path_and_filename_to_save.append(".dat");

    std::ifstream training_file(path_and_filename_to_save.c_str());
    std::ifstream temp_training_file(path_and_filename_to_save.c_str());

    std::cout<<"boundary file to training data..."<<std::endl;

    //string is just for testing stuff
    std::string teststring;
    temp_training_file>>teststring;
	if(!training_file)
	{
		//VECTOR TRAINING DATA == EMPTY(!) TRAINING DATA
		vector_training_data_boundaries.resize(arcs.size(),10);
	}
	else  //FILE IS EXISTEND
	{
		if(teststring.size()!=0)
		{
			//FILE IS NOT EMPTY
			int i=0;
			//temp file is used to get one line  out of the training file
			std::vector<int> temp;
			temp.resize(2);
			while(!training_file.eof())
			{
                //the vector of the training data is
				training_file>>temp[0]>>temp[1];
                if (temp[0]>arcs.size()-1)
                {
                   std::cout<<"error in file_to_training_data(..)"<<std::endl;
                   std::cout<<"label for arc nr "<<temp[0]<<" exists. nr of arcs: "<<arcs.size()<<std::endl;
                }
                else vector_training_data_boundaries[temp[0]]=temp[1];
				i++;
			}
		}
		if(teststring.size()==0)
		{
			//VECTOR TRAINING DATA == EMTY(!) TRAINING DATA
			std::vector<int>  temp;
			vector_training_data_boundaries.resize(arcs.size(),10);
	    }
	}
	std::cout<<"...done"<<std::endl;

    //JUNCTIONS
    path_and_filename_to_save.resize(path_and_filename_to_save.size()-4);
    path_and_filename_to_save.append(".junctions.dat");

    std::ifstream training_file_junctions(path_and_filename_to_save.c_str());
    std::ifstream temp_training_file_junctions(path_and_filename_to_save.c_str());

    std::cout<<"junction file to training data..."<<std::endl;

    teststring.resize(0);
    temp_training_file_junctions>>teststring;
	if(!training_file_junctions)
	{
		//VECTOR TRAINING DATA == EMPTY(!) TRAINING DATA
		vector_training_data_junctions.resize(junctions.size(),10);
	}
	else  //FILE IS EXISTEND
	{
		if(teststring.size()!=0)
		{
			//FILE IS NOT EMPTY
			int i=0;
			//temp file is used to get one line  out of the training file
			std::vector<int> temp;
			temp.resize(2);
			while(!training_file_junctions.eof())
			{
                //the vector of the training data is
				training_file_junctions>>temp[0]>>temp[1];
                if (temp[0]>junctions.size()-1)
                {
                   std::cout<<"error in file_to_training_data(..)"<<std::endl;
                   std::cout<<"label for junctions nr "<<temp[0]<<" exists. nr of junctions: "<<junctions.size()<<std::endl;
                }
                else vector_training_data_junctions[temp[0]]=temp[1];
				i++;
			}
		}
		if(teststring.size()==0)
		{
			//VECTOR TRAINING DATA == EMPTY(!) TRAINING DATA
			vector_training_data_junctions.resize(junctions.size(),10);
	    }
	}
	std::cout<<"...done"<<std::endl;
}

cimg_library::CImg<unsigned char> Boundary_training_image::select_image_section(cimg_library::CImg<unsigned char> image,int display_x,int display_y)
{
    // canvas for our gui to draw on
    // CImg<type> name(dimx, dimy, dimz, colors)
    cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);

    // draw the current selection on canvas
    cimg_forXY(canvas,x,y)
    {
        if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
            canvas(x,y,0,0) = canvas(x,y,0,1) = canvas(x,y,0,2) = image(x + posx,y + posy);
    }

    // show the current canvas
    cimg_library::CImgDisplay selectionDisplay(canvas,"Select image section for labeling!");

    bool selected = false;

    // this loop allows to change the relative position and size of the image part which we are
    // selecting
    while(!selected)
    {
        canvas.assign( selectionDisplay.dimx(), selectionDisplay.dimy(), 1, 3);
        cimg_forXY(canvas,x,y)
        {
            if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
                canvas(x,y,0,0) = canvas(x,y,0,1) = canvas(x,y,0,2) = image(x+posx,y+posy);
        }

        canvas.display(selectionDisplay);

        selectionDisplay.wait();

        // move selection
        if( selectionDisplay.is_keyARROWUP )
        {
            posy -= 100;
            if (posy<0) posy=0;
        }
        if( selectionDisplay.is_keyARROWDOWN )
        {
            posy += 100;
            if (posy>(image.dimy()-canvas.dimy())) posy=image.dimy()-canvas.dimy();
        }
        if( selectionDisplay.is_keyARROWRIGHT )
        {
            posx += 100;
            if (posx>(image.dimx()-canvas.dimx())) posx=image.dimx()-canvas.dimx();
        }
        if( selectionDisplay.is_keyARROWLEFT )
        {
            posx -= 100;
            if (posx<0) posx=0;
        }

        // press "ENTER" if you are content with your selection
        if( selectionDisplay.is_keyENTER || selectionDisplay.is_closed)
        {
            selected = true;
        }

        // resize selection
        if( selectionDisplay.is_resized )
        {
            selectionDisplay.resize();
            if (selectionDisplay.dimx()>image.dimx()) selectionDisplay.resize(image.dimx(),selectionDisplay.dimy());
            if (selectionDisplay.dimy()>image.dimy()) selectionDisplay.resize(selectionDisplay.dimx(),image.dimy());
        }
    }

    // image part to process
    cimg_library::CImg<unsigned char> selected_image( selectionDisplay.dimx(), selectionDisplay.dimy(), 1, 3);
    cimg_forXY(selected_image, x, y)
    {
        if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
                selected_image(x,y,0,0) = selected_image(x,y,0,1) = selected_image(x,y,0,2) = image(x+posx,y+posy);
    }

    return selected_image;
}

void Boundary_training_image::do_training(void)
{
    /*
     *Creates GUI
     */
    typedef unsigned char color_type[3];
    color_type * color = new color_type[13];
    color[0][0] = 0;    color[0][1] = 255;  color[0][2] = 0;
    color[1][0] = 255;  color[1][1] = 0;    color[1][2] = 0;
    color[2][0] = 0;    color[2][1] = 100;  color[2][2] = 255;
    color[3][0] = 255;  color[3][1] = 255;  color[3][2] = 0;
    color[4][0] = 255;  color[4][1] = 0;    color[4][2] = 255;
    color[5][0] = 0;    color[5][1] = 255;  color[5][2] = 255;
    color[6][0] = 127;  color[6][1] = 255;  color[6][2] = 0;
    color[7][0] = 255;  color[7][1] = 127;  color[7][2] = 0;
    color[8][0] = 127;  color[8][1] = 0;    color[8][2] = 255;
    color[9][0] = 127;  color[9][1] = 127;  color[9][2] = 0;
    color[10][0] = 255;  color[10][1] = 255;  color[10][2] = 255;//undecided
    color[11][0] = 0;    color[11][1] = 255;  color[11][2] = 255;//junctions
    color[12][0] = 255;  color[12][1] = 255;  color[12][2] = 255;//white

    int x_mouse=0 ;
    int y_mouse=0 ;
    posx = 0;
    posy = 0;

    //wenn end of training true ist ist das training auf der aktuellen datei zuende
    //is_saved is used to remember if image has been saved after the last point added
    bool end_of_training=false;
    int  border=1;
    bool is_saved=true;
    bool modify_section=false;

    cimg_library::CImg<unsigned char> selected_image=image_to_train_on;
    cimg_library::CImg<unsigned char> men=menu;
    int factor=10;

    // some important parameters should not need to be changed in the source code
    // so let's load them from a parameter file
    ParameterFile paramFile;

    if( !paramFile.load(filepath_param_file.c_str()))
    {
        std::cout<<"Error: Parameter file could not be found!"<<std::endl;
        exit(-1);
    }

	Parameter<int> display_x;
	display_x.assign("", "display_x", 900);
    display_x.load(paramFile,"config");

	Parameter<int> display_y;
	display_y.assign("", "display_y", 600);
    display_y.load(paramFile,"config");

    while(end_of_training==false)
	{
        modify_section=false;

        if(image_to_train_on.dimx()>display_x || image_to_train_on.dimy()>display_y)
        {
            display_x=std::min(display_x(),image_to_train_on.dimx());
            display_y=std::min(display_y(),image_to_train_on.dimy());
            selected_image=select_image_section(image_to_train_on,display_x,display_y);
        }
        cimg_library::CImg<unsigned char> img=selected_image;

        	cimg_library::CImgDisplay  main_menu(menu,"Menu");
        cimg_library::CImgDisplay  zoom_disp(250,250,"ZOOM",0);
        	cimg_library::CImgDisplay  main_disp(selected_image,"Click a boundary");

        	main_menu.move(0,0);
        	main_disp.move(0,main_menu.dimy()+60);
        	zoom_disp.move(main_disp.dimx()+10,main_menu.dimy()+60);

        bool black=true;
        bool disp_b=true;
        bool disp_j=false;
        bool boundary_mode=true;
        bool junction_mode=false;        

        	while (modify_section==false && end_of_training==false)
	    {
            	//this disables to change the window sizes
        		if (main_disp.is_resized)
		    {
                main_disp.resize(main_disp);
		    }
        		if (main_menu.is_resized)
		    {
                main_menu.resize(main_menu);
		    }
        		if (zoom_disp.is_resized)
		    {
                zoom_disp.resize(zoom_disp);
		    }

	        //ZOOM IMAGE STUFF
	        //DEFINES THE ZOOM REGION TO GET A BIGGER ZOOM REGION YOU CAN MAKE factor bigger than 15
	        int x_0=x_mouse-factor;
	        int y_0=y_mouse-factor;
	        int x_1=x_mouse+factor;
	        int y_1=y_mouse+factor;
	        cimg_library::CImg<unsigned char> visu;

	        //VISU IS THE ZOOMED IMAGE
            if(border==0) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[0],0.5f);
        	    if(border==1) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[1],0.5f);
            if(border==2) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[2],0.5f);
            if(border==3) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[3],0.5f);
            if(border==4) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[4],0.5f);
            if(border==5) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[5],0.5f);
            if(border==6) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[6],0.5f);
            if(border==7) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[7],0.5f);
            if(border==8) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[8],0.5f);
            if(border==9) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[9],0.5f);

	        if(border==10) visu =img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[10],0.5f);

	        zoom_disp.display(visu);

	        //MENU STUFF
	        men=menu;
         	men.draw_text(2,30,"Coords (%d,%d)",color[12],0,1,11,x_mouse+posx,y_mouse+posy);
	        if(is_saved==true)
	        {
	            men.draw_text(2,0,"File is saved",color[0],0,1,11,1,1);
	        }
	        if(is_saved==false)
	        {
	            men.draw_text(2,0,"File is not saved",color[1],0,1,11,1,1);
	        }

	        if(boundary_mode==true)
	        {
	            men.draw_text(2,10,"Boundary",color[border],0,1,11,1,1);
	        }
	        else
	        {
	            men.draw_text(2,10,"Junction",color[border],0,1,11,1,1);
	        }

	        if(border==0) men.draw_text(2,20,"Mode=0",color[0],0,1,11,1,1);
		    if(border==1) men.draw_text(2,20,"Mode=1",color[1],0,1,11,1,1);
		    if(border==2) men.draw_text(2,20,"Mode=2",color[2],0,1,11,1,1);
		    if(border==3) men.draw_text(2,20,"Mode=3",color[3],0,1,11,1,1);
		    if(border==4) men.draw_text(2,20,"Mode=4",color[4],0,1,11,1,1);
		    if(border==5) men.draw_text(2,20,"Mode=5",color[5],0,1,11,1,1);
		    if(border==6) men.draw_text(2,20,"Mode=6",color[6],0,1,11,1,1);
		    if(border==7) men.draw_text(2,20,"Mode=7",color[7],0,1,11,1,1);
		    if(border==8) men.draw_text(2,20,"Mode=8",color[8],0,1,11,1,1);
		    if(border==9) men.draw_text(2,20,"Mode=9",color[9],0,1,11,1,1);
		    if(border==10) men.draw_text(2,20,"Mode=Rubber",color[10],0,1,11,1,1);

            {			
		        std::ostringstream Str;
		        Str << get_nr_of_boundaries();
		        std::string temp_string(Str.str());
		        std::string name="Boundaries: ";
		        name.append(temp_string);
		        men.draw_text(125,0,name.c_str(),color[12],0,1,11,1,1);
            }

            for (int i=0;i<nr_of_classes;i++)
            {
		        std::ostringstream Str;
		        Str << get_nr_of_class_boundaries(i);
		        std::string temp_string(Str.str());
                std::string name;
		        if (i==0) name="Boundaries Class 0: ";
		        if (i==1) name="Boundaries Class 1: ";
		        if (i==2) name="Boundaries Class 2: ";
		        if (i==3) name="Boundaries Class 3: ";
		        if (i==4) name="Boundaries Class 4: ";
		        if (i==5) name="Boundaries Class 5: ";
		        if (i==6) name="Boundaries Class 6: ";
		        if (i==7) name="Boundaries Class 7: ";
		        if (i==8) name="Boundaries Class 8: ";
		        if (i==9) name="Boundaries Class 9: ";
		        name.append(temp_string);
		        men.draw_text(125,10*(i+1),name.c_str(),color[i],0,1,11,1,1);
            }

            {			
		        std::ostringstream Str;
		        Str << get_nr_of_junctions();
		        std::string temp_string(Str.str());
		        std::string name="Junctions: ";
		        name.append(temp_string);
		        men.draw_text(300,0,name.c_str(),color[12],0,1,11,1,1);
            }

            for (int i=0;i<2;i++)
            {
		        std::ostringstream Str;
		        Str << get_nr_of_class_junctions(i);
		        std::string temp_string(Str.str());
                std::string name;
		        if (i==0) name="Junctions Class 0: ";
		        if (i==1) name="Junctions Class 1: ";
		        name.append(temp_string);
		        men.draw_text(300,10*(i+1),name.c_str(),color[i],0,1,11,1,1);
            }

		    main_menu.display(men);

		    //loop to display all boundaries if selected
		    int i=0;
		    img=selected_image;
		    if(disp_b)//TODO a boolen could be insertet to avoid the "display all points--while loop" if there are no new points to add
            {
                while(i<(int)arcs.size())//for each training data
                {
                    int size_of_this_arc=arcs[i].size();
                    if(vector_training_data_boundaries[i]==0)//A BOUNDARY CLASS 0
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[0],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==1)//A BOUNDARY CLASS 1
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[1],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==2)//A BOUNDARY CLASS 2
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[2],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==3)//A BOUNDARY CLASS 3
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[3],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==4)
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[4],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==5)
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[5],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==6)
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[6],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==7)
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[7],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==8)
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[8],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==9)
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[9],1);
                        }
                    }
                    else if(vector_training_data_boundaries[i]==10) //AN UNDECIDED BOUNDARY
                    {
                        //DRAW BOUNDARY
                        //each point of this arc
                        for(int j=0;j<size_of_this_arc;j++)
                        {
                            img.draw_point(arcs[i][j].x-posx,arcs[i][j].y-posy,color[10],0.5f);
                        }
                    }
                    if(vector_training_data_boundaries[i]>10 || (vector_training_data_boundaries[i]>(nr_of_classes-1) &&
                        vector_training_data_boundaries[i]!=10) || vector_training_data_boundaries[i]<0)
                    {
                        std::cout<<"ERROR BOUNDARY HAS BEEN CLASSIFIED WRONG"<<std::endl;
                    }
                    i++;
                }
            }
            
		    //loop to display all junctions if selected
            if(disp_j)
            {
                while(i<(int)junctions.size())//for each training data
                {
                    if(vector_training_data_junctions[i]==0)//A JUNCTION CLASS 0
                    {
                        //DRAW JUNCTION
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y-1-posy,color[0],1);
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y-posy,color[0],1);
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y+1-posy,color[0],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y-1-posy,color[0],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y-posy,color[0],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y+1-posy,color[0],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y-1-posy,color[0],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y-posy,color[0],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y+1-posy,color[0],1);
                    }
                    else if(vector_training_data_junctions[i]==1)//A JUNCTION CLASS 1
                    {
                        //DRAW JUNCTION
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y-1-posy,color[1],1);
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y-posy,color[1],1);
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y+1-posy,color[1],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y-1-posy,color[1],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y-posy,color[1],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y+1-posy,color[1],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y-1-posy,color[1],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y-posy,color[1],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y+1-posy,color[1],1);
                    }
                    else if(vector_training_data_junctions[i]==10) //AN UNDECIDED JUNCTION
                    {
                        //DRAW JUNCTION
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y-1-posy,color[10],1);
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y-posy,color[10],1);
                        img.draw_point(junctions[i].x-1-posx,junctions[i].y+1-posy,color[10],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y-1-posy,color[10],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y-posy,color[10],1);
                        img.draw_point(junctions[i].x-posx,junctions[i].y+1-posy,color[10],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y-1-posy,color[10],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y-posy,color[10],1);
                        img.draw_point(junctions[i].x+1-posx,junctions[i].y+1-posy,color[10],1);
                    }
                    if(vector_training_data_junctions[i]>10 || (vector_training_data_junctions[i]>2 && vector_training_data_junctions[i]!=10) ||
                       vector_training_data_junctions[i]<0)
                    {
                        std::cout<<"ERROR JUNCTION HAS BEEN CLASSIFIED WRONG"<<std::endl;
                    }
                    i++;
                }
            }
            
		    //NOW ALL BORDER AND NOW BORDER FILES ARE IN THE IMAGE "img" WE CAN DISPLAY IT
		    main_disp.display(img);
		    main_disp.wait();

    		//this disables to change the window sizes
    		if (main_disp.is_resized)
			{
    			main_disp.resize(main_disp);
			}
    		if (main_menu.is_resized)
			{
    			main_menu.resize(main_menu);
			}
    		if (zoom_disp.is_resized)
			{
    			zoom_disp.resize(zoom_disp);
			}

            //0
            if(main_disp.key==48)
		    {
		        border=0;
                end_of_training=false;
		    }

		    //1
            if(main_disp.key==49)
		    {
        		border=1;
        		end_of_training=false;
		    }

		    //2
            if(main_disp.key==50 && boundary_mode)
		    {
        		border=2;
		        end_of_training=false;
		    }

		    //3
            if(main_disp.key==51 && boundary_mode)
		    {
        		border=3;
        		end_of_training=false;
		    }

            //4
            if(main_disp.key==52 && boundary_mode)
		    {
        		border=4;
        		end_of_training=false;
		    }

		    //5
            if(main_disp.key==53 && boundary_mode)
		    {
		        border=5;
        		end_of_training=false;
		    }

		    //6
            if(main_disp.key==54 && boundary_mode)
		    {
        		border=6;
        		end_of_training=false;
		    }

		    //7
            if(main_disp.key==55 && boundary_mode)
		    {
        		border=7;
        		end_of_training=false;
		    }

		    //8
            if(main_disp.key==56 && boundary_mode)
		    {
        		border=8;
        		end_of_training=false;
		    }
	
        	//9
            if(main_disp.key==57 && boundary_mode)
		    {
        		border=9;
		        end_of_training=false;
		    }

    		//IF C (CLEAR) KEY IS PRESSED
    		if(main_disp.key==99||main_disp.key==67)
    		{
        		border=10;
		        end_of_training=false;
    		}

    		//TESTS IF Q-KEY (QUIT) IS PRESSED OR WINDOWS IS CLOSED
    		if(main_disp.key==113||main_disp.key==81 || main_disp.is_closed)
    		{
        		end_of_training=true;
    		}

            //TEST IF M-KEY (MODIFY SECTION) IS PRESSED
    		if(main_disp.key==109||main_disp.key==77)
    		{
    			training_data_to_file();
    			is_saved=true;
    			main_disp.close();
    			main_menu.close();
    			zoom_disp.close();
        		modify_section=true;
    		}
            
    		//TESTS IF B-KEY (BOUNDARY) IS PRESSED
    		if(main_disp.key==98||main_disp.key==66)
			{
			    bool twice=false;
                if (disp_b==false&&twice==false)
                {
                    disp_b=true;
                    twice=true;
                }
                if (disp_b==true&&twice==false)
                {
                    disp_b=false;
                    twice=true;
                }

                if(junction_mode) border=0;

                disp_j=false;
                boundary_mode=true;//Boundary mode
                junction_mode=false;
			}

    		//TESTS IF J-KEY (Junction) IS PRESSED
    		if(main_disp.key==106||main_disp.key==74)
			{
			    bool twice=false;
                if (disp_j==false&&twice==false)
                {
                    disp_j=true;
                    twice=true;
                }
                if (disp_j==true&&twice==false)
                {
                    disp_j=false;
                    twice=true;
                }

                if(boundary_mode) border=0;

                disp_b=false;
                boundary_mode=false;
                junction_mode=true;//Junction mode
			}
            
    		//TESTS IF N-KEY ( NO-BOUNDARY) IS PRESSED
    		if(main_disp.key==110||main_disp.key==78)
			{
			    bool twice=false;
                if (black==false&&twice==false)
                {
                    black=true;
                    twice=true;
                    color[10][0]=0;
                    color[10][1]=0;
                    color[10][2]=0;
                }
                if (black==true&&twice==false)
                {
                    black=false;
                    twice=true;
                    color[10][0]=255;
                    color[10][1]=255;
                    color[10][2]=255;
                }
			}

		    //TESTS IF S-KEY (Save Training data) IS PRESSED
		    if(main_disp.key==115||main_disp.key==83)
			{
			    training_data_to_file();
			    is_saved=true;
			    end_of_training=false;
			}

		    //Test if the mouse button is clicked on the image area */
     		if (main_disp.mouse_y>=0 && main_disp.mouse_x>=0)
			{
                //mouse stored to  x_mouse and y_mouse to use it out of the loop
			    x_mouse = main_disp.mouse_x;
			    y_mouse = main_disp.mouse_y;
			    if(main_disp.button||main_disp.key==97||main_disp.key==65)//mouse clicked or A-KEY (add) pressed
				{
                    if(boundary_mode)
                    {
				        bool found_point=false;
				        int found_at=-1;
				        //HERE WE HAVE TO CHECK IF THE POINT (X,Y) IS A POINT OF AN ARC
				        for(int arc_count=0;arc_count<(int)arcs.size() && found_point==false;arc_count++)
				        {
				            for(int point_count=0;point_count<(int)arcs[arc_count].size() && found_point==false;point_count++)
				            {
				                 if(arcs[arc_count][point_count].x==x_mouse+posx && arcs[arc_count][point_count].y==y_mouse+posy)
				                 {
				                     found_point=true;
				                     found_at=arc_count;
                                     /*
				                     std::cout<<"arc found: arc nr="<<arc_count<<std::endl;
				                     std::cout<<"points in this arc:"<<std::endl;
				                     for(int p_c=0; p_c<(int)arcs[arc_count].size() ;p_c++)
				                     {
				                         std::cout<<"x:"<<arcs[arc_count][p_c].x<<" y:"<<arcs[arc_count][p_c].y<<std::endl;
				                     }
                                     */
				                 }
				                 else
				                 {
				                     found_point=false;
				                 }
                            }
				        }

				        if(found_point==true)
				        {
					        vector_training_data_boundaries[found_at]=border;
                            //is_saved is used to see if the image must be saved or not
                            is_saved=false;
				        }
                    }
                    else
                    {
				        bool found_point=false;
				        int found_at=-1;
				        //HERE WE HAVE TO CHECK IF THE POINT (X,Y) IS A JUNCTION
				        for(int junc_count=0;junc_count<(int)junctions.size() && found_point==false;junc_count++)
				        {
			                 if(junctions[junc_count].x-1==x_mouse+posx && junctions[junc_count].y-1==y_mouse+posy ||
                                junctions[junc_count].x-1==x_mouse+posx && junctions[junc_count].y==y_mouse+posy ||
                                junctions[junc_count].x-1==x_mouse+posx && junctions[junc_count].y+1==y_mouse+posy ||
                                junctions[junc_count].x==x_mouse+posx && junctions[junc_count].y-1==y_mouse+posy ||
                                junctions[junc_count].x==x_mouse+posx && junctions[junc_count].y==y_mouse+posy ||
                                junctions[junc_count].x==x_mouse+posx && junctions[junc_count].y+1==y_mouse+posy ||
                                junctions[junc_count].x+1==x_mouse+posx && junctions[junc_count].y-1==y_mouse+posy ||
                                junctions[junc_count].x+1==x_mouse+posx && junctions[junc_count].y==y_mouse+posy ||
                                junctions[junc_count].x+1==x_mouse+posx && junctions[junc_count].y+1==y_mouse+posy)
			                 {
			                     found_point=true;
			                     found_at=junc_count;
                                 /*
			                     std::cout<<"junction found: junction nr="<<junc_count<<std::endl;
			                     std::cout<<"x:"<<junctions[junc_count].x<<" y:"<<junctions[junc_count].y<<std::endl;
                                 */
			                 }
			                 else
			                 {
			                     found_point=false;
			                 }
				        }

				        if(found_point==true)
				        {
                            //check if this junction is capable as training data
                            size_t nr_of_boundarys_crossing;
                            if(one_boundings(found_at,3) == 0) 
                            {
                                nr_of_boundarys_crossing = 3;
                            }
                            else
                            {
                                nr_of_boundarys_crossing = 4;
                            }

                            // now check the length of adjacent boundaries
                            for(size_t x=0; x<nr_of_boundarys_crossing; x++)
                            {
                                if (arcs[one_boundings(found_at,x)-1].size()<60)
                                {
                                    found_point=false;
                                    
                                }
                            }

    				        if(found_point==false) std::cout<<"this junction is not capable as training data, one arc is too short"<<std::endl;
                            else//capable
    				        {
                                vector_training_data_junctions[found_at]=border;
                                //is_saved is used to see if the image must be saved or not
                                is_saved=false;
                            }
				        }
                    }

				    main_disp.display(img);
                    end_of_training=false;
				}
			}
		}
	}
}

int Boundary_training_image::get_nr_of_boundaries(void)
{
    int i=0;
    int nr=0;
    while(i<(int)vector_training_data_boundaries.size())
    {
        if(vector_training_data_boundaries[i]!=10)
        {nr++;}
        i++;
    }
    return nr;
}

int Boundary_training_image::get_nr_of_class_boundaries(int class_nr)
{
    int i=0;
    int class_counter=0;
    while(i<(int)vector_training_data_boundaries.size())
	{
	    if(vector_training_data_boundaries[i]==class_nr)
	    {
		    class_counter++;
	    }
	    i++;
	}
    return class_counter;
}

int Boundary_training_image::get_nr_of_junctions(void)
{
    int i=0;
    int nr=0;
    while(i<(int)vector_training_data_junctions.size())
    {
        if(vector_training_data_junctions[i]!=10)
        {nr++;}
        i++;
    }
    return nr;
}

int Boundary_training_image::get_nr_of_class_junctions(int class_nr)
{
    int i=0;
    int class_counter=0;
    while(i<(int)vector_training_data_junctions.size())
	{
	    if(vector_training_data_junctions[i]==class_nr)
	    {
		    class_counter++;
	    }
	    i++;
	}
    return class_counter;
}
