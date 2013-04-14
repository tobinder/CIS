/*! \file pixel_training.h
 *  \brief Pixel training.
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
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <fstream>
#include <vector>
#include <algorithm>
#include <dirent.h>

#include "boundary_data_structure.h"
#include "cgp_structure.hxx"

//using namespace cimg_library;

/*
 *Pixel_training_image is used to Get training data for one(!) Image
 */
class Pixel_training_image
{
    private:

    cimg_library::CImg<unsigned char> image_to_train_on;
    cimg_library::CImg<unsigned char> as_image;
    cimg_library::CImg<unsigned char> prob_map;
    
    cimg_library::CImg<unsigned char> b_image;                  //Bubble image and 
    bool                bubble;                                 //bubble variable indicating the existence of a bubble image
    //TODO image aus trainingsdaten erstellen
    cimg_library::CImg<unsigned char> image_of_training_data;
    cimg_library::CImg<unsigned char> menu;
    std::string         training_data;
    std::vector< std::vector<int> > vector_training_data;
    std::string         filename_of_training_image;
    std::string         path_to_save_label_file;
    std::string         filepath_param_file;
    int                 nr_of_points;
    int                 nr_of_boundary_points;
    int                 nr_of_no_boundary_points;
    int                 posx;
    int                 posy;
    bool                prob;

    //colors
    std::vector<unsigned char> color_boundary;
    std::vector<unsigned char> color_no_boundary;
    std::vector<unsigned char> color_clear;
    //Formatierung der Datei::
    std::string   data_ending;   //WILL BE USED IN A LATER RELEASE TO SE WHAT IMAGE FORMAT IS USED
    bool          simple_mode;   //WILL BE USED IN A LATER RELEASE TO SET THE FORMAT OF THE Saved Training data

    //Private methods:
    cimg_library::CImg<unsigned char> select_image_section(cimg_library::CImg<unsigned char> image,int display_x,int display_y);

    //friend methods:
    friend bool x_and_y(std::vector<int> i,std::vector<int> j);
    friend void add_training_data(std::vector< std::vector<int> >& scr,int x_pixel,int y_pixel,bool border);
    friend bool sort_it(std::vector<int> i,std::vector<int> j );
    friend void remove_training_data(std::vector< std::vector<int> >& scr,int x_pixel,int y_pixel);
    
    public:
    
    Pixel_training_image(std::string filename, std::string filepath_as_image, std::string dest_path, std::string param_file="parameters.txt",
        std::string filepath_prob_map="", std::string filepath_segmentation="");

    void training_data_to_file(void);
    void file_to_training_data();
    //do_training is the main feature of the class,complete classification is in this method
    void do_training(void);
    void do_analysis(void);
    void draw_large_labels(vigra::IRGBImage &, int, int, const unsigned char *);
    std::vector< std::vector<int> >  get_vector_training_data(void);
    void set_vector_training_data(std::vector<std::vector<int> > new_vector_training_data);
    int get_nr_of_points(void);
    int get_nr_of_boundary_points(void);
    int get_nr_of_no_boundary_points(void);
};

/*
 *Implementation of class Pixel_training_image
 */
//friendmethoden
//TODO Y_Sort
bool sort_it(std::vector<int> i,std::vector<int> j)
{
    if(i[2]==j[2])
    {
        if(i[0]!=j[0]) return (i[0]<j[0]);
        else return (i[1]<j[1]);
    }
    else return i[2]>j[2];
}

bool x_and_y(std::vector<int> i,std::vector<int> j)
{
    return( (i[0]==j[0]) && i[1]==j[1] );
}

/*
 * add_training_data
*/
void add_training_data(std::vector<std::vector<int> >& scr ,int x_pixel,int y_pixel,int border)
{
    //NEW TRAINING DATA IS STORED IN A VECTOR
    std::vector<int> new_training_data;
    new_training_data.push_back(x_pixel);
    new_training_data.push_back(y_pixel);
    new_training_data.push_back(border);  //REMEMBER: BORDER==1 if Pixel is a Border, BORDER==0 if not
    //Seach if the pixel (x_pixel/y_pixel) has already been classified
    bool done=false;
    int i=0;
    while(i<(int)scr.size())
    {
        //checks IF x_pixel AND y_pixel are equal to the new x and y pixels
        // (classification can be different(!), will be checked a few lines later
        if( x_and_y(scr[i],new_training_data)==true)
        {
            //IF SAME PIXEL IS classified 2times with the same classes
            if(scr[i][2]==new_training_data[2])
            {
                //leave while loop
                i=scr.size();
                done=true;
            }
            //IF SAME PIXEL IS classified 2times with different(!) classes
            else
            {
                scr[i]=new_training_data;
                //leave while loop
                i=scr.size();
                done=true;
            }
        }
        //GO TO THE NEXT TRAINING Data Element
        else{i++;}
    }
    //IF done == false new_training data is not yet in scr => so we can add it and sort scr.
    if(done==false)
    {
        scr.push_back(new_training_data);
        std::sort(scr.begin(), scr.end(), sort_it);
    }
}

/*
 * removes_training_data removes one element with [i][0]==x_pixel && [i][0]==y_pixel
*/
void remove_training_data(std::vector< std::vector<int> >& scr,int x_pixel,int y_pixel)
{
    int i=0;
    while(i<(int)scr.size())
    {
        if(scr[i][0]==x_pixel && scr[i][1]==y_pixel)
        {
                scr.erase(scr.begin()+i);
                i=scr.size();
        }
        else
        {
                i++;
        }
    }
}

/*
 * training_data_to_file
 */
void Pixel_training_image::training_data_to_file(void)
{
    std::string filename_and_path_of_the_image=filename_of_training_image;
    std::string filename;
    std::string path_and_filename_to_save=path_to_save_label_file;
    filename=get_filename(filename_and_path_of_the_image);
    path_and_filename_to_save.append(filename);
    path_and_filename_to_save.append(".dat");

    std::ofstream training_file(path_and_filename_to_save.c_str());

    int i=0;
    while(i<(int)vector_training_data.size())
    {
        training_file <<vector_training_data[i][0]<<" "<<vector_training_data[i][1]<<" "<<vector_training_data[i][2]<< "\n"; //TODO STÖRT HIER DAS \n???
        i++;
    }

    training_file.close();
}

void Pixel_training_image::file_to_training_data(void)
{
    std::string filename_and_path=filename_of_training_image;
    std::string filename=get_filename(filename_and_path);

    std::string path_and_filename_to_save=path_to_save_label_file;

    DIR *pDir;
    pDir = opendir (path_and_filename_to_save.c_str());
    if (pDir == NULL)
    {
        std::cout<<"Error: Path "<<path_and_filename_to_save<<" does not exist!"<<std::endl;
        exit(-1);
    }

    path_and_filename_to_save.append(filename);
    path_and_filename_to_save.append(".dat");

    std::ifstream training_file(path_and_filename_to_save.c_str());
    std::ifstream temp_training_file(path_and_filename_to_save.c_str());

    //string is just for testing stuff
    std::string teststring;
    temp_training_file>>teststring;
    if(!training_file)
    {
        //VECTOR TRAINING DATA == EMPTY(!) TRAINING DATA
        std::vector<std::vector<int> > temp;
        vector_training_data=temp;
    }
    else  //FILE IS EXISTEND
    {
        if(teststring.size()!=0)
        {
            //FILE IS NOT EMPTY
            int i=0;
            //temp file is used to get one line  out of the training file
            std::vector<int> temp;
            temp.resize(3);
            while(!training_file.eof())
            {
                training_file>>temp[0]>>temp[1]>>temp[2];
                add_training_data(vector_training_data,temp[0],temp[1],temp[2]);
                i++;
            }
        }
        if(teststring.size()==0)
        {
            //VECTOR TRAINING DATA == EMPTY(!) TRAINING DATA
            std::vector<std::vector<int> > temp;
            vector_training_data=temp;
        }

        training_file.close();
        temp_training_file.close();
    }
}

Pixel_training_image::Pixel_training_image(std::string filename,std::string filepath_as_image,std::string dest_path,std::string param_file,
    std::string filepath_prob_map, std::string filepath_segmentation)
{
    cimg_library::CImg<unsigned char> temp_image_to_train_on(filename.c_str());
    cimg_library::CImg<unsigned char> temp_as_image(filepath_as_image.c_str());    
    if (filepath_prob_map == "") prob=false;
    else
    {
        std::ifstream file(filepath_prob_map.c_str());

        if(file)
        {
            std::cout<<"Probability map found!"<<std::endl;
            file.close();
            prob=true;
            cimg_library::CImg<unsigned char> temp_prob_map(filepath_prob_map.c_str());
            prob_map=temp_prob_map;

            filepath_segmentation.append(".h5");
            std::ifstream file2(filepath_segmentation.c_str());

            if(file2)
            {
                std::cout<<"Segmentation found!"<<std::endl;
                file2.close();

                //IMPORT RESULTS FROM HDF5 file
                std::vector< std::vector<point> > arcs;
                std::vector<point> junctions;

                marray::Marray<unsigned int> one_boundings;
                marray::Marray<unsigned int> two_boundings;
             
                vigra::BasicImage<unsigned int> ws_region_image;
                int dim_x, dim_y;

                load_cgp_data_structure(filepath_segmentation,
                                        ws_region_image,
                                        one_boundings,
                                        two_boundings,
                                        arcs,
                                        junctions,
                                        dim_x,
                                        dim_y,
                                        true);

                for(int i=0; i<arcs.size(); i++)
                {
                    //draw boundary
                    for(int j=0; j<arcs[i].size(); j++)
                    {
                        int x=arcs[i][j].x;
                        int y=arcs[i][j].y;
                        prob_map(x,y,0,0)=temp_prob_map(x,y);
                        prob_map(x,y,0,1)=temp_prob_map(x,y);
                        prob_map(x,y,0,2)=255;
                    }
                }
            }
        }
        else
        {
            prob=false;
            prob_map=temp_as_image;
        }
    }

    //Check if there is a bubble image
    /*
    filename.resize(filename.size() - 3);
    std::string filepath_b_image = filename.append("b.bmp");
    std::ifstream file_b(filepath_b_image.c_str());*/
    std::string filepath_b_image = filename;
    filepath_b_image.resize(filepath_b_image.size()-3);
    filepath_b_image.append("b.bmp");
    std::ifstream file_b(filepath_b_image.c_str());
    
    bubble = false;
    
    if(file_b)
    {
        std::cout << "Bubble image found!" << std::endl;
        cimg_library::CImg<unsigned char> temp_bubble_image(filepath_b_image.c_str());
        b_image = temp_bubble_image; 
        bubble = true;   
    }
    
    image_to_train_on=temp_image_to_train_on;
    as_image=temp_as_image;

    // Initialization of menu
    //temp CImg file is used to avoid an error
    //cimg_library::CImg<unsigned char> temp_menu(280,56,1,3,0);
    cimg_library::CImg<unsigned char> temp_menu(320,46,1,3,0);
    menu=temp_menu;

    //Set Filename
    filename_of_training_image=filename;
    path_to_save_label_file=dest_path;

    //Set parameters file
    filepath_param_file=param_file;

    // COLORS
    //RED
    color_boundary.push_back(255);
    color_boundary.push_back(0);
    color_boundary.push_back(0);
    //GREEN
    color_no_boundary.push_back(0);
    color_no_boundary.push_back(255);
    color_no_boundary.push_back(0);
    //BLUE
    color_clear.push_back(0);
    color_clear.push_back(0);
    color_clear.push_back(255);

    /*DO SOMETHING*/
}

cimg_library::CImg<unsigned char> Pixel_training_image::select_image_section(cimg_library::CImg<unsigned char> image,int display_x,int display_y)
{
    // canvas for our gui to draw on
    // CImg<type> name(dimx, dimy, dimz, colors)
    cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);
    if (canvas.dimx()>image.dimx()) canvas.resize(image.dimx(),canvas.dimy());
    if (canvas.dimy()>image.dimy()) canvas.resize(canvas.dimx(),image.dimy());

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

void Pixel_training_image::do_training(void)
{
    /*
    *Creates GUI
    */
    const unsigned char color1[] = { color_boundary[0],color_boundary[1],color_boundary[2] };
    const unsigned char color0[] = { color_no_boundary[0],color_no_boundary[1],color_no_boundary[2] };
    const unsigned char colorc[] = { color_clear[0],color_clear[1],color_clear[2]};

    const unsigned char color_white []={255,255,255};
    int x_mouse=0 ;
    int y_mouse=0 ;
    posx = 0;
    posy = 0;

    //if end_of_training is true (Q pressed or window closed) training an actual file is finished
    //is_saved is used to remember if image has been saved after the last point added
    bool end_of_training=false;
    //bool point_to_add=true;
    int  border=1;
    bool clear=false;
    bool is_saved=true;
    bool modify_section=false;
    bool prob_large=false;

    cimg_library::CImg<unsigned char> selected_image=image_to_train_on;
    cimg_library::CImg<unsigned char> selected_as_image=as_image;
    cimg_library::CImg<unsigned char> selected_prob_map=prob_map;
    cimg_library::CImg<unsigned char> selected_b_image = b_image;
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

    /*
    int grayvalue_choice=255;
 
    float grayvalue_histogram[256];
    for (int i=0;i<256;i++) grayvalue_histogram[i]=0;
 
    cimg_forXY(image_to_train_on, x, y)
    {
        grayvalue_histogram[image_to_train_on(x,y)]++;
    }
 
    for (int i=0;i<256;i++)
    {
        grayvalue_histogram[i]=grayvalue_histogram[i]/(float)(image_to_train_on.dimx()*image_to_train_on.dimy()/1000);
 
        float average=0.0f;
        for (int j=0;j<=i;j++) average+=(grayvalue_histogram[j]/(i+1));
 
        if (grayvalue_histogram[i]>average*5.0f && grayvalue_histogram[i]>5.0f && i>50)
        {
            grayvalue_choice=i;
            break;
        }
    }
    */
    int gray_threshold=255;//grayvalue_choice;
    /*
    float feature_threshold=4.0f;

    int old_block_grayvalue=-1;
    int block_grayvalue=0;
    int blocksize=100;

    for (int x=image_to_train_on.dimx()/2; x>blocksize; x-=blocksize)
    {
        for (int y=0; y<image_to_train_on.dimy(); y++)
        {
            for (int line=x; line>x-blocksize; line--)
                block_grayvalue+=image_to_train_on(line,y);
        }
        block_grayvalue=block_grayvalue/(blocksize*image_to_train_on.dimy());
        //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
        if (fabs(block_grayvalue-old_block_grayvalue)>10 && old_block_grayvalue>-1)
        {
            std::cout<<"low x found!"<<std::endl;
            if (x + display_x <= image_to_train_on.dimx()) posx=x;
            break;
        }
        old_block_grayvalue=block_grayvalue;
        block_grayvalue=0;
    }

    old_block_grayvalue=-1;

    for (int x=image_to_train_on.dimx()/2; x<image_to_train_on.dimx()-blocksize; x+=blocksize)
    {
        for (int y=0; y<image_to_train_on.dimy(); y++)
        {
            for (int line=x; line<x+blocksize; line++)
                block_grayvalue+=image_to_train_on(line,y);
        }
        block_grayvalue=block_grayvalue/(blocksize*image_to_train_on.dimy());
        //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
        if (fabs(block_grayvalue-old_block_grayvalue)>10 && old_block_grayvalue>-1)
        {
            std::cout<<"high x found!"<<std::endl;
            break;
        }
        old_block_grayvalue=block_grayvalue;
        block_grayvalue=0;
    }

    old_block_grayvalue=-1;

    for (int y=image_to_train_on.dimy()/2; y>blocksize; y-=blocksize)
    {
        for (int x=0; x<image_to_train_on.dimx(); x++)
        {
            for (int line=y; line>y-blocksize; line--)
                block_grayvalue+=image_to_train_on(x,line);
        }
        block_grayvalue=block_grayvalue/(blocksize*image_to_train_on.dimx());
        //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
        if (fabs(block_grayvalue-old_block_grayvalue)>10 && old_block_grayvalue>-1)
        {
            std::cout<<"low y found!"<<std::endl;
            if (y + display_y <= image_to_train_on.dimy()) posy=y;
            break;
        }
        old_block_grayvalue=block_grayvalue;
        block_grayvalue=0;
    }

    old_block_grayvalue=-1;

    for (int y=image_to_train_on.dimy()/2; y<image_to_train_on.dimy()-blocksize; y+=blocksize)
    {
        for (int x=0; x<image_to_train_on.dimx(); x++)
        {
            for (int line=y; line<y+blocksize; line++)
                block_grayvalue+=image_to_train_on(x,line);
        }
        block_grayvalue=block_grayvalue/(blocksize*image_to_train_on.dimx());
        //std::cout<<block_grayvalue<<" "<<old_block_grayvalue<<std::endl;
        if (fabs(block_grayvalue-old_block_grayvalue)>10 && old_block_grayvalue>-1)
        {
            std::cout<<"high y found!"<<std::endl;
            break;
        }
        old_block_grayvalue=block_grayvalue;
        block_grayvalue=0;
    }
    */

    while (end_of_training==false)
    {
        modify_section=false;

        //vigra::FImage feature_image;

        if(image_to_train_on.dimx()>display_x || image_to_train_on.dimy()>display_y)
        {
            selected_image=select_image_section(image_to_train_on,display_x,display_y);
            selected_as_image.resize(selected_image.dimx(),selected_image.dimy());
            selected_prob_map.resize(selected_image.dimx(),selected_image.dimy());
            selected_b_image.resize(selected_image.dimx(), selected_image.dimy());
            //feature_image.resize(selected_image.dimx(),selected_image.dimy());

            for (int x=0; x<selected_image.dimx(); x++)
                for (int y=0; y<selected_image.dimy(); y++)
                    if( x + posx >= 0 && x + posx <  as_image.dimx() && y + posy >= 0 && y + posy < as_image.dimy() )
                    {
                        selected_as_image(x,y,0,0) = selected_as_image(x,y,0,1) = selected_as_image(x,y,0,2) = as_image(x+posx,y+posy);                        
                        if (prob)
                        {
                            selected_prob_map(x,y,0,0) = prob_map(x+posx,y+posy,0,0);
                            selected_prob_map(x,y,0,1) = prob_map(x+posx,y+posy,0,1);
                            selected_prob_map(x,y,0,2) = prob_map(x+posx,y+posy,0,2);
                        }
                        if(bubble)
                        {
                            selected_b_image(x,y,0,0) = b_image(x+posx,y+posy,0,0);
                            selected_b_image(x,y,0,1) = b_image(x+posx,y+posy,0,1);
                            selected_b_image(x,y,0,2) = b_image(x+posx,y+posy,0,2);
                        }
                        //feature_image(x,y) = as_image(x+posx,y+posy);
                    }
        }
        /*
        else
        {
            feature_image.resize(image_to_train_on.dimx(),image_to_train_on.dimy());
            for (int x=0; x<selected_image.dimx(); x++)
                for (int y=0; y<selected_image.dimy(); y++)
                {
                    feature_image(x,y) = as_image(x,y);
                }
        }
        */

        cimg_library::CImg<unsigned char> img=selected_image;
        cimg_library::CImg<unsigned char> img2=selected_image;
        cimg_library::CImg<unsigned char> as_img=selected_as_image;
        cimg_library::CImg<unsigned char> prob_img=selected_prob_map;
        cimg_library::CImg<unsigned char> b_img=selected_b_image;
        vigra::IRGBImage output_image(image_to_train_on.dimx(),image_to_train_on.dimy());

        cimg_library::CImg<bool> temp_img(img.dimx(),img.dimy());

        /*
        //calculate feature eigenvalue hessian matrix
        vigra::FImage stxx(img.dimx(),img.dimy()), stxy(img.dimx(),img.dimy()), styy(img.dimx(),img.dimy());
        hessianMatrixOfGaussian(srcImageRange(feature_image),destImage(stxx), destImage(stxy), destImage(styy), 1.5);

        for(int y=0;y<img.dimy();y++)
        {
            for(int x=0;x<img.dimx();x++)
            {
                float lambda;
                float d=((stxx(x,y)-styy(x,y))/2);
                lambda= ((stxx(x,y)+styy(x,y))/2) - sqrt(d*d+ stxy(x,y)*stxy(x,y));
                feature_image(x,y)=fabs(lambda);
                temp_img(x,y)=false;
            }
        }
        */

        cimg_forXY(img, x, y)
        {
            //if (feature_image(x,y)>feature_threshold || selected_image(x,y)<gray_threshold)
            if (selected_image(x,y)<gray_threshold)
            {
                for (int x_fill=std::max(0,x-6); x_fill<std::min(img.dimx(),x+6); x_fill++)
                    for (int y_fill=std::max(0,y-6); y_fill<std::min(img.dimy(),y+6); y_fill++)
                        temp_img(x_fill,y_fill)=true;

            }
        }

        cimg_forXY(img, x, y)
        {
            if (temp_img(x,y)==false)
            {
                img(x,y,0,0)=img(x,y,0,1)=img(x,y,0,2)=255;
            }
            temp_img(x,y)=false;
        }

        cimg_library::CImgDisplay main_menu(menu,"Menu");
        cimg_library::CImgDisplay zoom_disp(250,250,"ZOOM",0);
        cimg_library::CImgDisplay zoom_as_disp(250,250,"ZOOM AS",0);
        cimg_library::CImgDisplay prob_disp(250,250,"PROB MAP",0);
        cimg_library::CImgDisplay zoom_b_disp(250,184,"ZOOM BUBBLES",0);
        cimg_library::CImgDisplay main_disp(selected_image,"Click a point");        

        main_menu.move(0,0);            
        main_disp.move(0,main_menu.dimy()+64);
        zoom_disp.move(main_disp.dimx()+80,main_menu.dimy()+64);
        zoom_as_disp.move(main_disp.dimx()+80,main_menu.dimy()+zoom_disp.dimy()+108);        
        if (prob) prob_disp.move(main_disp.dimx()+80,main_menu.dimy()+zoom_disp.dimy()+zoom_as_disp.dimy()+152);
        else prob_disp.close();
        if(bubble)
        {
            if(!prob) zoom_b_disp.move(main_disp.dimx()+80, main_menu.dimy()+zoom_disp.dimy()+zoom_as_disp.dimy()+108);
            if(prob) zoom_b_disp.move(main_disp.dimx()+80, main_menu.dimy()+zoom_disp.dimy()+zoom_as_disp.dimy()+prob_disp.dimy()+152);
        }        
        else zoom_b_disp.close();

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
            if (zoom_as_disp.is_resized)
            {
                zoom_as_disp.resize(zoom_as_disp);
            }
            if (prob) if (prob_disp.is_resized)
            {
                prob_disp.resize(prob_disp);
            }

            //ZOOM IMAGE STUFF
            //DEFINES THE ZOOM REGION TO GET A BIGGER ZOOM REGION, YOU CAN MAKE factor bigger than 15
            int x_0=x_mouse-factor;
            int y_0=y_mouse-factor;
            int x_1=x_mouse+factor;
            int y_1=y_mouse+factor;
            cimg_library::CImg<unsigned char> visu;
            cimg_library::CImg<unsigned char> visu_as;
            cimg_library::CImg<unsigned char> visu_prob;
            cimg_library::CImg<unsigned char> visu_b;

            if(border==1 && clear==false)
            {
                //VISU IS THE ZOOMED IMAGE
                visu=img2.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color1,0.5f);
                visu_as=as_img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color1,0.5f);
                visu_prob=prob_img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color1,0.5f);
                visu_b=b_img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color1,0.5f);
            }
            if(border==0 && clear==false)
            {
                visu=img2.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color0,0.5f);
                visu_as=as_img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color0,0.5f);
                visu_prob=prob_img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color0,0.5f);
                visu_b=b_img.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color0,0.5f);
            }
            if(clear==true)
            {
                //DRAW POINTS IN A 3x3 NEIGHBOURHOOD
                visu=img2.get_crop(x_0,y_0,x_1,y_1);
                //9 POINTS MUST BE DRAWN
                visu.draw_point(x_mouse-x_0+1,y_mouse-y_0,colorc,0.5f);
                visu.draw_point(x_mouse-x_0,y_mouse-y_0,colorc,0.5f);
                visu.draw_point(x_mouse-x_0-1,y_mouse-y_0,colorc,0.5f);
                visu.draw_point(x_mouse-x_0,y_mouse-y_0+1,colorc,0.5f);
                visu.draw_point(x_mouse-x_0,y_mouse-y_0-1,colorc,0.5f);
                visu.draw_point(x_mouse-x_0+1,y_mouse-y_0+1,colorc,0.5f);
                visu.draw_point(x_mouse-x_0-1,y_mouse-y_0+1,colorc,0.5f);
                visu.draw_point(x_mouse-x_0+1,y_mouse-y_0-1,colorc,0.5f);
                visu.draw_point(x_mouse-x_0-1,y_mouse-y_0-1,colorc,0.5f);

                //DRAW POINTS IN A 3x3 NEIGHBOURHOOD
                visu_as=as_img.get_crop(x_0,y_0,x_1,y_1);
                //9 POINTS MUST BE DRAWN
                visu_as.draw_point(x_mouse-x_0+1,y_mouse-y_0,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0,y_mouse-y_0,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0-1,y_mouse-y_0,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0,y_mouse-y_0+1,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0,y_mouse-y_0-1,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0+1,y_mouse-y_0+1,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0-1,y_mouse-y_0+1,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0+1,y_mouse-y_0-1,colorc,0.5f);
                visu_as.draw_point(x_mouse-x_0-1,y_mouse-y_0-1,colorc,0.5f);

                //DRAW POINTS IN A 3x3 NEIGHBOURHOOD
                visu_prob=prob_img.get_crop(x_0,y_0,x_1,y_1);
                
                //DRAW POINTS IN A 3x3 NEIGHBOURHOOD
                visu_b.get_crop(x_0,y_0,x_1,y_1);
                //9 POINTS MUST BE DRAWN
                visu_b.draw_point(x_mouse-x_0+1,y_mouse-y_0,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0,y_mouse-y_0,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0-1,y_mouse-y_0,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0,y_mouse-y_0+1,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0,y_mouse-y_0-1,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0+1,y_mouse-y_0+1,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0-1,y_mouse-y_0+1,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0+1,y_mouse-y_0-1,colorc,0.5f);
                visu_b.draw_point(x_mouse-x_0-1,y_mouse-y_0-1,colorc,0.5f);
            }

            zoom_disp.display(visu);
            zoom_as_disp.display(visu_as);
            if (prob) prob_disp.display(visu_prob);
            if (bubble) zoom_b_disp.display(visu_b);

            //MENU STUFF
            men=menu;
            men.draw_text(2,20,"Coords (%d,%d)",color_white,0,1,11,x_mouse+posx,y_mouse+posy);
            men.draw_text(2,30,"Gray Threshold: %d",color1,0,1,11,gray_threshold);
            //men.draw_text(125,30,"Calculated: %d",color1,0,1,11,grayvalue_choice);
            //men.draw_text(2,40,"Feature Threshold: %.1f",color_white,0,1,11,feature_threshold);

            if(is_saved==true)
            {
                men.draw_text(1,0,"File is saved",color0,0,1,11,1,1);
            }
            if(is_saved==false)
            {
                men.draw_text(1,0,"File is not saved",color1,0,1,11,1,1);
            }
            if(clear==false)
            {
                if(border==true)
                {
                    men.draw_text(2,10,"Mode=Boundary",color1,0,1,11,1,1);
                }
                if(border==false)
                {
                    men.draw_text(2,10,"Mode=No Boundary",color0,0,1,11,1,1);
                }
            }
            if(clear==true)
            {
                    men.draw_text(2,10,"Mode=Rubber",colorc,0,1,11,1,1);
            }

            {
                std::ostringstream Str;
                Str << get_nr_of_points();
                std::string temp_string(Str.str());
                std::string name="Points: ";
                name.append(temp_string);
                men.draw_text(125,0,name.c_str(),color_white,0,1,11,1,1);
            }

            {
                std::ostringstream Str;
                Str << get_nr_of_boundary_points();
                std::string temp_string(Str.str());
                std::string name="Boundary Points: ";
                name.append(temp_string);
                men.draw_text(125,10,name.c_str(),color1,0,1,11,1,1);
            }

            {
                std::ostringstream Str;
                Str << get_nr_of_no_boundary_points();
                std::string temp_string(Str.str());
                std::string name="Non-Boundary Points: ";
                name.append(temp_string);
                men.draw_text(125,20,name.c_str(),color0,0,1,11,1,1);
            }

            main_menu.display(men);

            if (prob_large) img=selected_prob_map;
            else img=selected_image;
            img2=selected_image;
            as_img=selected_as_image;

            cimg_forXY(img, x, y)
            {
                //if (feature_image(x,y)>feature_threshold || selected_image(x,y)<gray_threshold)
                if (selected_image(x,y)<gray_threshold)
                {
                    temp_img(x,y)=true;
                }
            }

            cimg_forXY(img, x, y)
            {
                if (temp_img(x,y)==false)
                {
                    img(x,y,0,0)=img(x,y,0,1)=img(x,y,0,2)=255;
                }
                temp_img(x,y)=false;
            }

            if(true)//TODO a boolen could be insertet to avoid the "display all points--while loop" if there are no new points to add
            {
                //loop to display all points:
                int i=0;

                while(i<(int)vector_training_data.size())//for each training data
                {
                    if(vector_training_data[i][2]==1) //A BOUNDARY POINT/PIXEL
                    {
                        //DRAW BOUNDARYS POINTS/PIXEL
                        if (vector_training_data[i][0]-posx>=0 && vector_training_data[i][1]-posy>=0 &&
                            vector_training_data[i][0]-posx<img.dimx() && vector_training_data[i][1]-posy<img.dimy())
                        {
                            img.draw_point(vector_training_data[i][0]-posx,vector_training_data[i][1]-posy,color1,1);
                            img2.draw_point(vector_training_data[i][0]-posx,vector_training_data[i][1]-posy,color1,1);
                            as_img.draw_point(vector_training_data[i][0]-posx,vector_training_data[i][1]-posy,color1,1);
                        }
                    }

                    if(vector_training_data[i][2]==0) //NOT A BOUNDARY POINT/PIXEL
                    {
                        //DRAW NO BOUNDARYS POINTS/PIXEL
                        if (vector_training_data[i][0]-posx>=0 && vector_training_data[i][1]-posy>=0 &&
                        vector_training_data[i][0]-posx<img.dimx() && vector_training_data[i][1]-posy<img.dimy())
                        {
                            img.draw_point(vector_training_data[i][0]-posx,vector_training_data[i][1]-posy,color0,1);
                            img2.draw_point(vector_training_data[i][0]-posx,vector_training_data[i][1]-posy,color0,1);
                            as_img.draw_point(vector_training_data[i][0]-posx,vector_training_data[i][1]-posy,color0,1);
                        }
                    }

                    //IF vector_training_data[2]!=1 and !=0 something has gone wrong
                    //1=BOUNDARY 0=NO BOUNDARY !!!
                    if(vector_training_data[i][2]!=0 && vector_training_data[i][2]!=1)
                    {
                        std::cout<<"ERROR PIXEL HAS BEEN CLASSIFIED WRONG"<<std::endl;
                    }

                    i++;
                }
            }

            //NOW ALL BORDER AND NO BORDER POINTS ARE IN THE IMAGE "img" WE CAN DISPLAY IT
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
            if (zoom_as_disp.is_resized)
            {
                zoom_as_disp.resize(zoom_as_disp);
            }
            if (prob) if (prob_disp.is_resized)
            {
                prob_disp.resize(prob_disp);
            }

            //TESTS IF Q-KEY (QUIT) IS PRESSED OR WINDOW IS CLOSED
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
                zoom_as_disp.close();                    
                if (prob) prob_disp.close();
                if(bubble) zoom_b_disp.close();
                modify_section=true;
            }

            //TESTS IF C-KEY (CLEAR) IS PRESSED
            if(main_disp.key==99||main_disp.key==67)
            {
                clear=true;
            }

            //TESTS IF B-KEY (BOUNDARY) IS PRESSED
            if(main_disp.key==98||main_disp.key==66)
            {
                border=1;
                clear=false;
                end_of_training=false;
            }

            //TESTS IF N-KEY (NO-BOUNDARY) IS PRESSED
            if(main_disp.key==110||main_disp.key==78)
            {
                border=0;
                clear=false;
                end_of_training=false;
            }

            //TESTS IF S-KEY (Save Training data) IS PRESSED
            if(main_disp.key==115||main_disp.key==83)
            {
                training_data_to_file();
                is_saved=true;
                end_of_training=false;
            }

            /*
            //TESTS IF X-KEY (higher feature_threshold) IS PRESSED
            if(main_disp.key==120||main_disp.key==88)
            {
                feature_threshold+=0.1;
                end_of_training=false;
            }

            //TESTS IF Y-KEY (lower feature_threshold) IS PRESSED
            if(main_disp.key==121||main_disp.key==89)
            {
                feature_threshold-=0.1;
                end_of_training=false;
            }
            */
            //TESTS IF W-KEY (lower gray_threshold by 1) IS PRESSED
            if(main_disp.key==119||main_disp.key==87)
            {
                gray_threshold-=1;
                end_of_training=false;
            }

            //TESTS IF E-KEY (higher gray_threshold by 1) IS PRESSED
            if(main_disp.key==101||main_disp.key==69)
            {
                gray_threshold+=1;
                end_of_training=false;
            }

            //TESTS IF R-KEY (lower gray_threshold by 10) IS PRESSED
            if(main_disp.key==114||main_disp.key==82)
            {
                gray_threshold-=10;
                end_of_training=false;
            }

            //TESTS IF T-KEY (higher gray_threshold by 10) IS PRESSED
            if(main_disp.key==116||main_disp.key==84)
            {
                gray_threshold+=10;
                end_of_training=false;
            }

            //TESTS IF P-KEY (print trainings data) IS PRESSED
            if(main_disp.key==112||main_disp.key==80)
            {
                std::string filename_and_path=filename_of_training_image;
                std::string filename=get_filename(filename_and_path);

                std::string filepath_output=path_to_save_label_file;
                filepath_output.append(filename);

                cimg_forXY(image_to_train_on, x, y)
                {
                    output_image(x,y)[0]=image_to_train_on(x,y,0,0);
                    output_image(x,y)[1]=image_to_train_on(x,y,0,1);
                    output_image(x,y)[2]=image_to_train_on(x,y,0,2);
                }

                //loop to display all points:
                int i=0;

                while(i<(int)vector_training_data.size())//for each training data
                {
                    if(vector_training_data[i][2]==1) //A BOUNDARY POINT/PIXEL
                    {
                        output_image(vector_training_data[i][0],vector_training_data[i][1])[0]=color1[0];
                        output_image(vector_training_data[i][0],vector_training_data[i][1])[1]=color1[1];
                        output_image(vector_training_data[i][0],vector_training_data[i][1])[2]=color1[2];

                        draw_large_labels(output_image, vector_training_data[i][0], vector_training_data[i][1], color1);
                    }

                    if(vector_training_data[i][2]==0) //NOT A BOUNDARY POINT/PIXEL
                    {
                        output_image(vector_training_data[i][0],vector_training_data[i][1])[0]=color0[0];
                        output_image(vector_training_data[i][0],vector_training_data[i][1])[1]=color0[1];
                        output_image(vector_training_data[i][0],vector_training_data[i][1])[2]=color0[2];

                        draw_large_labels(output_image, vector_training_data[i][0], vector_training_data[i][1], color0);
                    }

                    i++;
                }

                std::cout<<"Export labels to image: "<<filepath_output<<std::endl;
                exportImage(srcImageRange(output_image), vigra::ImageExportInfo(filepath_output.c_str()));
                    end_of_training=false;
            }

            //TESTS IF O-KEY (OTHER IMAGE) IS PRESSED
            if (main_disp.key==111||main_disp.key==79)
            {
                if (prob_large || !prob) prob_large=false;
                else prob_large=true;
                    clear=false;
                    end_of_training=false;
            }

            //Test if the mouse button is clicked on the image area
            if (main_disp.mouse_y>=0 && main_disp.mouse_x>=0)
            {
                //mouse stored to x_mouse and y_mouse to use it out of the loop
                x_mouse = main_disp.mouse_x;
                y_mouse = main_disp.mouse_y;
                if(main_disp.button||main_disp.key==97||main_disp.key==65)//mouse clicked or A-KEY (add) pressed
                {
                    if(clear==false)
                    {
                        add_training_data(vector_training_data,x_mouse+posx,y_mouse+posy,border);
                    }
                    if(clear==true) //"RUBBER MODE"
                    {
                        //REMOVE ALL TRAINING DATA IN A 3X3 NEIGHBOURHOOD
                        remove_training_data(vector_training_data,x_mouse+1+posx,y_mouse+1+posy);
                        remove_training_data(vector_training_data,x_mouse+posx,y_mouse+1+posy);
                        remove_training_data(vector_training_data,x_mouse-1+posx,y_mouse+1+posy);

                        remove_training_data(vector_training_data,x_mouse+1+posx,y_mouse+posy);
                        remove_training_data(vector_training_data,x_mouse+posx,y_mouse+posy);
                        remove_training_data(vector_training_data,x_mouse-1+posx,y_mouse+posy);

                        remove_training_data(vector_training_data,x_mouse+1+posx,y_mouse-1+posy);
                        remove_training_data(vector_training_data,x_mouse+posx,y_mouse-1+posy);
                        remove_training_data(vector_training_data,x_mouse-1+posx,y_mouse-1+posy);
                    }

                    //is_saved is used to see if the image must be saved or not
                    is_saved=false;
                    main_disp.display(img);
                       end_of_training=false;
                }
            }
        }
    }
}

int Pixel_training_image::get_nr_of_points(void)
{
    return vector_training_data.size();
}

int Pixel_training_image::get_nr_of_boundary_points(void)
{
    int i=0;
    int boundary_counter=0;
    while(i<(int)vector_training_data.size())
    {
        if(vector_training_data[i][2]==1)
        {
            boundary_counter++;
        }
        i++;
    }
    return boundary_counter;
}

int Pixel_training_image::get_nr_of_no_boundary_points(void)
{
    int i=0;
    int no_boundary_counter=0;
    while(i<(int)vector_training_data.size())
    {
        if(vector_training_data[i][2]==0)
        {
            no_boundary_counter++;
        }
        i++;
    }
    return no_boundary_counter;
}

std::vector<std::vector<int> >  Pixel_training_image::get_vector_training_data(void)
{
    return vector_training_data;
}

void Pixel_training_image::set_vector_training_data(std::vector<std::vector<int> > new_vector_training_data)
{
    vector_training_data=new_vector_training_data;
}

void import_and_interpolate_handmade_segmentatations(std::string segmentation_image_filepath ,std::string path_to_training_data)
{
    vigra::ImageImportInfo info(segmentation_image_filepath.c_str());

    int dim_x=info.width();
    int dim_y=info.height();

    //REZISIE IMAGE AND GRADIENT IMAGE TO THE SIZE
    vigra::BImage segmentation(dim_x,dim_y);
    importImage(info, destImage(segmentation));

    std::vector<std::vector<int> > vector_training_data;
    //FIND TRAINING DATA (dont search at the border)
    //AND INTERPOLATE IT
    for(int y=3;y<dim_y-3;y++)
    {
        for(int x=3;x<dim_x-3;x++)
        {
            if(segmentation(x,y)<=2 )
            {
                std::vector<int> one_training_point(3);
                one_training_point[0]=x;
                one_training_point[1]=y;
                one_training_point[2]=1;
                vector_training_data.push_back(one_training_point);
            }

        }
    }

    std::string filename_and_path_of_the_image=segmentation_image_filepath;
    std::string filename;
    std::string path_and_filename_to_save=path_to_training_data;
    filename=get_filename(filename_and_path_of_the_image);
    path_and_filename_to_save.append(filename);
    path_and_filename_to_save.append(".dat");

    std::ofstream training_file(path_and_filename_to_save.c_str());

    int i=0;
    while(i<(int)vector_training_data.size())
    {
        training_file <<vector_training_data[i][0]<<" "<<vector_training_data[i][1]<<" "<<vector_training_data[i][2]<< "\n"; //TODO STÖRT HIER DAS \n???
        i++;
    }

    training_file.close();
}

void Pixel_training_image::do_analysis(void)
{
    float grayvalue_histogram[256];
    for (int i=0;i<256;i++) grayvalue_histogram[i]=0;

    cimg_forXY(image_to_train_on, x, y)
    {
        grayvalue_histogram[image_to_train_on(x,y)]++;
    }

    std::string filename_and_path_of_the_image=filename_of_training_image;
    std::string filename;
    std::string path_and_filename_to_save=path_to_save_label_file;
    filename=get_filename(filename_and_path_of_the_image);
    path_and_filename_to_save.append(filename);
    path_and_filename_to_save.append(".dat");

    std::ofstream histogram_file(path_and_filename_to_save.c_str());

    for (int i=0;i<256;i++)
    {
        grayvalue_histogram[i]=grayvalue_histogram[i]/(float)(image_to_train_on.dimx()*image_to_train_on.dimy()/1000);

        float average=0.0f;
        for (int j=0;j<=i;j++) average+=(grayvalue_histogram[j]/(i+1));

        histogram_file <<i<<" "<<grayvalue_histogram[i]<<" "<<average<<" "<<grayvalue_histogram[i]/average<<"\n"; //TODO STÖRT HIER DAS \n???
    }

    histogram_file.close();
}

void Pixel_training_image::draw_large_labels(vigra::IRGBImage & output_image, int x, int y, const unsigned char * color)
{
    if (x>0)
    {
        output_image(x-1,y)[0]=color[0];
        output_image(x-1,y)[1]=color[1];
        output_image(x-1,y)[2]=color[2];
    }
    if (x>1)
    {
        output_image(x-2,y)[0]=color[0];
        output_image(x-2,y)[1]=color[1];
        output_image(x-2,y)[2]=color[2];
    }
    if (x+1<output_image.width())
    {
        output_image(x+1,y)[0]=color[0];
        output_image(x+1,y)[1]=color[1];
        output_image(x+1,y)[2]=color[2];
    }
    if (x+2<output_image.width())
    {
        output_image(x+2,y)[0]=color[0];
        output_image(x+2,y)[1]=color[1];
        output_image(x+2,y)[2]=color[2];
    }

    if (y>0)
    {  
        output_image(x,y-1)[0]=color[0];
        output_image(x,y-1)[1]=color[1];
        output_image(x,y-1)[2]=color[2];

        if (x>0)
        {
            output_image(x-1,y-1)[0]=color[0];
            output_image(x-1,y-1)[1]=color[1];
            output_image(x-1,y-1)[2]=color[2];
        }
        if (x>1)
        {
            output_image(x-2,y-1)[0]=color[0];
            output_image(x-2,y-1)[1]=color[1];
            output_image(x-2,y-1)[2]=color[2];
        }
        if (x+1<output_image.width())
        {
            output_image(x+1,y-1)[0]=color[0];
            output_image(x+1,y-1)[1]=color[1];
            output_image(x+1,y-1)[2]=color[2];
        }
        if (x+2<output_image.width())
        {
            output_image(x+2,y-1)[0]=color[0];
            output_image(x+2,y-1)[1]=color[1];
            output_image(x+2,y-1)[2]=color[2];
        }
    }

    if (y>1)
    {  
        output_image(x,y-2)[0]=color[0];
        output_image(x,y-2)[1]=color[1];
        output_image(x,y-2)[2]=color[2];

        if (x>0)
        {
            output_image(x-1,y-2)[0]=color[0];
            output_image(x-1,y-2)[1]=color[1];
            output_image(x-1,y-2)[2]=color[2];
        }
        if (x>1)
        {
            output_image(x-2,y-2)[0]=color[0];
            output_image(x-2,y-2)[1]=color[1];
            output_image(x-2,y-2)[2]=color[2];
        }
        if (x+1<output_image.width())
        {
            output_image(x+1,y-2)[0]=color[0];
            output_image(x+1,y-2)[1]=color[1];
            output_image(x+1,y-2)[2]=color[2];
        }
        if (x+2<output_image.width())
        {
            output_image(x+2,y-2)[0]=color[0];
            output_image(x+2,y-2)[1]=color[1];
            output_image(x+2,y-2)[2]=color[2];
        }
    }

    if (y+1<output_image.height())
    {  
        output_image(x,y+1)[0]=color[0];
        output_image(x,y+1)[1]=color[1];
        output_image(x,y+1)[2]=color[2];

        if (x>0)
        {
            output_image(x-1,y+1)[0]=color[0];
            output_image(x-1,y+1)[1]=color[1];
            output_image(x-1,y+1)[2]=color[2];
        }
        if (x>1)
        {
            output_image(x-2,y+1)[0]=color[0];
            output_image(x-2,y+1)[1]=color[1];
            output_image(x-2,y+1)[2]=color[2];
        }
        if (x+1<output_image.width())
        {
            output_image(x+1,y+1)[0]=color[0];
            output_image(x+1,y+1)[1]=color[1];
            output_image(x+1,y+1)[2]=color[2];
        }
        if (x+2<output_image.width())
        {
            output_image(x+2,y+1)[0]=color[0];
            output_image(x+2,y+1)[1]=color[1];
            output_image(x+2,y+1)[2]=color[2];
        }
    }

    if (y+2<output_image.height())
    {  
        output_image(x,y+2)[0]=color[0];
        output_image(x,y+2)[1]=color[1];
        output_image(x,y+2)[2]=color[2];

        if (x>0)
        {
            output_image(x-1,y+2)[0]=color[0];
            output_image(x-1,y+2)[1]=color[1];
            output_image(x-1,y+2)[2]=color[2];
        }
        if (x>1)
        {
            output_image(x-2,y+2)[0]=color[0];
            output_image(x-2,y+2)[1]=color[1];
            output_image(x-2,y+2)[2]=color[2];
        }
        if (x+1<output_image.width())
        {
            output_image(x+1,y+2)[0]=color[0];
            output_image(x+1,y+2)[1]=color[1];
            output_image(x+1,y+2)[2]=color[2];
        }
        if (x+2<output_image.width())
        {
            output_image(x+2,y+2)[0]=color[0];
            output_image(x+2,y+2)[1]=color[1];
            output_image(x+2,y+2)[2]=color[2];
        }
    }
}
