/*! \file boundary_data_structure.h
 * \brief Data structure for boundaries.
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

#include <vector>
#include <cgp/cgp_hdf5.hxx>
#include <vigra/convolution.hxx>

#define PI 3.14159

/*! \struct point
 * \brief A 2D point.
 */
struct point
{
    int x; /**< x-coordinate */
    int y; /**< y-coordinate */
};

/*! \fn calculatePhi(std::vector<point>* list, std::vector<float>* phi, int numSmoothing)
 * \brief Calculate from a list of points the angle of each point's normal.
 * \param list List of points
 * \param phi Vector the angle information is written into
 * \param numSmoothing No use?
 */
// Berechne aus einer parametrisierten Kurve (in form von Punkten in einer Liste) den Winkel
// der lokalen Normalen für jeden Kurvenpunkt
void calculatePhi(std::vector<point>* list, std::vector<float>* phi, int numSmoothing)
{
    int precision=2;//combine two pixels to one pixel

    if(list->size() < 2*precision) exit(-1); // necessary for further calculation

    //adapt to defined precision
    std::vector<point> new_list;
    for(int j=0; j<=list->size()-precision; j+=precision)
    {
        point av_point;
        av_point.x=0;
        av_point.y=0;

        for(int i=0; i<precision; i++)
        {
            av_point.x+=(*list)[j+i].x;
            av_point.y+=(*list)[j+i].y;
        }    

        av_point.x=av_point.x/precision;
        av_point.y=av_point.y/precision;

        new_list.push_back(av_point);
    }

    int size = new_list.size();

    if(size < 5)
    {
        //if the size of the reduced arc is smaller than 5, use a simpler approximation of the derivative

        float phi0 = atan2(new_list[1].y - new_list[0].y, new_list[1].x - new_list[0].x );
        (*phi).push_back( phi0 + PI/2.0f );

        for(int j=1; j<new_list.size()-1; j++)
        {
            float phi_temp = atan2(new_list[j+1].y - new_list[j-1].y, new_list[j+1].x - new_list[j-1].x );
            (*phi).push_back( phi_temp + PI/2.0f );
        }

        float phi_end = atan2(new_list.back().y - new_list[new_list.size()-2].y, new_list.back().x - new_list[new_list.size()-2].x );
        (*phi).push_back( phi_end + PI/2.0f );
    }
    else
    {
        (*phi).push_back(0.0f);
        (*phi).push_back(0.0f);

        std::vector<float> dx(size);
        std::vector<float> dy(size);

        // Calculate dx and dy using 5 pixels approximation
        for(int j=2; j<size-2; j++)
        {
            dx[j] = ( new_list[j+2].x - new_list[j-2].x + 8*new_list[j+1].x - 8*new_list[j-1].x )/12.0f;
            dy[j] = ( new_list[j+2].y - new_list[j-2].y + 8*new_list[j+1].y - 8*new_list[j-1].y )/12.0f;
        }

        dx[0]=dx[2];
        dx[1]=dx[2];
        dx[size-2]=dx[size-3];
        dx[size-1]=dx[size-3];

        dy[0]=dy[2];
        dy[1]=dy[2];
        dy[size-2]=dy[size-3];
        dy[size-1]=dy[size-3];

        // Smoothing dx and dy: Convolution with binomial mask
        std::vector<float> temp(size);
        temp[0]=dx[0];
        temp[1]=dx[1];
        temp[size-2]=dx[size-2];
        temp[size-1]=dx[size-1];

        for(int a=0; a<10; a++)
        {
            for(int j=2; j<size-2; j++)
                temp[j] = ( dx[j+2]+ dx[j+1]*4.0f+ dx[j]*6.0f+ dx[j-1]*4.0f+ dx[j-2] ) / 16.0f;

            for(int j=2; j<size-2; j++) 
                dx[j] = ( temp[j+2]+ temp[j+1]*4.0f+ temp[j]*6.0f+ temp[j-1]*4.0f+ temp[j-2] ) / 16.0f;
        }

        temp[0]=dy[0];
        temp[1]=dy[1];
        temp[size-2]=dy[size-2];
        temp[size-1]=dy[size-1];

        for(int a=0; a<10; a++)
        {
            for(int j=2; j<size-2; j++)
                temp[j] = ( dy[j+2]+ dy[j+1]*4.0f+ dy[j]*6.0f+ dy[j-1]*4.0f+ dy[j-2] ) / 16.0f;

            for(int j=2; j<size-2; j++) 
                dy[j] = ( temp[j+2]+ temp[j+1]*4.0f+ temp[j]*6.0f+ temp[j-1]*4.0f+ temp[j-2] ) / 16.0f;
        }

        //calculate orientation angle phi
        for(int j=2; j<size-2; j++)
        {
            float phi_temp=atan2(dy[j],dx[j]);
            (*phi).push_back(phi_temp);
            //(*phi).push_back(dx[j]);
        }

        float pback = (*phi).back();

        (*phi).push_back(pback);
        (*phi).push_back(pback);

        (*phi)[0] = (*phi)[2];
        (*phi)[1] = (*phi)[2];

        std::vector<float> old;

        //look for 2 pi flips
        for(int i=0; i<size; i++)
        {
            bool not_found=true;

            //compare phi with (up to 5) previous values
            for(int m=0; m<old.size() && not_found; m++)
            {
                if ( fabs(old[m] - (*phi)[i]) > fabs(old[m] - (*phi)[i] + 2.0f*PI) )
                {
                    not_found=false;
                    (*phi)[i]-=2.0f*PI;
                }
                else if ( fabs(old[m] - (*phi)[i]) > fabs(old[m] - (*phi)[i] - 2.0f*PI) )
                {
                    not_found=false;
                    (*phi)[i]+=2.0f*PI;
                }
            }

            old.push_back((*phi)[i]);
            if (old.size()>5) old.erase(old.begin());
        }

        (*phi)[0] = (*phi)[2];
        (*phi)[1] = (*phi)[2];

        (*phi)[size-1] = (*phi)[size-3];
        (*phi)[size-2] = (*phi)[size-3];

        // Smoothing orientation angle phi: Convolution with binomial mask
        temp[0]=(*phi)[0];
        temp[1]=(*phi)[1];
        temp[size-2]=(*phi)[size-2];
        temp[size-1]=(*phi)[size-1];

        for(int a=0; a<10; a++)
        {
            for(int j=2; j<size-2; j++)
                temp[j] = ( (*phi)[j+2]+ (*phi)[j+1]*4.0f+ (*phi)[j]*6.0f+ (*phi)[j-1]*4.0f+ (*phi)[j-2] ) / 16.0f;

            for(int j=2; j<size-2; j++) 
                (*phi)[j] = ( temp[j+2]+ temp[j+1]*4.0f+ temp[j]*6.0f+ temp[j-1]*4.0f+ temp[j-2] ) / 16.0f;
        }

        //now resize to original arc length
        int j=0;
        while(phi->size()<2*size-1)
        {
            float mean=((*phi)[j]+(*phi)[j+1])/2.0f;
            j++;
            phi->insert(phi->begin()+j, mean);
            j++;
        }

        pback = (*phi).back();

        while(phi->size()<list->size()) (*phi).push_back(pback);
    }
}

/*! \fn calculateCurv( std::vector<float>* phi, std::vector<float>* curv )
 * \brief Calculate the local curvature from the local normals.
 * \param phi Vector containing the angle information
 * \param curv Vector the local curvature is written into
 */
// calculate the local curvature from the local normals
void calculateCurv( std::vector<float>* phi, std::vector<float>* curv )
{
    curv->push_back( 0.0f );

    for(int i = 1; i < phi->size()-1; i++)
    {
        float temp=((*phi)[i+1] - (*phi)[i-1])/2.0f;
        if (temp<-PI) temp+=PI*2.0f;
        if (temp>PI) temp-=PI*2.0f;
        curv->push_back(fabs(temp));
    }

    float cback = (*curv).back();

    curv->push_back(cback);
    (*curv)[0] = (*curv)[1];

    int size=curv->size();

    // Smoothing curvature: Convolution with binomial mask
    std::vector<float> temp(size);
    temp[0]=(*curv)[0];
    temp[1]=(*curv)[1];
    temp[size-2]=(*curv)[size-2];
    temp[size-1]=(*curv)[size-1];

    for(int a=0; a<10; a++)
    {
        for(int j=2; j<size-2; j++)
            temp[j] = ( (*curv)[j+2]+ (*curv)[j+1]*4.0f+ (*curv)[j]*6.0f+ (*curv)[j-1]*4.0f+ (*curv)[j-2] ) / 16.0f;

        for(int j=2; j<size-2; j++) 
            (*curv)[j] = ( temp[j+2]+ temp[j+1]*4.0f+ temp[j]*6.0f+ temp[j-1]*4.0f+ temp[j-2] ) / 16.0f;
    }
}

/*! \fn calculateCrossSection( std::vector<point>* arc, std::vector<float>* phi, vigra::BImage* image, std::vector<std::vector<float> >* cross )
 * \brief Calculate the cross-section pixels for one arc.
 * \param arc Vector containing the arcs
 * \param phi Vector containing the angle information
 * \param image Image providing width and height information
 * \param cross Vector the cross section information is written into
 */
// calculate the cross-section pixels for one arc
void calculateCrossSection( std::vector<point>* arc, std::vector<float>* phi, vigra::BImage* image, std::vector<std::vector<float> >* cross )
{
    // loop through the pixels of the arc
    for(int i=0; i<(*arc).size(); i++)
    {
        // sample the cross section using information on the local normal
        std::vector<float> crossSection;
        for(int j=-10; j<11; j++)
        {
            int x =  (int)((*arc)[i].x + cos((*phi)[i])*j);
            int y =  (int)((*arc)[i].y + sin((*phi)[i])*j);

            if( 0 <= x && (*image).width() > x + 1 && 0 <= y && (*image).height() > y + 1 )
                crossSection.push_back((float)(*image)(x,y));
        }	
        (*cross).push_back(crossSection);
    }
}

/*! \fn sort_two_cell(std::vector<point> & this_arc, std::vector<point> & this_arc_sorted,int x_center,int y_center)
 * \brief Sort two cells.
 */
void sort_two_cell(std::vector<point> & this_arc, std::vector<point> & this_arc_sorted,int x_center,int y_center)
{
    //first we search the point next to the center point
    size_t x_start=x_center;
    size_t y_start=y_center;
    bool found_start_arc=false;
    bool start_arc_front=true;
    point start_arc;
    for(size_t k=0; k<this_arc.size() && found_start_arc==false; ++k)
    {
        //search starting arc
        size_t xx=this_arc[k].x;
        size_t yy=this_arc[k].y;

        //left
        if(yy==y_start && xx==(x_start-1) && found_start_arc==false)
        {
            found_start_arc=true;
        }
        //right
        if(yy==y_start && xx==(x_start+1) && found_start_arc==false)
        {
            found_start_arc=true;
        }
        //up
        if(yy==(y_start+1) && xx==x_start && found_start_arc==false)
        {
            found_start_arc=true;
        }
        //down
        if(yy==(y_start-1) && xx==x_start && found_start_arc==false)
        {
            found_start_arc=true;
        }

        if(found_start_arc==true)
        {
            start_arc.x=xx;
            start_arc.y=yy;
            this_arc_sorted.push_back(start_arc);
            if (k>this_arc.size()/2) start_arc_front=false;
        }
    }

    size_t nr_of_arc_elements_found=0;
    if(found_start_arc==true)
    {
        nr_of_arc_elements_found=1;

        point last_arc_found=this_arc_sorted[0];
        size_t xx_last=last_arc_found.x;
        size_t yy_last=last_arc_found.y;

        //now he have to search for the next point
        //for schleife die solange läuft bis tatsächlich alle gefunden sind
        //bool v_arc=false;
        //bool h_arc=false;

        while(nr_of_arc_elements_found!=this_arc.size())
        {
            //std::cout<<nr_of_arc_elements_found<<" of "<<this_arc.size()<<std::endl;
            //each arc has 6 possible neigbours
            //but first we must know ,if the last arc which has been found , is horiontal or verical

            bool found_next=false;
            point next_arc;

            bool visited_vertical_or_horizontal=false;
            //vertical arc
            if(xx_last%2!=0 && yy_last%2==0 && visited_vertical_or_horizontal==false)
            {
                visited_vertical_or_horizontal=true;

                size_t k_sort;
                if (start_arc_front) k_sort=0;
                else k_sort=this_arc.size()-1;

                while(k_sort<this_arc.size() && k_sort>=0 && found_next==false)
                {
                    //search starting arc
                    size_t xx=this_arc[k_sort].x;
                    size_t yy=this_arc[k_sort].y;

                    //we have to make sure that we dont add the last arc agaiN!!
                    bool is_the_last_one=false;
                    found_next=false;

                    if(xx==last_arc_found.x && yy==last_arc_found.y)
                    {
                        is_the_last_one=true;
                    }

                    //now we have to search wich of the 6 possible neigbours is the real one

                    //up left
                    if(xx==xx_last-1 && yy==yy_last+1 && found_next==false && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //up
                    if(xx==xx_last && yy==yy_last+2 && found_next==false && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //up right
                    if(xx==xx_last+1 && yy==yy_last+1 && found_next==false && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //down right
                    if(xx==xx_last+1 && yy==yy_last-1 && found_next==false && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //down
                    if(xx==xx_last && yy==yy_last-2 && found_next==false && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //down left
                    if(xx==xx_last-1 && yy==yy_last-1 && found_next==false && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    if(found_next==true)
                    {
                        this_arc_sorted.push_back(next_arc);
                        last_arc_found.x=this_arc_sorted[this_arc_sorted.size()-2].x;
                        last_arc_found.y=this_arc_sorted[this_arc_sorted.size()-2].y;
                        xx_last=next_arc.x;
                        yy_last=next_arc.y;

                        nr_of_arc_elements_found++;
                        this_arc[k_sort].x=-10;
                    }

                    if (start_arc_front) ++k_sort;
                    else --k_sort;
                }
            }

            //horizontal
            if(xx_last%2==0 && yy_last%2!=0 && visited_vertical_or_horizontal==false)
            {
                visited_vertical_or_horizontal=true;
                //std::cout<<"horizontal"<<std::endl;

                size_t k_sort;
                if (start_arc_front) k_sort=0;
                else k_sort=this_arc.size()-1;

                while(k_sort<this_arc.size() && k_sort>=0 && found_next==false)
                {
                    //search starting arc
                    size_t xx=this_arc[k_sort].x;
                    size_t yy=this_arc[k_sort].y;

                    //we have to make sure that we dont add the last arc agaiN!!
                    found_next=false;
                    bool is_the_last_one=false;
                    if(xx==last_arc_found.x && yy==last_arc_found.y )
                    {
                        is_the_last_one=true;
                        //std::cout<<"found last one"<<std::endl;
                    }
                    //now we have to search wich of the 6 possible neigbours is the real one

                    //up left
                    if(xx==xx_last-1 && yy==yy_last+1 && found_next==false  && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //left
                    if(xx==xx_last-2 && yy==yy_last && found_next==false  && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //up right
                    if(xx==xx_last+1 && yy==yy_last+1 && found_next==false  && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //down right
                    if(xx==xx_last+1 && yy==yy_last-1 && found_next==false  && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //right
                    if(xx==xx_last+2 && yy==yy_last && found_next==false  && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    //down left
                    if(xx==xx_last-1 && yy==yy_last-1 && found_next==false  && is_the_last_one==false)
                    {
                        found_next=true;
                        next_arc.x=xx;
                        next_arc.y=yy;
                    }

                    if(found_next==true)
                    {
                        this_arc_sorted.push_back(next_arc);
                        last_arc_found.x=this_arc_sorted[this_arc_sorted.size()-2].x;
                        last_arc_found.y=this_arc_sorted[this_arc_sorted.size()-2].y;
                        xx_last=next_arc.x;
                        yy_last=next_arc.y;

                        nr_of_arc_elements_found++;
                        this_arc[k_sort].x=-10;
 
                    }

                    if (start_arc_front) ++k_sort;
                    else --k_sort;
                }
            }

            if(visited_vertical_or_horizontal==false)
            {
                std::cout<<"error in sort"<<std::endl;
            }
        }

        //NOW WE SHOULD HAVE A SORTED VECTOR OF THE ARC!

        //DEBUG CHECK IF IT WORKED FINE

        if(this_arc_sorted.size()==this_arc.size())
        {
            //std::cout<<"size is the same"<<std::endl;
        }
        else
        {
            std::cout<<"size is  NOT the same"<<std::endl;
        }

    }
    else
    {
        std::cout<<"could not find the begin of the arc"<<std::endl;
    }
}
