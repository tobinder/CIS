/*! \file cgp_structure.hxx
 * \brief CGP data structure.
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

#include "boundary_data_structure.h"
#include "marray.hxx"
#include <sys/time.h>
#include <omp.h>

#if HAVE_MPI
#include <cgp/CgpxMaster.hxx>
#include <cgp/CgpxWorker.hxx>
#endif

typedef unsigned int label_type;
typedef int coordinate_type;

//3-cell is a region of a image
//2-cell is a boundary of an image (where 2 regions meet/touch ech other)
//1-cell is a "junction" of a image (where 3 or 4 boundarys meet)
//0-cell are not existend in a regular 2D-Image! (because 2 junctions never "cross")

/*! \fn compute_and_save_cgp_data_structure(std::string dest_path, int size, int rank)
 * \brief Compute and save the CGP data structure.
 * \param dest_path Destination path
 * \param size Size of the CGP structure
 * \param rank Rank of the CGP structure
 */
void compute_and_save_cgp_data_structure(std::string dest_path, int size, int rank)
{
    std::string segFilename = dest_path.c_str();
    segFilename.append(".h5");

    std::string datasetName = "ws_image";
    
    std::string topologyFilename = dest_path.c_str();
    topologyFilename.append(".grid.h5");

	std::string geometryFilename = dest_path.c_str();
    geometryFilename.append(".objects.h5");
    
    // open segmentation file and dataset
    hid_t segFileHandle = H5Fopen(segFilename.c_str(),
        H5F_ACC_RDONLY, H5P_DEFAULT);
    if(segFileHandle < 0) {
        std::cerr << "Error opening: \"" << segFilename << "\"" << std::endl;
        exit(-1);
    }

    // Open an existing dataset
    hid_t dataset_id = H5Dopen(segFileHandle, datasetName.c_str(), H5P_DEFAULT);

    // Get filespace handle
    hid_t filespace = H5Dget_space(dataset_id);

    // Get dimensions
    hsize_t dims[3];
    H5Sget_simple_extent_dims(filespace, dims, NULL);

    int dim_y=dims[0];
    int dim_x=dims[1];

    // Close the filespace
    H5Sclose(filespace);

    // Close the dataset
    H5Dclose(dataset_id);

    marray::Vector<coordinate_type> blockShape(3);
    blockShape(0) = 2000;
    blockShape(1) = 2000;
    blockShape(2) = 1;

    #if (!HAVE_MPI)

        typedef cgp::hdf5::TopologyBlockwiseMaster<label_type, coordinate_type> Master;
        typedef cgp::hdf5::TopologyBlockwiseWorker<label_type, coordinate_type> Worker;
        Master master(blockShape, segFileHandle, datasetName, topologyFilename, true);
        Worker worker(blockShape, segFileHandle, datasetName, topologyFilename, true);

        while(master.pending()) {
            int j = master.next();
            worker.doJob(j);
        }
        master.postprocess();
        
    #else /* MPI available */

        if (size==1)//started as one process
        {
            typedef cgp::hdf5::TopologyBlockwiseMaster<label_type, coordinate_type> Master;
            typedef cgp::hdf5::TopologyBlockwiseWorker<label_type, coordinate_type> Worker;
            Master master(blockShape, segFileHandle, datasetName, topologyFilename, true);
            Worker worker(blockShape, segFileHandle, datasetName, topologyFilename, true);

            while(master.pending()) {
                int j = master.next();
                worker.doJob(j);
            }
            master.postprocess();
        }
        else
        {
            typedef CgpxMaster<label_type, coordinate_type> Master;
            typedef CgpxWorker<label_type, coordinate_type> Worker;

            if(rank==0)
            {
                Master master(size-1 /*number workers*/, MPI_COMM_WORLD,
                              blockShape, segFileHandle, datasetName,
                              topologyFilename);
                master.run();
            }
            else
            {
                Worker worker(blockShape, segFileHandle, datasetName,
                              topologyFilename);
                worker.run();
            }
            MPI_Finalize();
        }

    #endif // HAVE_MPI

    if (rank==0)
    {
        cgp::hdf5::GeometrySuppression s= cgp::hdf5::SUPPRESS_3_SETS;

	    size_t numberOfBinningGroups = 1000;

	    // open topology file
	    hid_t topologyFileHandle = H5Fopen(topologyFilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	    if(topologyFileHandle < 0)
        {
		    std::cout << "Error opening: " << topologyFilename << std::endl;
		    exit(-1);
	    }

	    // create geometry file
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
	    hid_t geometryFileHandle = H5Fcreate(geometryFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
	    if(geometryFileHandle < 0)
        {
		    std::cout << "Error creating file: " << geometryFilename << std::endl;
		    H5Fclose(topologyFileHandle);		
		    exit(-1);
	    }

	    // extract topology and geometry
	    cgp::hdf5::geometry3blockwise<label_type, coordinate_type>(topologyFileHandle, geometryFileHandle, numberOfBinningGroups, s, true);
	
	    // clean-up
	    H5Fclose(topologyFileHandle);		
        remove(topologyFilename.c_str());

        std::cout<<"sort 2-sets..."<<std::endl;

        //open needed structures from file to sort two cells
        marray::Marray<unsigned int> one_boundings;
        std::vector< std::vector<point> > arcs;
        std::vector<point> junctions;        

        htri_t check = H5Lexists(geometryFileHandle, "/neighborhood-1", H5P_DEFAULT);

        if(check)
        {
            // Open an existing dataset
            hid_t dataset_id = H5Dopen(geometryFileHandle, "/neighborhood-1", H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            size_t size[] = {dims[1],dims[0]};
            one_boundings.resize(size,size+2);
            
            //loop over columns in one_boundings
            for(int i= 0; i<dims[1]; i++)
            {
                // Datastructure to read in
                unsigned int * column_values = new unsigned int[dims[1]];

                // dataspace for one column
                hsize_t column[2];
                column[0] = dims[0];
                column[1] = 1;  
                hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

                // select file hyperslab 
                hsize_t start_h[2];// start of hyperslab
                hsize_t count[2];// block count
                
                count[0]  = dims[0]; 
                count[1]  = 1;

                start_h[0]  = 0; 
                start_h[1]  = i;

                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                H5Dread(dataset_id, H5T_NATIVE_UINT, mdataspace_id, filespace, H5P_DEFAULT, column_values);

                for (int x=0;x<dims[0];x++)
                {            
                    unsigned int value=(unsigned int)column_values[x];
                    one_boundings(i,x)=value;
                }
                // Close the memoryspace
                H5Sclose(mdataspace_id);

                delete column_values;
            }

            // Close the filespace
            H5Sclose(filespace);

            // Close the dataset
            H5Dclose(dataset_id);
        }
        /*
        //KURZE AUSGABE ZUM TESTEN der 1 boundings
        for(int y=0;y<(int)one_boundings.shape(0);y++)
        {
            for(int x=0;x<(int)one_boundings.shape(1);x++)
            {
                std::cout<<one_boundings(y,x)<<"   ";
            }
            std::cout<<std::endl;
        }
        */
        else
        {
            std::cout << "Segmentation has no junctions!" << std::endl;
            H5Fclose(geometryFileHandle);
            exit(-1);
        }

        unsigned int max_labels[4];

        {
            // Open an existing dataset
            hid_t dataset_id = H5Dopen(geometryFileHandle, "/max-labels", H5P_DEFAULT);

            H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, max_labels);

            // Close the dataset
            H5Dclose(dataset_id);
        }

        //std::cout<<max_labels[0]<<" "<<max_labels[1]<<" "<<max_labels[2]<<" "<<max_labels[3]<<std::endl;

        // Open an existing group
        hid_t one_sets_group_id = H5Gopen(geometryFileHandle, "/1-sets", H5P_DEFAULT);

        for (size_t i=1; i<max_labels[1]+1; i++)
        {
            std::stringstream datasetNameStream;
            datasetNameStream << i << "-1";
            std::stringstream binningGroupName;
            binningGroupName << "bin-" << i % numberOfBinningGroups;
            hid_t h = H5Gopen(one_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);

            // Open an existing dataset
            hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);

            //Datastructure to read in
            int junction_pos[3];

            H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, junction_pos);

            point junction_p;
            junction_p.x=junction_pos[1];
            junction_p.y=junction_pos[0];

            junctions.push_back(junction_p);

            // Close the dataset
            H5Dclose(dataset_id);

            H5Gclose(h);
        }

        H5Gclose(one_sets_group_id);

        long unsigned int * arc_part_counters= new long unsigned int [max_labels[2]];

        {
            // Open an existing dataset
            hid_t dataset_id = H5Dopen(geometryFileHandle, "/parts-counters-2", H5P_DEFAULT);

            H5Dread(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, arc_part_counters);

            // Close the dataset
            H5Dclose(dataset_id);
        }

        // Open an existing group
        hid_t two_sets_group_id = H5Gopen(geometryFileHandle, "/2-sets", H5P_DEFAULT);

        arcs.resize(max_labels[2]);

        /*
        int arc_progress=0;
        int last_printout=0;

        double temp = max_labels[2]/100;
        round(temp);
        int arc_number_rounded = temp*100;
        */

        timeval start, end;
        gettimeofday(&start, 0);

        #pragma omp parallel
        {
            #pragma omp for
            for (int i=1; i<max_labels[2]+1; i++)
            {
                std::vector<point> this_arc;
                std::vector<point> this_arc_sorted;

                std::stringstream binningGroupName;
                binningGroupName << "bin-" << i % numberOfBinningGroups;

                #pragma omp critical
                {
                    //open bin-group
                    hid_t h = H5Gopen(two_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);

                    for (size_t part=1; part<arc_part_counters[i-1]+1; part++)
                    {
                        std::stringstream datasetNameStream;
                        datasetNameStream << i <<"-"<<part;//check for -2 and -3 ... and combine arc pixels
                        
                        // Open an existing dataset
                        hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);

                        // Get filespace handle
                        hid_t filespace = H5Dget_space(dataset_id);

                        // Get dimensions
                        hsize_t dims[2];
                        H5Sget_simple_extent_dims(filespace, dims, NULL);

                        // Datastructure to read in
                        int * data_out_x = new int[dims[1]];
                        int * data_out_y = new int[dims[1]];

                        //dataspace for one coordinate
                        hsize_t row[2];
                        row[0] = 1;
                        row[1] = dims[1];  
                        hid_t mdataspace_id = H5Screate_simple(2, row, NULL);

                        // select file hyperslab 
                        hsize_t start_h[2];// start of hyperslab
                        hsize_t count[2];// block count
                        
                        count[0]  = 1; 
                        count[1]  = dims[1];

                        start_h[0]  = 0; 
                        start_h[1]  = 0;

                        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                        H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_y);

                        start_h[0]  = 1;

                        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                        H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_x);

                        for (int j=0;j<dims[1];j++)
                        {
                            point push;
                            push.x=abs(data_out_x[j]);
                            push.y=abs(data_out_y[j]);
                            this_arc.push_back(push);
                        }

                        // Close the memoryspace
                        H5Sclose(mdataspace_id);

                        // Close the filespace
                        H5Sclose(filespace);

                        // Close the dataset
                        H5Dclose(dataset_id);

                        delete data_out_x;
                        delete data_out_y;

                    }

                    //close bin-group
                    H5Gclose(h);
                }
                #pragma omp end critical

                int x_center=-1, y_center=-1;

                for(int y=0;y<(int)one_boundings.shape(0);y++)
                {
                    for(int x=0;x<(int)one_boundings.shape(1);x++)
                    {
                        if (one_boundings(y,x)==i)
                        {
                            x_center=junctions[y].x;
                            y_center=junctions[y].y;
                        }
                    }
                }

                //arc has no junction, start point for sorting is needed
                if (x_center<0)
                {
                    bool found_border_pixel=false;
    
                    //if there is connectivity with image border, sorting is not relevant, area will belong to outside label
                    for (int p=0; p<this_arc.size() && !found_border_pixel; p++)
                    {
                        if (this_arc[p].x==0 || this_arc[p].x/2==dim_x-1 || this_arc[p].y==0 || this_arc[p].y/2==dim_y-1)
                        {
                            this_arc_sorted=this_arc;
                            found_border_pixel=true;
                        }
                    }

                    //take first arc point as start point
                    if(!found_border_pixel)
                    {
                        x_center=this_arc[0].x+1;
                        y_center=this_arc[0].y;
                        sort_two_cell(this_arc,this_arc_sorted,x_center,y_center);
                    }
                }
                else sort_two_cell(this_arc,this_arc_sorted,x_center,y_center);

                //now write sorted arc back

                #pragma omp critical
                {
                    //open bin-group
                    hid_t h = H5Gopen(two_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);

                    for (size_t part=1; part<arc_part_counters[i-1]+1; part++)
                    {
                        std::stringstream datasetNameStream;
                        datasetNameStream << i <<"-"<<part;//check for -2 and -3 ... and combine arc pixels
                        
                        // Open an existing dataset

                        hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);

                        // Get filespace handle
                        hid_t filespace = H5Dget_space(dataset_id);

                        // Get dimensions
                        hsize_t dims[2];
                        H5Sget_simple_extent_dims(filespace, dims, NULL);

                        // Datastructure to write
                        int * data_out_x = new int[dims[1]];
                        int * data_out_y = new int[dims[1]];

                        for (int j=0;j<dims[1];j++)
                        {
                            data_out_x[j]=this_arc_sorted[0].x;
                            data_out_y[j]=this_arc_sorted[0].y;
                            this_arc_sorted.erase(this_arc_sorted.begin());
                        }
                    
                        //dataspace for one coordinate
                        hsize_t row[2];
                        row[0] = 1;
                        row[1] = dims[1];  
                        hid_t mdataspace_id = H5Screate_simple(2, row, NULL);

                        // select file hyperslab 
                        hsize_t start_h[2];// start of hyperslab
                        hsize_t count[2];// block count
                        
                        count[0]  = 1; 
                        count[1]  = dims[1];

                        start_h[0]  = 0; 
                        start_h[1]  = 0;

                        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                        H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_y);

                        start_h[0]  = 1;

                        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                        H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_x);

                        // Close the memoryspace
                        H5Sclose(mdataspace_id);

                        // Close the filespace
                        H5Sclose(filespace);

                        // Close the dataset
                        H5Dclose(dataset_id);

                        delete data_out_x;
                        delete data_out_y;

                    }

                    //close bin-group
                    H5Gclose(h);
                }

                /*
                //nr of arcs processed
                arc_progress++;

                //a percent step is reached
                if ( ((100*arc_progress)%arc_number_rounded) == 0.0f )
                {
                    //progress in percent
                    int progress_ratio=100*arc_progress/max_labels[2]+1;

                    if (progress_ratio!=last_printout)
                    {
                        if (progress_ratio==1)
                        {
                            gettimeofday(&end, 0);
                            int min=0;

                            //multiplied by 100 to get total runtime
                            int sec=100*(end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/10000;
                            while (sec>59)
                            {
                                min=min+1;
                                sec=sec-60;
                            }

                            std::cout << "Total runtime approximately: "<<min<<" min "<<sec<<" s"<< std::endl;
                        }

                        if ((progress_ratio)%10 > 0) std::cout<<".";
                        else std::cout<<progress_ratio<<"% done ("<<arc_progress<<"/"<<max_labels[2]<<")"<<std::endl;
                        last_printout=progress_ratio;
                    }

                }
                */

                #pragma omp end critical

            }

        }//end of parallelisation

        delete arc_part_counters;

        H5Gclose(two_sets_group_id);

        // Close the file
        H5Fclose(geometryFileHandle);
    }

    std::cout<<"...done"<<std::endl;

}

/*! \fn load_cgp_data_structure(std::string filepath_to_ws_region_image,
                             vigra::BasicImage<unsigned int> & ws_region_image,
                             marray::Marray<unsigned int> & one_boundings,
                             marray::Marray<unsigned int> & two_boundings,
                             std::vector< std::vector<point> > & arcs,
                             std::vector<point> & junctions,
                             int & dim_x,
                             int & dim_y,
                             bool two_pixel_boundary)
 * \brief Load the CGP data structure into the data structures given as arguments.
 * \param filepath_to_ws_region_image Filepath to the watershed segmentation image.
 * \param ws_region_image Watershed segmentation image
 * \param one_boundings Marray representing one boundings
 * \param two_boundings Marray representing two boundings
 * \param arcs Vector containing the segmentation arcs
 * \param junctions Vector containing the segmentation junctions
 * \param dim_x Image width
 * \param dim_y Image height
 * \param two_pixel_boundary Flag to activate two pixel wide boundaries
 */
void load_cgp_data_structure(std::string filepath_to_ws_region_image,
                             vigra::BasicImage<unsigned int> & ws_region_image,
                             marray::Marray<unsigned int> & one_boundings,
                             marray::Marray<unsigned int> & two_boundings,
                             std::vector< std::vector<point> > & arcs,
                             std::vector<point> & junctions,
                             int & dim_x,
                             int & dim_y,
                             bool two_pixel_boundary)
{
    const size_t numberOfBinningGroups = 1000;

    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_to_ws_region_image<<std::endl;

    // Open an existing file
    hid_t file_id0 = H5Fopen(filepath_to_ws_region_image.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_id0, "/ws_image", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[3];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        dim_y=dims[0];
        dim_x=dims[1];

        ws_region_image.resize(dims[1],dims[0]);

        //loop over rows in ws_image
        for(int i= 0; i<dims[0]; i++)
        {
            // Datastructure to read in
            unsigned int * row_values = new unsigned int[dims[1]];

            // dataspace for one row
            hsize_t row[3];
            row[0] = 1;
            row[1] = dims[1];
            row[2] = 1;  
            hid_t mdataspace_id = H5Screate_simple(3, row, NULL);

            // select file hyperslab 
            hsize_t start_h[3];// start of hyperslab
            hsize_t count[3];// block count
            
            count[0]  = 1; 
            count[1]  = dims[1];
            count[2]  = 1; 

            start_h[0]  = i; 
            start_h[1]  = 0;
            start_h[2]  = 0;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_UINT, mdataspace_id, filespace, H5P_DEFAULT, row_values);

            for (int x=0;x<dims[1];x++)
            {            
                unsigned int value=row_values[x];
                ws_region_image(x,i)=value;
            }
            // Close the memoryspace
            H5Sclose(mdataspace_id);

            delete row_values;
        }

        // Close the filespace
        H5Sclose(filespace);

        // Close the dataset
        H5Dclose(dataset_id);
    }

    // Close the file
    H5Fclose(file_id0);

    std::cout<<"...done"<<std::endl;

    filepath_to_ws_region_image.resize(filepath_to_ws_region_image.size()-3);
    filepath_to_ws_region_image.append(".objects.h5");

    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_to_ws_region_image<<std::endl;

    // Open an existing file
    hid_t file_id = H5Fopen(filepath_to_ws_region_image.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_id, "/neighborhood-2", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        size_t size[] = {dims[1],dims[0]};
        two_boundings.resize(size,size+2);

        //loop over columns in two_boundings
        for(int i= 0; i<dims[1]; i++)
        {
            // Datastructure to read in
            unsigned int * column_values = new unsigned int[dims[1]];

            // dataspace for one column
            hsize_t column[2];
            column[0] = dims[0];
            column[1] = 1;  
            hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

            // select file hyperslab 
            hsize_t start_h[2];// start of hyperslab
            hsize_t count[2];// block count
            
            count[0]  = dims[0]; 
            count[1]  = 1;

            start_h[0]  = 0; 
            start_h[1]  = i;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_UINT, mdataspace_id, filespace, H5P_DEFAULT, column_values);

            for (int x=0;x<dims[0];x++)
            {            
                unsigned int value=(unsigned int)column_values[x];
                two_boundings(i,x)=value;
            }
            // Close the memoryspace
            H5Sclose(mdataspace_id);

            delete column_values;
        }

        // Close the filespace
        H5Sclose(filespace);

        // Close the dataset
        H5Dclose(dataset_id);
    }
    /*
    //KURZE AUSGABE ZUM TESTEN der 2 boundings
    for(int y=0;y<(int)two_boundings.shape(0);y++)
    {
        for(int x=0;x<(int)two_boundings.shape(1);x++)
        {
            std::cout<<two_boundings(y,x)<<"   ";
        }
        std::cout<<std::endl;
    }
    */
    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_id, "/neighborhood-1", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        size_t size[] = {dims[1],dims[0]};
        one_boundings.resize(size,size+2);

        //loop over columns in one_boundings
        for(int i= 0; i<dims[1]; i++)
        {
            // Datastructure to read in
            unsigned int * column_values = new unsigned int[dims[1]];

            // dataspace for one column
            hsize_t column[2];
            column[0] = dims[0];
            column[1] = 1;  
            hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

            // select file hyperslab 
            hsize_t start_h[2];// start of hyperslab
            hsize_t count[2];// block count
            
            count[0]  = dims[0]; 
            count[1]  = 1;

            start_h[0]  = 0; 
            start_h[1]  = i;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_UINT, mdataspace_id, filespace, H5P_DEFAULT, column_values);

            for (int x=0;x<dims[0];x++)
            {            
                unsigned int value=(unsigned int)column_values[x];
                one_boundings(i,x)=value;
            }
            // Close the memoryspace
            H5Sclose(mdataspace_id);

            delete column_values;
        }

        // Close the filespace
        H5Sclose(filespace);

        // Close the dataset
        H5Dclose(dataset_id);
    }
    /*
    //KURZE AUSGABE ZUM TESTEN der 1 boundings
    for(int y=0;y<(int)one_boundings.shape(0);y++)
    {
        for(int x=0;x<(int)one_boundings.shape(1);x++)
        {
            std::cout<<one_boundings(y,x)<<"   ";
        }
        std::cout<<std::endl;
    }
    */

    unsigned int max_labels[4];

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_id, "/max-labels", H5P_DEFAULT);

        H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, max_labels);

        // Close the dataset
        H5Dclose(dataset_id);
    }

    //std::cout<<max_labels[0]<<" "<<max_labels[1]<<" "<<max_labels[2]<<" "<<max_labels[3]<<std::endl;

    // Open an existing group
    hid_t one_sets_group_id = H5Gopen(file_id, "/1-sets", H5P_DEFAULT);

    for (size_t i=1; i<max_labels[1]+1; i++)
    {
        std::stringstream datasetNameStream;
        datasetNameStream << i << "-1";
        std::stringstream binningGroupName;
        binningGroupName << "bin-" << i % numberOfBinningGroups;
        hid_t h = H5Gopen(one_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);

        // Open an existing dataset
        hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);

        //Datastructure to read in
        int junction_pos[3];

        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, junction_pos);

        point junction_p;
        junction_p.x=junction_pos[1]/2;
        junction_p.y=junction_pos[0]/2;

        junctions.push_back(junction_p);

        // Close the dataset
        H5Dclose(dataset_id);

        H5Gclose(h);
    }

    H5Gclose(one_sets_group_id);

    long unsigned int * arc_part_counters= new long unsigned int [max_labels[2]];

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_id, "/parts-counters-2", H5P_DEFAULT);

        H5Dread(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, arc_part_counters);

        // Close the dataset
        H5Dclose(dataset_id);
    }

    // Open an existing group
    hid_t two_sets_group_id = H5Gopen(file_id, "/2-sets", H5P_DEFAULT);

    arcs.resize(max_labels[2]);

    for (size_t i=1; i<max_labels[2]+1; i++)
    {
        std::vector<point> this_arc;

        std::stringstream binningGroupName;
        binningGroupName << "bin-" << i % numberOfBinningGroups;
        hid_t h = H5Gopen(two_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);

        for (size_t part=1; part<arc_part_counters[i-1]+1; part++)
        {
            std::stringstream datasetNameStream;
            datasetNameStream << i <<"-"<<part;//check for -2 and -3 ... and combine arc pixels
            
            // Open an existing dataset
            hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            // Datastructure to read in
            int * data_out_x = new int[dims[1]];
            int * data_out_y = new int[dims[1]];

            //dataspace for one coordinate
            hsize_t row[2];
            row[0] = 1;
            row[1] = dims[1];  
            hid_t mdataspace_id = H5Screate_simple(2, row, NULL);

            // select file hyperslab 
            hsize_t start_h[2];// start of hyperslab
            hsize_t count[2];// block count
            
            count[0]  = 1; 
            count[1]  = dims[1];

            start_h[0]  = 0; 
            start_h[1]  = 0;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_y);

            start_h[0]  = 1;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_x);

            for (int j=0;j<dims[1];j++)
            {
                point push;
                push.x=data_out_x[j];
                push.y=data_out_y[j];
                this_arc.push_back(push);
            }

            // Close the memoryspace
            H5Sclose(mdataspace_id);

            // Close the filespace
            H5Sclose(filespace);

            // Close the dataset
            H5Dclose(dataset_id);

            delete data_out_x;
            delete data_out_y;

        }

        for(int k=0;k<(int)this_arc.size();k++)
        {
            if(two_pixel_boundary==true)
            {
                point p1;
                point p2;
                //y is .5
                if(this_arc[k].x% 2 ==0 && this_arc[k].y% 2 !=0 )
                {
                    p1.x=(this_arc[k].x)/2;
                    p1.y=(this_arc[k].y)/2;

                    p2.x=(this_arc[k].x)/2;
                    p2.y=((this_arc[k].y)/2)+1;
                }
                //x is .5
                if(this_arc[k].x% 2 !=0 && this_arc[k].y% 2 ==0 )
                {
                    p1.x=(this_arc[k].x)/2;
                    p1.y=(this_arc[k].y)/2;

                    p2.x=(this_arc[k].x/2)+1;
                    p2.y=(this_arc[k].y)/2;
                }

                arcs[i-1].push_back(p1);
                arcs[i-1].push_back(p2);
            }
            if(two_pixel_boundary==false)
            {
                point p;

                p.x=(this_arc[k].x)/2;
                p.y=(this_arc[k].y)/2;

                arcs[i-1].push_back(p);

            }
        }

        H5Gclose(h);

    }

    H5Gclose(two_sets_group_id);

    // Close the file
    H5Fclose(file_id);

    std::cout<<"...done"<<std::endl;
}
