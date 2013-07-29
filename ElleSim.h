#pragma once

#include <iostream>     
#include <fstream>   
#include <sstream>

//A 2D point in an Elle data set
struct EllePoint2D
{
    int number;
    double x;
    double y;
    vigra::Point2D rasterized;
};

//A flynn in an Elle data set
struct ElleFlynn
{
    int number;
    std::vector<EllePoint2D> indices;
};

//An Elle data set
class ElleDataSet
{
    private:
    //Flynns
    std::vector<ElleFlynn> flynns;
    //Locations
    std::vector<EllePoint2D> locations;
    //Raw strings from the data set
    std::vector<std::string> inputLines;

    //Bounding box width/height
    double bBoxWidth;
    double bBoxHeight;

    //Maximum node separation
    double maxNodeSeparation;

    //Minimum node separation
    double minNodeSeparation;
    
    public:
    //Constructor with filepath argument
    ElleDataSet(std::string filepath);
    //Loads a data set
    void loadDataSet(std::string filepath);
    //Draws a .png image from the dataset with the given dimensions
    void exportToImage(int dim_x, int dim_y, std::string filepath);
};

ElleDataSet::ElleDataSet(std::string filepath)
{
    loadDataSet(filepath);
}

//----------------------------------------------------------------------
void ElleDataSet::loadDataSet(std::string filepath)
{
    std::cout << "Loading Elle simulation data set from file '" << filepath << "'" << std::endl;
    std::ifstream input(filepath.c_str());
    std::vector<std::string> flynnsLines;

    if(input)
    {
        std::cout << "Loading...";
        std::string inputLine;
        while(std::getline(input, inputLine))
        {
            inputLines.push_back(inputLine);
        }
        
        std::cout << "processing..." << std::endl;
        //The line numbers of the flynns and locations
        int flynnsLine = -1;
        int mineralLine = -1;
        int locationsLine = -1;
        int unodesLine = -1;
        int optionsLine = -1;
        for(int i = 0; i < inputLines.size(); i++)
        {
            if(inputLines[i].compare("OPTIONS"))
            {
                optionsLine = i;
                break;
            }
        }        
        for(int i = 0; i < inputLines.size(); i++)
        {
            if(inputLines[i].compare("FLYNNS") == 0)
            {
                flynnsLine = i;
                break;
            }
        }

        for(int i = flynnsLine; i < inputLines.size(); i++)
        {
            if(inputLines[i].compare("MINERAL") == 0)
            {
                mineralLine = i;
                break;
            }
        }
        for(int i = mineralLine; i < inputLines.size(); i++)
        {
            if(inputLines[i].compare("LOCATION") == 0)
            {
                locationsLine = i;
                break;
            }
        }
        for(int i = locationsLine; i < inputLines.size(); i++)
        {
            if(inputLines[i].compare("UNODES") == 0)
            {
                unodesLine = i;
                break;
            }
        }
        if(flynnsLine == -1 || locationsLine == -1)
        {
            std::cout << "Error reading the data set!" << std::endl;
            return;
        }

        std::cout << "getting options (maximum/minimum node separation, bounding box dimensions)..." << std::endl;
        //Bounding box edge points
        EllePoint2D bBox_LL;    //Lower left
        EllePoint2D bBox_LR;    //Lower right
        EllePoint2D bBox_UR;    //Upper right
        EllePoint2D bBox_UL;    //Upper left
        int bBoxLine = -1;
        for(int i = optionsLine + 1; i < flynnsLine; i++)
        {
            //Tokenize the line and store the substrings in a vector
            std::string line = inputLines[i];
            std::stringstream lineSS(line);
            
            std::istream_iterator<std::string> lineIt(lineSS);
            std::istream_iterator<std::string> lineEnd;
            std::vector<std::string> lineResults(lineIt, lineEnd);
            
            for(int j = 0; j < lineResults.size(); j++)
            {
                if(lineResults[j].compare("MaxNodeSeparation") == 0)
                {
                    std::istringstream(lineResults[j+1]) >> maxNodeSeparation;
                    std::cout << "Maximum node separation = " << maxNodeSeparation << std::endl;
                    break;
                }
                if(lineResults[j].compare("MinNodeSeparation") == 0)
                {
                    std::istringstream(lineResults[j+1]) >> minNodeSeparation;
                    std::cout << "Minimum node separation = " << minNodeSeparation << std::endl;
                    break;
                }
                if(lineResults[j].compare("CellBoundingBox") == 0)
                {
                    std::istringstream(lineResults[j+1]) >> bBox_LL.x;
                    std::istringstream(lineResults[j+2]) >> bBox_LL.y;
                    bBoxLine = i;
                    for(int k = i + 1; k <= i + 3; k++)
                    {
                        //Tokenize the line and store the substrings in a vector
                        std::string linek = inputLines[k];
                        std::stringstream lineSSk(linek);
            
                        std::istream_iterator<std::string> lineItk(lineSSk);
                        std::istream_iterator<std::string> lineEndk;
                        std::vector<std::string> lineResultsk(lineItk, lineEndk);
                        
                        if(k == i + 1)
                        {
                            std::istringstream(lineResultsk[0]) >> bBox_LR.x;
                            std::istringstream(lineResultsk[1]) >> bBox_LR.y;
                        }
                        if(k == i + 2)
                        {
                            std::istringstream(lineResultsk[0]) >> bBox_UR.x;
                            std::istringstream(lineResultsk[1]) >> bBox_UR.y;
                        }
                        if(k == i + 3)
                        {
                            std::istringstream(lineResultsk[0]) >> bBox_UL.x;
                            std::istringstream(lineResultsk[1]) >> bBox_UL.y;
                        }
                    }
                    bBoxWidth = bBox_LR.x - bBox_LL.x;
                    bBoxHeight = bBox_UL.y - bBox_LL.y;
                    std::cout << "Bounding box dimensions: " << bBoxWidth << " (horizontal), " << bBoxHeight << " (vertical)" << std::endl;
                    i = i + 3;
                    break;
                }
            }
        }
        
        std::cout << "loading data into structures..." << std::endl;
        //Get the locations line per line
        for(int i = locationsLine + 1; i < unodesLine; i++)
        {
            EllePoint2D p_temp;
            std::string line = inputLines[i];
            std::stringstream lineSS(line);
    
            std::istream_iterator<std::string> lineIt(lineSS);
            std::istream_iterator<std::string> lineEnd;
            std::vector<std::string> lineResults(lineIt, lineEnd);
            
            std::istringstream(lineResults[0]) >> p_temp.number;
            std::istringstream(lineResults[1]) >> p_temp.x;
            std::istringstream(lineResults[2]) >> p_temp.y;
            
            locations.push_back(p_temp);
        }
        
        //Get the flynns line per line
        for(int i = flynnsLine + 1; i < mineralLine; i++)
        {
            ElleFlynn flynn;
            //Tokenize the line and store the substrings in a vector
            std::string line = inputLines[i];
            std::stringstream lineSS(line);
            
            std::istream_iterator<std::string> lineIt(lineSS);
            std::istream_iterator<std::string> lineEnd;
            std::vector<std::string> lineResults(lineIt, lineEnd);
            
            std::vector<EllePoint2D> pV;
            for(int j = 0; j < lineResults.size(); j++)
            {
                if(j == 0)
                {
                    std::istringstream(lineResults[j]) >> flynn.number;
                }
                else
                {
                    EllePoint2D p;
                    std::istringstream(lineResults[j]) >> p.number;
                    pV.push_back(p);
                }
            }
            //Get the locations belonging to the flynn nodes
            for(int j = 0; j < pV.size(); j++)
            {
                for(int k = 0; k < locations.size(); k++)
                {
                    if(pV[j].number == locations[k].number)
                    {
                        pV[j].x = locations[k].x;
                        pV[j].y = locations[k].y;
                    }
                }
            }
            flynn.indices = pV;
            flynns.push_back(flynn);
        }
        
        std::cout << "done! Loaded " << flynns.size() << " flynns and " << locations.size() << " locations." << std::endl;
    }
    else
    {
        std::cout << "Error loading the file!" << std::endl;
        return;
    }
}

//----------------------------------------------------------------------
void ElleDataSet::exportToImage(int dim_x, int dim_y, std::string filepath)
{
    //Calculate the scaling for the x- and y-axis
    double scalingX;
    double scalingY;
    //Image dimensions, either explicitly specified or implicity derived from the scaling
    int img_x;
    int img_y;

    //Only width in pixels specified
    if(dim_x > 0 && dim_y == -1)
    {
        scalingX = (double)dim_x/bBoxWidth;
        scalingY = bBoxHeight/bBoxWidth;
        img_x = dim_x;
        img_y = dim_x*scalingY;
        //Correct the scaling for the y-axis
        scalingY = (double)img_y/bBoxHeight;
    }
    //Only height in pixels specified
    if(dim_x == -1 && dim_y > 0)
    {
        scalingX = bBoxWidth/bBoxHeight;
        scalingY = (double)dim_y/bBoxHeight;
        img_x = dim_y*scalingX;
        img_y = dim_y;
        //Correct the scaling for the x-axis
        scalingX = (double)img_x/bBoxWidth;
    }
    //Both width and height in pixels specified
    if(dim_x > 0 && dim_y > 0)
    {
        scalingX = (double)dim_x/bBoxWidth;
        scalingY = (double)dim_y/bBoxHeight;
        img_x = dim_x;
        img_y = dim_y;
    }
    else if(dim_x == 0 || dim_y == 0 || dim_x < -1 || dim_y < -1)
    {
        std::cout << "Error: invalid dimensions (" << dim_x << ", " << dim_y << ")" << std::endl;
        return; 
    }
    //Calculate the optimal dimensions, so that no information from the dataset gets lost during the rasterization (start with minNodeSeparation = 2px on the x-axis)
    else if(dim_x == -1 && dim_y == -1)
    {
        dim_x = (double)2*bBoxWidth/minNodeSeparation;
        scalingX = (double)dim_x/bBoxWidth;
        scalingY = bBoxHeight/bBoxWidth;
        img_x = dim_x;
        img_y = dim_x*scalingY;
        //Correct the scaling for the y-axis
        scalingY = (double)img_y/bBoxHeight;
    }

    /*std::cout << "Resulting image dimensions: " << img_x << "x" << img_y << ", scaling: " << scalingX << " (horizontal), " << scalingY << " (vertical)" << std::endl;*/

    //Rasterize the node locations
    for(int i = 0; i < flynns.size(); i++)
    {
        //The flynn's number will be the label
        int label = flynns[i].number;

        //Iterate through all flynn nodes
        for(int j = 0; j < flynns[i].indices.size(); j++)
        {
            vigra::Point2D jPoint;
            jPoint.x = flynns[i].indices[j].x*scalingX;
            jPoint.y = (bBoxHeight - flynns[i].indices[j].y)*scalingY;

            flynns[i].indices[j].rasterized = jPoint;
        }
    }

    //Search for unconnected flynns
    int max_diff_x = 2*maxNodeSeparation*scalingX;
    int max_diff_y = 2*maxNodeSeparation*scalingY;
    std::vector<ElleFlynn> originalFlynns = flynns;
    for(int i = 0; i < flynns.size(); i++)
    {
        std::vector<EllePoint2D> originalPositions = originalFlynns[i].indices;
        //Iterate through all flynn nodes, get the rasterized positions
        for(int j = 1; j < originalPositions.size(); j++)
        {
            //Check if there are two points that are on opposite ends of the image (horizontal
            if(abs(originalPositions[j].rasterized.x - originalPositions[j-1].rasterized.x) >= img_x - max_diff_x)
            {
                //Naive approach; certainly, this can be improved, but at least it yields correct results
                for(int k = 0; k < originalPositions.size(); k++)
                {
                    if(originalPositions[k].rasterized.x + img_x >= 2*img_x - originalPositions[k].rasterized.x - max_diff_x)
                    {
                        //Do nothing for rasterized node positions that are already on the right-hand side of the image
                    }
                    else
                    {
                        flynns[i].indices[k].rasterized.x = originalPositions[k].rasterized.x + img_x;
                    }
                }
            }
            //Do the same for the vertical opposite ends
            if(abs(originalPositions[j].rasterized.y - originalPositions[j-1].rasterized.y) >= img_y - max_diff_y)
            {
                //Naive approach; certainly, this can be improved, but at least it yields correct results
                for(int k = 0; k < originalPositions.size(); k++)
                {
                    if(originalPositions[k].rasterized.y + img_y >= 2*img_y - originalPositions[k].rasterized.y - max_diff_y)
                    {
                        //Do nothing for rasterized node positions that are already on the right-hand side of the image
                    }
                    else
                    {
                        flynns[i].indices[k].rasterized.y = originalPositions[k].rasterized.y + img_y;
                    }
                }
            }
        }
    }
    
    //Remove outliers/artifacts by comparing the horizontal and vertical distance to the previous node
    int maxDistance = 3*maxNodeSeparation*scalingX*3*maxNodeSeparation*scalingX + 3*maxNodeSeparation*scalingY*3*maxNodeSeparation*scalingY;
    int artifactsCount = 0;
    std::cout << "Removing outliers/artifacts, maximum node separation: " << ceil(sqrt(maxDistance)) << "px"<< std::endl;  
    for(int i = 0; i < flynns.size(); i++)
    {
        std::vector<EllePoint2D> artifacts;
        bool artifactsFound = false;
        for(int j = 1; j < flynns[i].indices.size(); j++)
        {
            int diff_x = abs(flynns[i].indices[j].rasterized.x - flynns[i].indices[j-1].rasterized.x);
            int diff_y = abs(flynns[i].indices[j].rasterized.y - flynns[i].indices[j-1].rasterized.y);
            if(diff_x*diff_x + diff_y*diff_y > maxDistance && (diff_x < img_x - 1 || diff_y < img_y - 1))
            {
                artifactsFound = true;
                artifacts.push_back(flynns[i].indices[j-1]);
                artifactsCount++;
            }
        }
        //Find all artifacts in the flynn
        while(artifactsFound == true)
        {
            artifactsFound = false;
            for(int j = 1; j < flynns[i].indices.size(); j++)
            {
                int diff_x = abs(flynns[i].indices[j].rasterized.x - flynns[i].indices[j-1].rasterized.x);
                int diff_y = abs(flynns[i].indices[j].rasterized.y - flynns[i].indices[j-1].rasterized.y);
                if(diff_x*diff_x + diff_y*diff_y > maxDistance && (diff_x < img_x - 1 || diff_y < img_y - 1))
                {
                    bool doNotAdd = false;
                    for(int k = 0; k < artifacts.size(); k++)
                    {
                        if(artifacts[k].number == flynns[i].indices[j-1].number)
                        {
                            doNotAdd = true;
                        }
                    }
                    if(doNotAdd == false)
                    {
                          artifacts.push_back(flynns[i].indices[j-1]);
                          artifactsFound = true;
                          artifactsCount++;
                    }
                }
            }
        }
        //Delete the artifacts if at least one artifact has been found
        if(artifacts.size() > 0)
        {
            for(int k = 0; k < flynns[i].indices.size(); k++)
            {
                for(int l = 0; l < artifacts.size(); l++)
                {
                    if(flynns[i].indices[k].number == artifacts[l].number)
                    {
                        flynns[i].indices.erase(flynns[i].indices.begin() + k);
                        k--;
                    }
                }
            }
        }
    }
    std::cout << "Removed " << artifactsCount << " outliers/artifacts" << std::endl;

    //Temporary label image to find the width and height the image needs to be extended to in order to fit all grains in it
    vigra::IImage label_image_t(img_x*2, img_y*2);

    //Draw the flynns into the temporary label image
    for(int i = 0; i < flynns.size(); i++)
    {
        //Iterate through all flynn nodes
        for(int j = 0; j < flynns[i].indices.size(); j++)
        {
            vigra::Point2D jPoint = flynns[i].indices[j].rasterized;
            label_image_t(jPoint.x, jPoint.y) = 255;
        }
    }

    //Get the maximum width of the image
    int maxWidth = 0;
    for(int y = 0; y < 2*img_y; y++)
    {
        for(int x = img_x; x < 2*img_x; x++)
        {
            if(label_image_t(x, y) == 255)
            {
                if(maxWidth < x)
                {
                    maxWidth = x;
                }
            }
        }
    }
   
    //Get the maximum height of the image
    int maxHeight = 0;
    for(int y = img_y; y < 2*img_y; y++)
    {
        for(int x = 0; x < 2*img_x; x++)
        {
            if(label_image_t(x, y) == 255)
            {
                if(maxHeight < y)
                {
                    maxHeight = y;
                }
            }
        }
    }

    //The final label image with a 2 pixel safety margin
    vigra::IImage label_image(maxWidth + 2, maxHeight + 2);
    
    //Draw the flynns into the label image
    for(int i = 0; i < flynns.size(); i++)
    {
        //Iterate through all flynn nodes and draw splines between them
        for(int j = 1; j <= flynns[i].indices.size(); j++)
        {
            vigra::Point2D jPoint = flynns[i].indices[j%flynns[i].indices.size()].rasterized;
            vigra::Point2D jPointPrev = flynns[i].indices[(j-1)%flynns[i].indices.size()].rasterized;
            
            if(j == flynns[i].indices.size())
            {
                int maxDiffX = (maxNodeSeparation*scalingX) + 0.2*maxNodeSeparation*scalingX;
                int maxDiffY = (maxNodeSeparation*scalingY) + 0.2*maxNodeSeparation*scalingY;
                
                double dist = sqrt(abs(jPoint.x - jPointPrev.x)*abs(jPoint.x - jPointPrev.x) + abs(jPoint.y - jPointPrev.y)*abs(jPoint.y - jPointPrev.y));
                double distDiff = sqrt(maxDiffX*maxDiffX + maxDiffY*maxDiffY);
                if(ceil(dist) > ceil(distDiff))
                {
                    break;
                }
            }

            
            label_image(jPoint.x, jPoint.y) = 255;
            label_image(jPointPrev.x, jPointPrev.y) = 255;
            
            int diff_x;
            int diff_y;
                        
            if(jPoint.x > jPointPrev.x && abs(jPoint.y - jPointPrev.y) <= abs(jPoint.x - jPointPrev.x))
            {
                diff_x = jPoint.x - jPointPrev.x;
                diff_y = jPoint.y - jPointPrev.y;
                for(int x = jPointPrev.x; x <= jPoint.x; x++)
                {
                    int y = ceil(((double)diff_y/(double)diff_x)*x + jPointPrev.y - ((double)diff_y/(double)diff_x)*jPointPrev.x);
                    label_image(x, y) = 255;
                }
            }
            else if(jPoint.x < jPointPrev.x && abs(jPoint.y - jPointPrev.y) <= abs(jPoint.x - jPointPrev.x))
            {
                diff_x = jPointPrev.x - jPoint.x;
                diff_y = jPointPrev.y - jPoint.y;
                for(int x = jPoint.x; x <= jPointPrev.x; x++)
                {
                    int y = ceil(((double)diff_y/(double)diff_x)*x + jPoint.y - ((double)diff_y/(double)diff_x)*jPoint.x);
                    label_image(x, y) = 255;
                }
            }
            else if(jPoint.x == jPointPrev.x)
            {
                if(jPoint.y > jPointPrev.y)
                {
                    for(int y = jPointPrev.y; y <= jPoint.y; y++)
                    {
                        label_image(jPoint.x, y) = 255;
                    }
                }
                if(jPoint.y < jPointPrev.y)
                {
                    for(int y = jPoint.y; y <= jPointPrev.y; y++)
                    {
                        label_image(jPoint.x, y) = 255;
                    }
                }                
            }
            else if(abs(jPoint.y - jPointPrev.y) > abs(jPoint.x - jPointPrev.x) && jPoint.x != jPointPrev.x)
            {
                if(jPoint.y > jPointPrev.y)
                {
                    diff_x = jPoint.x - jPointPrev.x;
                    diff_y = jPoint.y - jPointPrev.y;
                    for(int y = jPointPrev.y; y <= jPoint.y; y++)
                    {
                        int x = ceil(((double)diff_x/(double)diff_y)*y + jPointPrev.x - ((double)diff_x/(double)diff_y)*jPointPrev.y);
                        label_image(x, y) = 255;
                    }
                }
                if(jPoint.y < jPointPrev.y)
                {
                    diff_x = jPoint.x - jPointPrev.x;
                    diff_y = jPoint.y - jPointPrev.y;
                    for(int y = jPoint.y; y <= jPointPrev.y; y++)
                    {
                        int x = ceil(((double)diff_x/(double)diff_y)*y + jPoint.x - ((double)diff_x/(double)diff_y)*jPoint.y);
                        label_image(x, y) = 255;
                    }
                }
            }
        }
    }
    
    //Export the label image
    std::cout << "Exporting image to '" << filepath << "'" << std::endl;
    exportImage(srcImageRange(label_image), vigra::ImageExportInfo(filepath.c_str()));
}