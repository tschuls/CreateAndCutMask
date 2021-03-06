/*=========================================================================
*
*  Copyright Insight Software Consortium
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImageFileWriter.h"

#include <vector>
#include "itksys/SystemTools.hxx"

#include "itkPointSet.h"

#include "itkImageRegionIterator.h"

#include "itkVector.h"
#include "itkMatrix.h"

using namespace std;

typedef signed short    PixelType;
const unsigned int      Dimension = 3;
typedef itk::Image< PixelType, Dimension >      ImageType;
typedef itk::ImageSeriesReader< ImageType >     ReaderType;
typedef ImageType::Pointer ImagePointerType;

typedef itk::PointSet< double, 3 >   PointSetType;
typedef PointSetType::PointType DoublePointType;
typedef PointSetType::PointsContainer PointsContainerType;
typedef PointsContainerType::Pointer PointsContainerPointer;

typedef itk::Vector<double, 3> VectorType;
typedef itk::Point<double,3> PointType;


PointsContainerPointer readLandmarksWithOrigin(string filename, ImagePointerType image) {
    if (filename==""){
        return NULL;
    }

    std::vector<double> extent(Dimension);
    ImageType::SpacingType spacing = image->GetSpacing();
    ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    for (int d=0;d<Dimension;++d) {
        extent[d] = spacing[d] * size[d];
    }


    PointSetType::Pointer  pointSet = PointSetType::New();
    PointsContainerPointer points=pointSet->GetPoints();

    ifstream ifs(filename.c_str());

    //read origin
    DoublePointType origin;
    for (int d=0;d<Dimension;++d) {
        ifs>>origin[d];
    }

    int i=0;
    while ( ! ifs.eof() ) {
        DoublePointType point;
        bool fullPoint=true;
        for (int d=0;d<Dimension;++d){
            ifs>>point[d];
            if (ifs.eof()){
                fullPoint=false;
                break;
            }

        }
        if (fullPoint){
          //do not move the point by the origin, since the DICOM
          //images read by this tool are still correctly aligned
          //and the stored landmarks are already correct!
          points->InsertElement(i, point);
          ++i;
        }
    } 
    return points;
}

bool isInside(PointType pointECS, double width, double length, double depth, double marginWidth, double marginLength, double marginDepth) {
    //check if x smaller than width+marginWidth
    //check if y smaller than length+marginLength
    //check if z smaller than depth+marginDepth
    if ((0-marginWidth) < pointECS[0] && pointECS[0] < (width+marginWidth)) {
        if ((0-marginLength) < pointECS[1] && pointECS[1] < (length+marginLength)) { 
            if ( -(depth+ marginDepth) < pointECS[2] && pointECS[2] < (0+ marginDepth)) { 
              return true;
            }
        }
    }
    return false;

}

void printPoint(const PointType& point){
    cout << "point: " << point[0] << " " << point[1] << " " << point[2] << endl;
}


int main( int argc, char* argv[] )
{
    if( argc < 4 )
    {
        std::cerr << "Usage: " << argv[0] <<
            "InputDicomDirectory  OutputDicomDirectory landmarkFile" << std::endl;
        std::cerr << "Format of the landmark file: " << std::endl;
        std::cerr << "    (origin)    e.g. -239.5 -239.5 -183.75 " << std::endl;
        std::cerr << "    (nose)      e.g. -4.26322 -127.767 70.1066  " << std::endl;
        std::cerr << "    (eye1)      e.g. -55.0607 -74.4404 92.9136  " << std::endl;
        std::cerr << "    (eye2)      e.g.  48.6076 -73.5048 92.9136 " << std::endl;
        std::cerr << "    (ganthion)  e.g. -7.37327 -108.12 -1.42449 " << std::endl;
        return EXIT_FAILURE;
    }
    
    
    //****************************************************************
    //****** read input image ****************************************
    //****************************************************************
    typedef itk::GDCMImageIO                        ImageIOType;
    typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

    ImageIOType::Pointer gdcmIO = ImageIOType::New();
    NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();

    namesGenerator->SetInputDirectory( argv[1] );
    const ReaderType::FileNamesContainer & filenames = namesGenerator->GetInputFileNames();

    const unsigned int numberOfFileNames =  filenames.size();
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetImageIO( gdcmIO );
    reader->SetFileNames( filenames );

    try
    {
        reader->Update();
    }
    catch (itk::ExceptionObject &excp)
    {
        std::cerr << "Exception thrown while writing the image" << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }

    ImageType::Pointer inputImage = reader->GetOutput();



    //****************************************************************
    //****** read landmarks file *************************************
    //****************************************************************
    PointsContainerPointer landmarks = readLandmarksWithOrigin(argv[3], inputImage);

    PointType eye1 = landmarks->GetElement(1);
    PointType eye2 = landmarks->GetElement(2);
    PointType gnathion = landmarks->GetElement(3);
    PointType nose = landmarks->GetElement(0);

    printPoint(eye1);
    printPoint(eye2);
    printPoint(gnathion);
    printPoint(nose);

    //****************************************************************
    //****** calculate values for inside/outside computation *********
    //****************************************************************    
    VectorType v1;
    v1 = eye2 - eye1;

    VectorType v2;
    v2 = gnathion - eye1;

    VectorType v3 = itk::CrossProduct(v1,v2);

    v2 = itk::CrossProduct(v3,v1);

    
    
    v1.Normalize();
    v2.Normalize();
    v3.Normalize();

    //create matrix
    typedef itk::Matrix<double, 3, 3> MatrixType;
    MatrixType M;
    M(0,0) = v1[0];
    M(1,0) = v1[1];
    M(2,0) = v1[2];
    M(0,1) = v2[0];
    M(1,1) = v2[1];
    M(2,1) = v2[2];
    M(0,2) = v3[0];
    M(1,2) = v3[1];
    M(2,2) = v3[2];
  
    MatrixType inverseM = M.GetInverse();

    PointType eye1ECS = inverseM * eye1;
    PointType eye2ECS = inverseM * eye2;
    PointType gnathionECS = inverseM * gnathion;
    PointType noseECS = inverseM * nose;
    
    cout << "ECS" << endl;
    printPoint(eye1ECS);
    printPoint(eye2ECS);
    printPoint(gnathionECS);
    printPoint(noseECS);

    double width = abs(eye2ECS[0] - eye1ECS[0]);
    double length = abs(gnathionECS[1] - eye1ECS[1]);
    double depth = abs(noseECS[2]- eye1ECS[2]);
    
    cout << "box has sizing: " << endl;
    cout << "  width: " << width << endl;
    cout << "  length: " << length << endl;
    cout << "  depth: " << depth << endl;


    printPoint(inputImage->GetOrigin());

    ImagePointerType maskImage = ImageType::New();
    maskImage->SetRegions(inputImage->GetLargestPossibleRegion());
    maskImage->SetOrigin(inputImage->GetOrigin());
    maskImage->SetSpacing(inputImage->GetSpacing());
    maskImage->SetDirection(inputImage->GetDirection());
    maskImage->Allocate();

    //****************************************************************
    //****** iterate over image and mask *****************************
    //****************************************************************      
    typedef itk::ImageRegionIterator<ImageType> ItType;
    ItType it( inputImage, inputImage->GetLargestPossibleRegion() );
    it.GoToBegin();

    ItType maskIt(maskImage, maskImage->GetLargestPossibleRegion());
    maskIt.GoToBegin();

    while (!it.IsAtEnd()) {
        PointType queryPoint;
        inputImage->TransformIndexToPhysicalPoint(it.GetIndex(),queryPoint);
        
        //printPoint(queryPoint);
        
        VectorType newPoint = queryPoint - eye1;
        VectorType resultPoint = inverseM * newPoint;

        //printPoint(newPoint);
        //printPoint(resultPoint);

        if (isInside(resultPoint, width, length, depth, 15.0,5.0,20.0)) {
          inputImage->SetPixel(it.GetIndex(), -1000);
          maskImage->SetPixel(maskIt.GetIndex(), 0);
        } else {
          maskImage->SetPixel(maskIt.GetIndex(), 1);
        }
        ++it;
        ++maskIt;
    }

    //****************************************************************
    //****** write mask out  *****************************************
    //**************************************************************** 
    typedef  itk::ImageFileWriter< ImageType  > WriterType;
    WriterType::Pointer maskWriter = WriterType::New();
    string name = "\\mask.vtk";
    string maskFilename = argv[2] + name;
    maskWriter->SetFileName(maskFilename);
    maskWriter->SetInput(maskImage);
    maskWriter->Update();

    //****************************************************************
    //****** write masked image out **********************************
    //**************************************************************** 
    const char * outputDirectory = argv[2];
    itksys::SystemTools::MakeDirectory( outputDirectory );
    typedef signed short    OutputPixelType;
    const unsigned int      OutputDimension = 2;
    typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
    typedef itk::ImageSeriesWriter<ImageType, Image2DType >  SeriesWriterType;
    SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

    //write identical headers also with original UID's -> besides the changed pixel values the files look the same!
    seriesWriter->SetInput( inputImage );

    gdcmIO->SetKeepOriginalUID(true);
    seriesWriter->SetImageIO( gdcmIO );
    namesGenerator->SetOutputDirectory( outputDirectory );
    seriesWriter->SetFileNames( namesGenerator->GetOutputFileNames() );
    seriesWriter->SetMetaDataDictionaryArray(reader->GetMetaDataDictionaryArray() );
    try
    {
        seriesWriter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << "Exception thrown while writing the series " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;

}

