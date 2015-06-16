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
            //move point by: (Landmark-Origin)-Extent/2
            point = point - origin;
            for (int d=0;d<Dimension;++d) {
                point[d] = point[d] - extent[d]/2;
            }
            points->InsertElement(i, point);
            ++i;
        }
    } 
    return points;
}

bool isInside(PointType point, double width, double height, double depth) {
    //check if x smaller than width
    //check if y smaller than height
    //check if z smaller than depth
    if (0 < point[0] && point[0] < width) { //x is width
        if (0 < point[1] && point[1] < height) { //z is height
            if (0 < point[2] && point[2] < depth) { //y is depth
                return true;
            }
        }
    }
    return false;

}

int main( int argc, char* argv[] )
{
    if( argc < 4 )
    {
        std::cerr << "Usage: " << argv[0] <<
            " DicomInputDirectory  OutputDicomDirectory landmarkFile" << std::endl;
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
    //std::cout << numberOfFileNames << std::endl;
    //for(unsigned int fni = 0; fni < numberOfFileNames; ++fni)
    //{
    //    std::cout << "filename # " << fni << " = ";
    //    std::cout << filenames[fni] << std::endl;
    //}
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
    PointType mandible = landmarks->GetElement(3);
    PointType nose = landmarks->GetElement(0);

    //****************************************************************
    //****** calculate values for inside/outside computation *********
    //****************************************************************    
    VectorType v1;
    v1 = eye2 - eye1;

    VectorType v2;
    v2 = mandible - eye1;

    VectorType v3 = itk::CrossProduct(v2,v1);

    double width = abs(eye2[0] - eye1[0]);
    double heigth = abs(mandible[2] - eye1[2]);
    double depth = abs(nose[1]- eye1[1]);
    
    cout << "box has sizing: " << endl;
    cout << "  width: " << width << endl;
    cout << "  height: " << heigth << endl;
    cout << "  depth: " << depth << endl;

    
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

    //****************************************************************
    //****** iterate over image **************************************
    //****************************************************************      
    typedef itk::ImageRegionIterator<ImageType> ItType;
    ItType it( inputImage, inputImage->GetLargestPossibleRegion() );
    it.GoToBegin();

    while (!it.IsAtEnd()) {
        PointType queryPoint;
        inputImage->TransformIndexToPhysicalPoint(it.GetIndex(),queryPoint);

        PointType newPoint = queryPoint - eye1;
        PointType resultPoint = inverseM * newPoint;

        if (isInside(resultPoint, width,heigth,depth)) {
            inputImage->SetPixel(it.GetIndex(), -1000);
        }
        ++it;
    }


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

