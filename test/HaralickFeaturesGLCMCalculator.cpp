#include "itkImage.h"
#include "itkGreyLevelCooccurrenceMatrix.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGLCMImageCalculator.h"
#include "itkHaralickFeaturesGLCMCalculator.h"
#include "itkImageFileWriter.h"

#include <iostream>

const unsigned int W = 20;
const unsigned int H = 20;
const unsigned int D = 2;

typedef itk::Image<unsigned char, D> ImageType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIterator;

typedef itk::Statistics::GreyLevelCooccurrenceMatrix< unsigned int > GLCMType;

typedef itk::Statistics::GLCMImageCalculator< ImageType, GLCMType > GLCMImageCalculatorType;

typedef itk::Statistics::HaralickFeaturesGLCMCalculator< GLCMType, double > FeaturesCalculatorType;

typedef itk::ImageFileWriter<ImageType> WriterType;

void CreateImage(ImageType::Pointer img);

int main(int argc, char** argv)
{
  ImageType::Pointer image = ImageType::New();
  CreateImage(image);

  GLCMImageCalculatorType::Pointer calculator = GLCMImageCalculatorType::New();
  calculator->SetImage(image);
  calculator->SetMatrixSize(4);

  GLCMImageCalculatorType::OffsetVectorPointer offsets = GLCMImageCalculatorType::OffsetVectorType::New();
  //GLCMImageCalculatorType::OffsetType offset1 = {{0, 1}};
  GLCMImageCalculatorType::OffsetType offset2 = {{0, 1}};
  //offsets->push_back(offset1);
  offsets->push_back(offset2);
  calculator->SetOffsets(offsets);

  calculator->SetRegion(image->GetLargestPossibleRegion());

  calculator->Compute();

  calculator->GetCooccurrenceMatrix()->Print(std::cout);

  FeaturesCalculatorType::Pointer featuresCalculator = FeaturesCalculatorType::New();
  featuresCalculator->SetCooccurrenceMatrix(calculator->GetCooccurrenceMatrix());

  featuresCalculator->Compute();

  std::cout << featuresCalculator->GetFeatures() << std::endl;
}

void CreateImage(ImageType::Pointer image)
{
  ImageType::RegionType region;

  {
    ImageType::IndexType index = {{ 0, 0 }};
    ImageType::SizeType imageSize = {{ W, H }};

    region.SetIndex( index );
    region.SetSize( imageSize );
  }

  image->SetRegions( region );
  image->Allocate();

  ImageIterator imageIt( image, region );
  ImageType::IndexType ind;
  for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt)
  {
    ind = imageIt.GetIndex();
	//if(ind[1] < (W/2))
	//    imageIt.Set(ind[1] % 2);
	//else
	//    imageIt.Set(ind[0] % 2);
	imageIt.Set(ind[1] % 2);
  }

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(image);
  writer->SetFileName("out.bmp");
  writer->Update();
}
