#include "itkImage.h"
#include "itkGreyLevelCooccurrenceMatrix.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGLCMImageCalculator.h"
#include "itkHaralickFeaturesGLCMCalculator.h"

#include <iostream>

const unsigned int W = 20;
const unsigned int H = 20;
const unsigned int D = 2;

typedef itk::Image<unsigned char, D> ImageType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIterator;

typedef itk::Statistics::GreyLevelCooccurrenceMatrix< unsigned int > GLCMType;

typedef itk::Statistics::GLCMImageCalculator< ImageType, GLCMType > GLCMImageCalculatorType;

void CreateImage(ImageType::Pointer img);

int main(int argc, char** argv)
{
  ImageType::Pointer image = ImageType::New();
  CreateImage(image);

  GLCMImageCalculatorType::Pointer calculator = GLCMImageCalculatorType::New();
  calculator->SetImage(image);
  calculator->SetMatrixSize(8);

  GLCMImageCalculatorType::OffsetVectorPointer offsets = GLCMImageCalculatorType::OffsetVectorType::New();
  GLCMImageCalculatorType::OffsetType offset1 = {{0, 1}};
  GLCMImageCalculatorType::OffsetType offset2 = {{1, 0}};
  offsets->push_back(offset1);
  offsets->push_back(offset2);
  calculator->SetOffsets(offsets);

  //calculator->DebugOn();

  calculator->SetRegion(image->GetLargestPossibleRegion());

  for(int i = 0; i < 1000; ++i)
  {
    calculator->ResetMatrix();
    calculator->Compute();
  }

  //calculator->GetCooccurrenceMatrix()->Print(std::cout);
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
    imageIt.Set(ind[1] % 2 + 1);
  }
}
