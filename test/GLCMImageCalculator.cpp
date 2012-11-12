#include "itkGLCMImageCalculator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGreyLevelCooccurrenceMatrix.h"

#include <iostream>

const unsigned int W = 200;
const unsigned int H = 200;
const unsigned int D = 2;

typedef itk::Image<unsigned char, D> Image2DType;
typedef itk::ImageRegionIteratorWithIndex< Image2DType > Image2DIterator;

typedef itk::Statistics::GreyLevelCooccurrenceMatrix< unsigned int > GLCMType;

typedef itk::GLCMImageCalculator< Image2DType, GLCMType > GLCMImageCalculatorType;

void CreateImage(Image2DType::Pointer img);

int main(int argc, char** argv)
{
  Image2DType::Pointer image = Image2DType::New();
  CreateImage(image);

  GLCMImageCalculatorType::Pointer calculator = GLCMImageCalculatorType::New();
  calculator->SetImage(image);
  calculator->SetMatrixSize(8);

  Image2DType::OffsetType offset1 = {{0, 1}};
  calculator->SetOffset(offset1);

  calculator->SetRegion(image->GetLargestPossibleRegion());

  calculator->Compute();
}

void CreateImage(Image2DType::Pointer image)
{
  Image2DType::RegionType region;

  {
    Image2DType::IndexType index = {{ 0, 0 }};
    Image2DType::SizeType imageSize = {{ W, H }};

    region.SetIndex( index );
    region.SetSize( imageSize );
  }

  image->SetRegions( region );
  image->Allocate();

  Image2DIterator imageIt( image, region );
  Image2DType::IndexType ind;
  for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt)
  {
    ind = imageIt.GetIndex();
    imageIt.Set(ind[1] % 2 + 1);
  }
}
