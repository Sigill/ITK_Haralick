#include "itkGreyLevelCooccurrenceMatrix.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGLCMImageCalculator.h"
#include "itkHaralickFeaturesGLCMCalculator.h"

#include <iostream>

const unsigned int W = 20;
const unsigned int H = 20;
const unsigned int D = 2;

typedef itk::Image<unsigned char, D> Image2DType;
typedef itk::ImageRegionIteratorWithIndex< Image2DType > Image2DIterator;

typedef itk::Statistics::GreyLevelCooccurrenceMatrix< unsigned int > GLCMType;

typedef itk::GLCMImageCalculator< Image2DType, GLCMType > GLCMImageCalculatorType;

typedef itk::HaralickFeaturesGLCMCalculator < GLCMType, double > FeaturesCalculatorType;

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

  calculator->DebugOn();

  calculator->SetRegion(image->GetLargestPossibleRegion());

  calculator->Compute();

  calculator->GetCooccurrenceMatrix()->Print(std::cout);

  FeaturesCalculatorType::Pointer featuresCalculator = FeaturesCalculatorType::New();
  featuresCalculator->SetCooccurrenceMatrix(calculator->GetCooccurrenceMatrix());
  featuresCalculator->Compute();

  std::cout << featuresCalculator->GetFeatures() << std::endl;
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
