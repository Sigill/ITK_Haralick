#include <iostream>

#include <itkImageRegionIteratorWithIndex.h>

#include "itkScalarImageToCooccurrenceMatrixFilter.h"

const unsigned int W = 200;
const unsigned int H = 200;
const unsigned int D = 2;

typedef itk::Image<unsigned char, D> Image2DType;
typedef itk::ImageRegionIteratorWithIndex< Image2DType > Image2DIterator;
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter< Image2DType > Image2DCoocurrenceMatrixComputer;

int main(int argc, char **argv)
{
  Image2DType::Pointer image = Image2DType::New();

  Image2DType::RegionType region;

  {
    Image2DType::IndexType index = {{ 0, 0 }};
    Image2DType::SizeType imageSize = {{ W, H }};

    region.SetIndex( index );
    region.SetSize( imageSize );
  }

  image->SetRegions( region );
  image->Allocate(); 

  //--------------------------------------------------------------------------
  //  1 2 1 2 1
  //  1 2 1 2 1
  //  1 2 1 2 1
  //  1 2 1 2 1
  //  1 2 1 2 1
  //--------------------------------------------------------------------------
  {
    Image2DIterator imageIt( image, region );
    Image2DType::IndexType ind;
    for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt)
    {
      ind = imageIt.GetIndex();
      imageIt.Set(ind[1] % 2 + 1);
    }
  }


  Image2DCoocurrenceMatrixComputer::Pointer coocurrenceMatrixComputer = Image2DCoocurrenceMatrixComputer::New();

  Image2DType::RegionType roi;
  Image2DType::IndexType rIndex = {{ 0, 0 }};
  Image2DType::SizeType rSize = {{ 32, 32 }};

  roi.SetSize( rSize );
  roi.SetIndex( rIndex );

  coocurrenceMatrixComputer->SetNumberOfBinsPerAxis(4); //reasonable number of bins

  coocurrenceMatrixComputer->SetRegionOfInterest(roi);

  coocurrenceMatrixComputer->SetInput(image);

  Image2DType::OffsetType offset1 = {{0, 1}};
  Image2DCoocurrenceMatrixComputer::OffsetVectorPointer offsetV = Image2DCoocurrenceMatrixComputer::OffsetVector::New();
  offsetV->push_back(offset1);
  coocurrenceMatrixComputer->SetOffsets(offsetV);

  coocurrenceMatrixComputer->Update();

  for(int i = 0; i < 1000; ++i) 
  {
    image->Modified();
    coocurrenceMatrixComputer->Update();
    //std::cout << "Requested Region:" << std::endl << image->GetRequestedRegion() << std::endl;
    //std::cout << "Requested Region:" << std::endl << image->GetRequestedRegion() << std::endl;
  }
  const Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType* cooc = coocurrenceMatrixComputer->GetOutput();

  typename Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType::ConstIterator coocIt = cooc->Begin(), coocBegin = cooc->Begin(), coocEnd = cooc->End();

  while(coocIt != coocEnd)
  {
    //std::cout << coocIt.GetIndex() << " -> " << coocIt.GetFrequency() << std::endl;
    std::cout << ( (unsigned int)(coocIt - coocBegin) ) << " -> " << (*coocIt) << std::endl;
    ++coocIt;
  }

  /*
  Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType::IndexType one_one( hist->GetMeasurementVectorSize() );
  Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType::IndexType one_two( hist->GetMeasurementVectorSize() );
  Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType::IndexType two_one( hist->GetMeasurementVectorSize() );
  Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType::IndexType two_two( hist->GetMeasurementVectorSize() );

  one_one[0] = 1;
  one_one[1] = 1;

  one_two[0] = 1;
  one_two[1] = 2;

  two_one[0] = 2;
  two_one[1] = 1;

  two_two[0] = 2;
  two_two[1] = 2;

  float ooF, otF, toF, ttF, totalF;
  ooF = hist->GetFrequency(one_one);
  otF = hist->GetFrequency(one_two);
  toF = hist->GetFrequency(two_one);
  ttF = hist->GetFrequency(two_two);
  totalF = hist->GetTotalFrequency();

  std::cout << ooF << std::endl;
  std::cout << otF << std::endl;
  std::cout << toF << std::endl;
  std::cout << ttF << std::endl;
  std::cout << totalF << std::endl;
  */

}
