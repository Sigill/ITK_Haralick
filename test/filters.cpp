#include <iostream>

#include <itkImageRegionIteratorWithIndex.h>

#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkCoocurrenceMatrixNormalizerFilter.h"
#include "itkCoocurrenceMatrixToHaralickTextureFeaturesFilter.h"

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
      /*
      if(ind[1] % 3 == 0)
        imageIt.Set(2);
      else 
        imageIt.Set(1);
      */
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

  const Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType* cooc = coocurrenceMatrixComputer->GetOutput();

  itk::Statistics::CoocurrenceMatrixNormalizerFilter< Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType, float >::Pointer normalizerFilter = 
    itk::Statistics::CoocurrenceMatrixNormalizerFilter< Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType, float >::New();

  normalizerFilter->SetInput(coocurrenceMatrixComputer->GetOutput());

  image->Modified();
  coocurrenceMatrixComputer->Update();
  normalizerFilter->Update();

  typedef typename itk::Statistics::CoocurrenceMatrixNormalizerFilter< Image2DCoocurrenceMatrixComputer::CoocurrenceMatrixType, float >::NormalizedCoocurrenceMatrixType NormalizedCoocurrenceMatrix;
  typedef typename itk::Statistics::CoocurrenceMatrixToHaralickTextureFeaturesFilter< NormalizedCoocurrenceMatrix > HaralickFeaturesComputer;
  HaralickFeaturesComputer::Pointer haralickFeaturesComputer = HaralickFeaturesComputer::New();
  haralickFeaturesComputer->SetInput(normalizerFilter->GetOutput());

  for(int i = 0; i < 1000; ++i) 
  {
    image->Modified();
    normalizerFilter->Update();
    haralickFeaturesComputer->Update();
  }

  typename NormalizedCoocurrenceMatrix::ConstPointer ncooc = normalizerFilter->GetOutput();
  typename NormalizedCoocurrenceMatrix::ConstIterator coocIt = ncooc->Begin(), coocBegin = ncooc->Begin(), coocEnd = ncooc->End();

  unsigned int i1, i2;
  while(coocIt != coocEnd)
  {
    //std::cout << coocIt.GetIndex() << " -> " << coocIt.GetFrequency() << std::endl;
    ncooc->GetIndexes((unsigned int)(coocIt - coocBegin), &i1, &i2);
    std::cout << "(" << i1 << "; " << i2 << ")" << " -> " << (*coocIt) << std::endl;
    ++coocIt;
  }
  std::cout << "Total Frequency : " << coocurrenceMatrixComputer->GetOutput()->GetTotalCount() << std::endl;


  std::cout << "Energy: " << haralickFeaturesComputer->GetEnergy() << std::endl;
  std::cout << "Entropy: " << haralickFeaturesComputer->GetEntropy() << std::endl;
  std::cout << "Correlation: " << haralickFeaturesComputer->GetCorrelation() << std::endl;
  std::cout << "InverseDifferenceMoment: " << haralickFeaturesComputer->GetInverseDifferenceMoment() << std::endl;
  std::cout << "Inertia: " << haralickFeaturesComputer->GetInertia() << std::endl;
  std::cout << "ClusterShade: " << haralickFeaturesComputer->GetClusterShade() << std::endl;
  std::cout << "ClusterProminence: " << haralickFeaturesComputer->GetClusterProminence() << std::endl;
  std::cout << "HaralickCorrelation: " << haralickFeaturesComputer->GetHaralickCorrelation() << std::endl;

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
