#include <iostream>

#include "itkScalarImageToHaralickTextureFeaturesFilter.h"

#include "itkImageRegionIteratorWithIndex.h"

#include <callgrind.h>

const unsigned int W = 200;
const unsigned int H = 200;
const unsigned int D = 2;

typedef itk::Image<unsigned char, D> ImageType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIterator;
typedef typename itk::Statistics::ScalarImageToHaralickTextureFeaturesFilter< ImageType, double  > ScalarImageToHaralickTextureFeaturesFilter;

int main(int argc, char **argv)
{
  ImageType::Pointer image = ImageType::New();

  ImageType::RegionType region;

  {
    ImageType::IndexType index = {{ 0, 0 }};
    ImageType::SizeType imageSize = {{ W, H }};

    region.SetIndex( index );
    region.SetSize( imageSize );
  }

  image->SetRegions( region );
  image->Allocate();

  {
    ImageIterator imageIt( image, region );
    ImageType::IndexType ind;
    for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt)
    {
      ind = imageIt.GetIndex();
      imageIt.Set(ind[1] % 2 + 1);
    }
  }

  typename ScalarImageToHaralickTextureFeaturesFilter::Pointer featuresComputer = ScalarImageToHaralickTextureFeaturesFilter::New();

  featuresComputer->SetNumberOfBinsPerAxis(4);

  featuresComputer->SetInput(image);

  typename ScalarImageToHaralickTextureFeaturesFilter::OffsetType offset1 = {{0, 1}};
  typename ScalarImageToHaralickTextureFeaturesFilter::OffsetVectorType::Pointer offsetV = ScalarImageToHaralickTextureFeaturesFilter::OffsetVectorType::New();
  offsetV->push_back(offset1);
  featuresComputer->SetOffsets(offsetV);

  typename ScalarImageToHaralickTextureFeaturesFilter::InputImageType::RegionType roi;
  typename ScalarImageToHaralickTextureFeaturesFilter::InputImageType::IndexType rIndex = {{ 1, 0 }};
  typename ScalarImageToHaralickTextureFeaturesFilter::InputImageType::SizeType rSize = {{ 7, 7 }};
  roi.SetSize( rSize );
  roi.SetIndex( rIndex );

  //CALLGRIND_START_INSTRUMENTATION 
  featuresComputer->SetRegionOfInterest(roi);

  featuresComputer->Update();

  /*
  rIndex[1] = 1;
  roi.SetIndex(rIndex);
  featuresComputer->SetRegionOfInterest(roi);

  featuresComputer->Update();
  */
  //CALLGRIND_STOP_INSTRUMENTATION 

  std::cout << "Energy: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Energy) << std::endl;
  std::cout << "Entropy: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Entropy) << std::endl;
  std::cout << "Correlation: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Correlation) << std::endl;
  std::cout << "InverseDifferenceMoment: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::InverseDifferenceMoment) << std::endl;
  std::cout << "Inertia: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Inertia) << std::endl;
  std::cout << "ClusterShade: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::ClusterShade) << std::endl;
  std::cout << "ClusterProminence: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::ClusterProminence) << std::endl;
  std::cout << "HaralickCorrelation: " << featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::HaralickCorrelation) << std::endl;

  return 0;
}
