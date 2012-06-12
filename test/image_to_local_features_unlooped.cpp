#include <iostream>

#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarImageToHaralickTextureFeaturesFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVariableLengthVector.h"

const unsigned int D = 2;
const unsigned int PosterizationLevel = 4;

typedef itk::Image<unsigned char, D> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilter;

typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIterator;
typedef typename itk::Statistics::ScalarImageToHaralickTextureFeaturesFilter< ImageType, double  > ScalarImageToHaralickTextureFeaturesFilter;

typedef ::itk::Size< D > RadiusType;

int main(int argc, char **argv)
{
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName("img.bmp");
  reader->Update();

  typename RescaleFilter::Pointer rescaler = RescaleFilter::New();
  rescaler->SetInput(reader->GetOutput());
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(PosterizationLevel);
  rescaler->Update();

  typename ScalarImageToHaralickTextureFeaturesFilter::Pointer featuresComputer = ScalarImageToHaralickTextureFeaturesFilter::New();
  featuresComputer->SetNumberOfBinsPerAxis(PosterizationLevel);
  featuresComputer->SetInput(rescaler->GetOutput());

  typename ScalarImageToHaralickTextureFeaturesFilter::OffsetType offset1 = {{0, 1}};
  typename ScalarImageToHaralickTextureFeaturesFilter::OffsetVectorType::Pointer offsetV = ScalarImageToHaralickTextureFeaturesFilter::OffsetVectorType::New();
  offsetV->push_back(offset1);
  featuresComputer->SetOffsets(offsetV);

  typename ScalarImageToHaralickTextureFeaturesFilter::InputImageType::RegionType windowRegion;
  typename ScalarImageToHaralickTextureFeaturesFilter::InputImageType::IndexType windowIndex;
  RadiusType windowRadius = {{ 3, 3 }};
  typename ScalarImageToHaralickTextureFeaturesFilter::InputImageType::SizeType windowSize;

  for(unsigned int i = 0; i < D; ++i)
    {
    windowSize.SetElement(i, (windowRadius.GetElement(i) << 1) + 1);
    }

  typename ImageType::RegionType requestedRegion = rescaler->GetOutput()->GetLargestPossibleRegion();

  typedef itk::VariableLengthVector<double> VariableVectorType;
  VariableVectorType features;
  features.SetSize(8);

  itk::ImageRegionConstIteratorWithIndex< ImageType > imageIterator(rescaler->GetOutput(), requestedRegion);
  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
  {
    windowIndex = imageIterator.GetIndex();
    windowIndex -= windowRadius;

    windowRegion.SetIndex(windowIndex);
    windowRegion.SetSize(windowSize);
    windowRegion.Crop(requestedRegion);

    featuresComputer->SetRegionOfInterest(windowRegion);
    featuresComputer->Update();

    features[0] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Energy);
    features[1] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Entropy);
    features[2] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Correlation);
    features[3] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::InverseDifferenceMoment);
    features[4] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::Inertia);
    features[5] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::ClusterShade);
    features[6] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::ClusterProminence);
    features[7] = featuresComputer->GetFeature(ScalarImageToHaralickTextureFeaturesFilter::HaralickFeaturesComputer::HaralickCorrelation);

    //std::cout << imageIterator.GetIndex() << features << std::endl;

    ++imageIterator;
  }

  return 0;
}
