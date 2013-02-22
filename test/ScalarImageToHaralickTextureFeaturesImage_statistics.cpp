#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarImageToHaralickTextureFeaturesImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include <vector>

const unsigned int D = 2;
const unsigned int PosterizationLevel = 4;

typedef itk::Image< unsigned char, D > InOutImage;
typedef itk::ImageFileReader< InOutImage > Reader;
typedef itk::RescaleIntensityImageFilter< InOutImage, InOutImage > ImageRescaler;
typedef typename itk::Statistics::ScalarImageToHaralickTextureFeaturesImageFilter< InOutImage, double  > ScalarImageToHaralickTextureFeaturesImageFilter;
typedef typename ScalarImageToHaralickTextureFeaturesImageFilter::OutputImageType FeaturesImage;
typedef itk::Image< double, D > FeatureImage;
typedef itk::VectorIndexSelectionCastImageFilter< FeaturesImage, FeatureImage > IndexSelectionFilter;
typedef itk::RescaleIntensityImageFilter<FeatureImage, FeatureImage> FeatureImageRescaler;
typedef itk::CastImageFilter< FeatureImage, InOutImage > CastFilter;
typedef itk::ImageFileWriter<InOutImage> Writer;
typedef itk::StatisticsImageFilter< FeatureImage > StatisticsComputer;

int main(int argc, char **argv)
{
  if(argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " <image>" << std::endl << std::endl;
    exit(0);
  }

  Reader::Pointer reader = Reader::New();
  reader->SetFileName(argv[1]);

  ImageRescaler::Pointer rescaler = ImageRescaler::New();
  rescaler->SetInput(reader->GetOutput());
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(PosterizationLevel - 1);

  ScalarImageToHaralickTextureFeaturesImageFilter::Pointer haralickComputer = ScalarImageToHaralickTextureFeaturesImageFilter::New();
  ScalarImageToHaralickTextureFeaturesImageFilter::RadiusType windowRadius; windowRadius.Fill(3);
  haralickComputer->SetInput(rescaler->GetOutput());
  haralickComputer->SetWindowRadius(windowRadius);
  haralickComputer->SetNumberOfBinsPerAxis(PosterizationLevel);

  ScalarImageToHaralickTextureFeaturesImageFilter::OffsetType offset1 = {{0, 1}};
  ScalarImageToHaralickTextureFeaturesImageFilter::OffsetVectorType::Pointer offsets = ScalarImageToHaralickTextureFeaturesImageFilter::OffsetVectorType::New();
  offsets->push_back(offset1);
  haralickComputer->SetOffsets(offsets);

  IndexSelectionFilter::Pointer indexSelector = IndexSelectionFilter::New();
  indexSelector->SetInput(haralickComputer->GetOutput());

  StatisticsComputer::Pointer statisticsComputer = StatisticsComputer::New();
  statisticsComputer->SetInput(indexSelector->GetOutput());

  FeatureImageRescaler::Pointer featureImageRescaler = FeatureImageRescaler::New();
  featureImageRescaler->SetInput(indexSelector->GetOutput());
  featureImageRescaler->SetOutputMinimum(0);
  featureImageRescaler->SetOutputMaximum(255);

  CastFilter::Pointer castFilter = CastFilter::New();
  castFilter->SetInput(featureImageRescaler->GetOutput());

  Writer::Pointer writer = Writer::New();
  writer->SetInput(castFilter->GetOutput());

  std::vector< std::string > features_name;
  features_name.push_back("angularSecondMoment");
  features_name.push_back("contrast");
  features_name.push_back("variance");
  features_name.push_back("inverse_difference_moment");
  features_name.push_back("sum_average");
  features_name.push_back("sum_variance");
  features_name.push_back("sum_entropy");
  features_name.push_back("entropy");
  features_name.push_back("differenceVariance");
  features_name.push_back("differenceEntropy");

  std::vector< std::string >::const_iterator names_it = features_name.begin(), names_end = features_name.end();
  unsigned int i = 0;
  while(names_it != names_end)
  {
    indexSelector->SetIndex(i);
    writer->SetFileName((*names_it) + ".png");
    writer->Update();

    statisticsComputer->Update();

    std::cout << (*names_it) << std::endl;
    std::cout << "\tMin: " << statisticsComputer->GetMinimum() << std::endl;
    std::cout << "\tMax: " << statisticsComputer->GetMaximum() << std::endl;
    std::cout << "\tMean: " << statisticsComputer->GetMean() << std::endl;
    std::cout << "\tStd.: " << statisticsComputer->GetSigma() << std::endl;

    ++i; ++names_it;
  }

  return 0;
}

