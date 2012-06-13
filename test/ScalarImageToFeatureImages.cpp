#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarImageToLocalHaralickTextureFeaturesFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkCastImageFilter.h"

const unsigned int D = 2;
const unsigned int PosterizationLevel = 4;

typedef itk::Image< unsigned char, D > InOutImage;
typedef itk::ImageFileReader< InOutImage > Reader;
typedef itk::RescaleIntensityImageFilter< InOutImage, InOutImage > ImageRescaler;
typedef typename itk::Statistics::ScalarImageToLocalHaralickTextureFeaturesFilter< InOutImage, double  > HaralickFeaturesFilter;
typedef typename HaralickFeaturesFilter::OutputImageType FeaturesImage;
typedef itk::Image< double, D > FeatureImage;
typedef itk::VectorIndexSelectionCastImageFilter< FeaturesImage, FeatureImage > IndexSelectionFilter;
typedef itk::RescaleIntensityImageFilter<FeatureImage, FeatureImage> FeatureImageRescaler;
typedef itk::CastImageFilter< FeatureImage, InOutImage > CastFilter;
typedef itk::ImageFileWriter<InOutImage> Writer;

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
  rescaler->SetOutputMaximum(PosterizationLevel);

  HaralickFeaturesFilter::Pointer haralickComputer = HaralickFeaturesFilter::New();
  HaralickFeaturesFilter::RadiusType windowRadius; windowRadius.Fill(3);
  haralickComputer->SetInput(rescaler->GetOutput());
  haralickComputer->SetWindowRadius(windowRadius);
  haralickComputer->SetNumberOfBinsPerAxis(PosterizationLevel);

  HaralickFeaturesFilter::OffsetType offset1 = {{0, 1}};
  HaralickFeaturesFilter::OffsetVectorType::Pointer offsets = HaralickFeaturesFilter::OffsetVectorType::New();
  offsets->push_back(offset1);
  haralickComputer->SetOffsets(offsets);

  IndexSelectionFilter::Pointer indexSelector = IndexSelectionFilter::New();
  indexSelector->SetInput(haralickComputer->GetOutput());

  FeatureImageRescaler::Pointer featureImageRescaler = FeatureImageRescaler::New();
  featureImageRescaler->SetInput(indexSelector->GetOutput());
  featureImageRescaler->SetOutputMinimum(0);
  featureImageRescaler->SetOutputMaximum(255);

  CastFilter::Pointer castFilter = CastFilter::New();
  castFilter->SetInput(featureImageRescaler->GetOutput());

  Writer::Pointer writer = Writer::New();
  writer->SetInput(castFilter->GetOutput());


  indexSelector->SetIndex(0);
  writer->SetFileName("energy.png");
  writer->Update();

  indexSelector->SetIndex(1);
  writer->SetFileName("entropy.png");
  writer->Update();

  indexSelector->SetIndex(2);
  writer->SetFileName("correlation.png");
  writer->Update();

  indexSelector->SetIndex(3);
  writer->SetFileName("inverse_difference_moment.png");
  writer->Update();

  indexSelector->SetIndex(4);
  writer->SetFileName("inertia.png");
  writer->Update();
  
  indexSelector->SetIndex(5);
  writer->SetFileName("cluster_shade.png");
  writer->Update();

  indexSelector->SetIndex(6);
  writer->SetFileName("cluster_prominence.png");
  writer->Update();

  indexSelector->SetIndex(7);
  writer->SetFileName("haralick_correlation.png");
  writer->Update();

  return 0;
}

