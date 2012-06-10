#include <iostream>

#include "itkScalarImageToLocalHaralickTextureFeaturesFilter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"

const unsigned int W = 50;
const unsigned int H = 50;
const unsigned int D = 2;

typedef itk::Image<unsigned char, D> ImageType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIterator;
typedef typename itk::Statistics::ScalarImageToLocalHaralickTextureFeaturesFilter< ImageType, double  > ScalarImageToLocalHaralickTextureFeaturesFilter;
typedef typename ScalarImageToLocalHaralickTextureFeaturesFilter::OutputImageType HaralickImageType;

int main(int argc, char **argv)
{
  /*
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
  */
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName("img.bmp");
  reader->Update();

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilter;

  typename RescaleFilter::Pointer rescaler = RescaleFilter::New();
  rescaler->SetInput(reader->GetOutput());
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(4);

  typename ScalarImageToLocalHaralickTextureFeaturesFilter::Pointer haralickComputer = ScalarImageToLocalHaralickTextureFeaturesFilter::New();
  typename ScalarImageToLocalHaralickTextureFeaturesFilter::RadiusType windowRadius; windowRadius.Fill(3);
  haralickComputer->SetInput(rescaler->GetOutput());
  haralickComputer->SetWindowRadius(windowRadius);
  haralickComputer->SetNumberOfBinsPerAxis(4);

  typename ScalarImageToLocalHaralickTextureFeaturesFilter::OffsetType offset1 = {{0, 1}};
  typename ScalarImageToLocalHaralickTextureFeaturesFilter::OffsetVectorType::Pointer offsetV = ScalarImageToLocalHaralickTextureFeaturesFilter::OffsetVectorType::New();
  offsetV->push_back(offset1);
  haralickComputer->SetOffsets(offsetV);

  haralickComputer->Update();

  typename HaralickImageType::PixelType pixel;
  typename HaralickImageType::IndexType index;
  itk::ImageRegionIteratorWithIndex< HaralickImageType > outputIterator(haralickComputer->GetOutput(),haralickComputer->GetOutput()->GetRequestedRegion());
  outputIterator.GoToBegin();
  while(!outputIterator.IsAtEnd())
  {
    index = outputIterator.GetIndex();
    pixel = outputIterator.Get();
    std::cout << index << pixel << std::endl;
    ++outputIterator;
  }

  return 0;
}

