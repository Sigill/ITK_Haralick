#include <iostream>
#include <itkVectorContainer.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>

const unsigned int W = 16;
const unsigned int H = 16;
const unsigned int D = 16;
const unsigned int Dim = 3;

typedef unsigned char PixelType;
typedef itk::Image<PixelType, Dim> ImageType;
typedef itk::ImageRegionConstIteratorWithIndex< ImageType > ConstIteratorWidx;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorWidx;
typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIterator;

typedef ImageType::OffsetType OffsetType;
typedef itk::VectorContainer< unsigned char, OffsetType > OffsetVector;
typedef typename OffsetVector::Pointer OffsetVectorPointer;

typedef typename NeighborhoodIterator::NeighborIndexType NeighborIndexType;
typedef itk::VectorContainer< unsigned char, NeighborIndexType > NeighborIndexVector;
typedef NeighborIndexVector::Pointer NeighborIndexVectorPointer;

int main(int argc, char **argv)
{
  ImageType::Pointer image = ImageType::New();

  ImageType::RegionType imageRegion;

  {
    ImageType::IndexType index = {{ 0, 0, 0 }};
    ImageType::SizeType imageSize = {{ W, H, D }};

    imageRegion.SetIndex( index );
    imageRegion.SetSize( imageSize );
  }

  image->SetRegions( imageRegion );
  image->Allocate();

  {
    IteratorWidx imageIt( image, imageRegion );
    ImageType::IndexType ind;
    for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt)
    {
      ind = imageIt.GetIndex();
      imageIt.Set(ind[1] % 2 + 1);
    }
  }

  OffsetVectorPointer offsets = OffsetVector::New();
  {
    OffsetType off;
    off[0] = 2; off[1] = 0; off[2] = 0;
    offsets->push_back(off);
    off[0] = 0; off[1] = 2;
    offsets->push_back(off);
  }

  NeighborhoodIterator::RadiusType offsetsRadius; offsetsRadius.Fill(0);
  {
    OffsetVector::ConstIterator off_it;
    for ( off_it = offsets->Begin(); off_it != offsets->End(); off_it++ )
    {
      for ( unsigned int i = 0; i < Dim; i++ )
      {
        unsigned int distance = vnl_math_abs(off_it.Value()[i]);
        if ( distance > offsetsRadius[i] )
        {
          offsetsRadius[i] = distance;
        }
      }
    }
  }

  std::cout << "Offsets minimal radius: " << offsetsRadius << std::endl;

  NeighborIndexVectorPointer offsetsIndexes = NeighborIndexVector::New();
  {
    NeighborhoodIterator::NeighborhoodType neighborhood;
    neighborhood.SetRadius(offsetsRadius);

    OffsetVector::ConstIterator off_it;
    for ( off_it = offsets->Begin(); off_it != offsets->End(); off_it++ )
    {
      offsetsIndexes->push_back(neighborhood.GetNeighborhoodIndex(off_it.Value()));
      std::cout << neighborhood.GetNeighborhoodIndex(off_it.Value()) << std::endl;
    }
  }

  ImageType::RegionType windowRegion;
  ImageType::IndexType windowIndex;
  ImageType::OffsetType windowRadius; windowRadius.Fill(3);
  ImageType::SizeType windowSize; windowSize.Fill((3 << 1) + 1);

  std::cout << "Window size: " << windowSize << std::endl;

  ImageType::IndexType pixelIndex;

  NeighborIndexVector::ConstIterator off_it;
  PixelType centerPixelIntensity, offsetPixelIntensity;
  NeighborIndexType offsetPixelNeighborIndex;
  ImageType::IndexType offsetPixelIndex;

  bool pixelInBounds;

  ConstIteratorWidx iit(image, imageRegion);
  iit.GoToBegin();
  while(!iit.IsAtEnd()) {
    windowIndex = iit.GetIndex();
    //std::cout << windowIndex << std::endl;
    windowIndex -= windowRadius;

    windowRegion.SetIndex(windowIndex);
    windowRegion.SetSize(windowSize);
    windowRegion.Crop(imageRegion);

    NeighborhoodIterator nit(offsetsRadius, image, windowRegion);
    nit.NeedToUseBoundaryConditionOff();
    nit.GoToBegin();
    while(!nit.IsAtEnd()) {
      centerPixelIntensity = nit.GetCenterPixel();
      //std::cout << "\t" << nit.GetIndex() << std::endl;
      for ( off_it = offsetsIndexes->Begin(); off_it != offsetsIndexes->End(); off_it++ )
      {
        offsetPixelNeighborIndex = off_it.Value();
        offsetPixelIndex = nit.GetIndex(offsetPixelNeighborIndex);

        if(!windowRegion.IsInside(offsetPixelIndex))
        {
          //std::cout << "\t\t" << offsetPixelIndex << std::endl;
          continue;
        }
        offsetPixelIntensity = image->GetPixel(offsetPixelIndex);
      }

      ++nit;
    }

    ++iit;
  }

  return 0;
}
