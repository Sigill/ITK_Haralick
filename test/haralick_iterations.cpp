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

  ImageType::OffsetType offsetsRadius; offsetsRadius.Fill(0);
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

  ImageType::RegionType windowRegion;
  ImageType::IndexType windowIndex;
  ImageType::OffsetType windowRadius; windowRadius.Fill(3);
  ImageType::SizeType windowSize; windowSize.Fill((3 << 1) + 1);

  std::cout << "Window size: " << windowSize << std::endl;

  ImageType::IndexType pixelIndex;

  OffsetVector::ConstIterator off_it, off_it_begin = offsets->Begin(), off_it_end = offsets->End();
  PixelType centerPixelIntensity, offsetPixelIntensity;
  ImageType::IndexType centerPixelIndex, offsetPixelIndex;

  ConstIteratorWidx iit(image, imageRegion);
  iit.GoToBegin();
  while(!iit.IsAtEnd()) {
    windowIndex = iit.GetIndex();
    //std::cout << windowIndex << std::endl;
    windowIndex -= windowRadius;

    windowRegion.SetIndex(windowIndex);
    windowRegion.SetSize(windowSize);
    windowRegion.Crop(imageRegion);

    ConstIteratorWidx wit(image, windowRegion);
    wit.GoToBegin();
    while(!wit.IsAtEnd())
    {
      centerPixelIndex = wit.GetIndex();

      for ( off_it = off_it_begin; off_it != off_it_end; ++off_it )
      {
        offsetPixelIndex = centerPixelIndex + off_it.Value();

        if(windowRegion.IsInside(offsetPixelIndex))
        {
          offsetPixelIntensity = image->GetPixel(offsetPixelIndex);
        }
      }
      ++wit;
    }

    ++iit;
  }

  return 0;
}
