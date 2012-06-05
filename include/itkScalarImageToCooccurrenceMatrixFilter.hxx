/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkScalarImageToCooccurrenceMatrixFilter_hxx
#define __itkScalarImageToCooccurrenceMatrixFilter_hxx

#include "itkScalarImageToCooccurrenceMatrixFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "vnl/vnl_math.h"

namespace itk
{
namespace Statistics
{
template< class TImageType >
ScalarImageToCooccurrenceMatrixFilter< TImageType >::ScalarImageToCooccurrenceMatrixFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(1);

  this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );

  //mask inside pixel value
  this->m_InsidePixelValue = NumericTraits< PixelType >::One;

  this->m_NumberOfBinsPerAxis = DefaultBinsPerAxis;
  this->m_Normalize = false;
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >
::SetOffset(const OffsetType offset)
{
  OffsetVectorPointer offsetVector = OffsetVector::New();

  offsetVector->push_back(offset);
  this->SetOffsets(offsetVector);
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >
::SetInput(const ImageType *image)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                    const_cast< ImageType * >( image ) );
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >
::SetMaskImage(const ImageType *image)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 1,
                                    const_cast< ImageType * >( image ) );
}

template< class TImageType >
const TImageType *
ScalarImageToCooccurrenceMatrixFilter< TImageType >
::GetInput() const
{
  return static_cast< const ImageType * >( this->GetPrimaryInput() );
}

template< class TImageType >
const TImageType *
ScalarImageToCooccurrenceMatrixFilter< TImageType >
::GetMaskImage() const
{
  return static_cast< const ImageType * >( this->ProcessObject::GetInput(1) );
}

template< class TImageType >
const typename ScalarImageToCooccurrenceMatrixFilter< TImageType >::CoocurrenceMatrixType *
ScalarImageToCooccurrenceMatrixFilter< TImageType >
::GetOutput() const
{
  const CoocurrenceMatrixType *output =
    static_cast< const CoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  return output;
}

template< class TImageType >
typename ScalarImageToCooccurrenceMatrixFilter< TImageType >::DataObjectPointer
ScalarImageToCooccurrenceMatrixFilter< TImageType >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
{
  typename CoocurrenceMatrixType::Pointer output = CoocurrenceMatrixType::New();
  return static_cast< DataObject * >( output.GetPointer() );
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >::GenerateData(void)
{
  CoocurrenceMatrixType *output =
    static_cast< CoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  const ImageType *input = this->GetInput();

  // At this point input must be non-NULL because the ProcessObject
  // checks the number of required input to be non-NULL pointers before
  // calling this GenerateData() method.

  // Next, find the minimum radius that encloses all the offsets.
  unsigned int minRadius = 0;
  typename OffsetVector::ConstIterator offsets;
  for ( offsets = m_Offsets->Begin(); offsets != m_Offsets->End(); offsets++ )
    {
    for ( unsigned int i = 0; i < offsets.Value().GetOffsetDimension(); i++ )
      {
      unsigned int distance = vnl_math_abs(offsets.Value()[i]);
      if ( distance > minRadius )
        {
        minRadius = distance;
        }
      }
    }

  RadiusType radius;
  radius.Fill(minRadius);

  const ImageType *maskImage = NULL;

  // Check if a mask image has been provided
  //
  if ( this->GetNumberOfIndexedInputs() > 1 )
    {
    maskImage = this->GetMaskImage();
    }

  // Now fill in the histogram
  if ( maskImage != NULL )
    {
    this->FillCoocurrenceMatrixWithMask(radius, input->GetRequestedRegion(), maskImage);
    }
  else
    {
    this->FillCoocurrenceMatrix( radius, input->GetRequestedRegion() );
    }

  // Normalizse the histogram if requested
  if ( m_Normalize )
    {
    this->NormalizeCoocurrenceMatrix();
    }
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >::FillCoocurrenceMatrix(RadiusType radius,
                                                                                     RegionType region)
{
  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  CoocurrenceMatrixType *output =
    static_cast< CoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  typedef ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  NeighborhoodIteratorType neighborIt;
  neighborIt = NeighborhoodIteratorType(radius, input, region);

  for ( neighborIt.GoToBegin(); !neighborIt.IsAtEnd(); ++neighborIt )
    {
    const PixelType centerPixelIntensity = neighborIt.GetCenterPixel();
    if ( centerPixelIntensity < 0 || centerPixelIntensity >= m_NumberOfBinsPerAxis )
      {
      continue; // don't put a pixel in the histogram if the value
                // is out-of-bounds.
      }

    typename OffsetVector::ConstIterator offsets;
    for ( offsets = m_Offsets->Begin(); offsets != m_Offsets->End(); offsets++ )
      {
      bool            pixelInBounds;
      const PixelType pixelIntensity =
        neighborIt.GetPixel(offsets.Value(), pixelInBounds);

      if ( !pixelInBounds )
        {
        continue; // don't put a pixel in the histogram if it's out-of-bounds.
        }

      if ( pixelIntensity < 0 || pixelIntensity >= m_NumberOfBinsPerAxis)
        {
        continue; // don't put a pixel in the histogram if the value
                  // is out-of-bounds.
        }

      // Now make both possible co-occurrence combinations and increment the
      // histogram with them.

      output->IncrementFrequency(centerPixelIntensity, pixelIntensity);
      output->IncrementFrequency(pixelIntensity, centerPixelIntensity);
      }
    }
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >::FillCoocurrenceMatrixWithMask(RadiusType radius,
                                                                                             RegionType region,
                                                                                             const ImageType *maskImage)
{
  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  CoocurrenceMatrixType *output =
    static_cast< CoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram
  typedef ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  NeighborhoodIteratorType neighborIt, maskNeighborIt;
  neighborIt = NeighborhoodIteratorType(radius, input, region);
  maskNeighborIt = NeighborhoodIteratorType(radius, maskImage, region);

  for ( neighborIt.GoToBegin(), maskNeighborIt.GoToBegin();
        !neighborIt.IsAtEnd(); ++neighborIt, ++maskNeighborIt )
    {
    if ( maskNeighborIt.GetCenterPixel() != m_InsidePixelValue )
      {
      continue; // Go to the next loop if we're not in the mask
      }

    const PixelType centerPixelIntensity = neighborIt.GetCenterPixel();

    if ( centerPixelIntensity < 0 || centerPixelIntensity >= m_NumberOfBinsPerAxis )
      {
      continue; // don't put a pixel in the histogram if the value
                // is out-of-bounds.
      }

    typename OffsetVector::ConstIterator offsets;
    for ( offsets = this->GetOffsets()->Begin(); offsets != this->GetOffsets()->End(); offsets++ )
      {
      if ( maskNeighborIt.GetPixel( offsets.Value() ) != m_InsidePixelValue )
        {
        continue; // Go to the next loop if we're not in the mask
        }

      bool            pixelInBounds;
      const PixelType pixelIntensity =
        neighborIt.GetPixel(offsets.Value(), pixelInBounds);

      if ( !pixelInBounds )
        {
        continue; // don't put a pixel in the histogram if it's out-of-bounds.
        }

      if ( pixelIntensity < 0 || pixelIntensity > m_NumberOfBinsPerAxis )
        {
        continue; // don't put a pixel in the histogram if the value
                  // is out-of-bounds.
        }

      // Now make both possible co-occurrence combinations and increment the
      // histogram with them.

      output->IncrementFrequency(centerPixelIntensity, pixelIntensity);
      output->IncrementFrequency(pixelIntensity, centerPixelIntensity);
      }
    }
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >::NormalizeCoocurrenceMatrix(void)
{
  CoocurrenceMatrixType *output =
    static_cast< CoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  output->Normalize();
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >::PrintSelf(std::ostream & os,
                                                                                 Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Offsets: " << this->GetOffsets() << std::endl;
  os << indent << "NumberOfBinsPerAxis: " << this->GetNumberOfBinsPerAxis() << std::endl;
  os << indent << "Normalize: " << this->GetNormalize() << std::endl;
  os << indent << "InsidePixelValue: " << this->GetInsidePixelValue() << std::endl;
}
} // end of namespace Statistics
} // end of namespace itk

#endif
