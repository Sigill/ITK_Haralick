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
#include "itkImageRegionConstIteratorWithIndex.h"
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

  output->Reset();

  const ImageType *input = this->GetInput();

  // At this point input must be non-NULL because the ProcessObject
  // checks the number of required input to be non-NULL pointers before
  // calling this GenerateData() method.

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
    this->FillCoocurrenceMatrixWithMask(maskImage);
    }
  else
    {
    this->FillCoocurrenceMatrix();
    }

  // Normalizse the histogram if requested
  if ( m_Normalize )
    {
    this->NormalizeCoocurrenceMatrix();
    }
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >::FillCoocurrenceMatrix(void)
{
  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  CoocurrenceMatrixType *output =
    static_cast< CoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  itk::ImageRegionConstIteratorWithIndex< ImageType > iterator(input, m_RegionOfInterest);
  PixelType centerPixelIntensity, offsetPixelIntensity;
  typename ImageType::IndexType centerPixelIndex, offsetPixelIndex;
  typename OffsetVector::ConstIterator off_it, off_it_begin = m_Offsets->Begin(), off_it_end = m_Offsets->End();

  iterator.GoToBegin();
  while(!iterator.IsAtEnd())
  {
    centerPixelIndex = iterator.GetIndex();
    centerPixelIntensity = input->GetPixel(centerPixelIndex);
    if ( centerPixelIntensity < 0 || centerPixelIntensity >= m_NumberOfBinsPerAxis )
      {
      continue; // don't put a pixel in the histogram if the value
                // is out-of-bounds.
      }

    for ( off_it = off_it_begin; off_it != off_it_end; ++off_it )
      {
        offsetPixelIndex = centerPixelIndex + off_it.Value();

        if(m_RegionOfInterest.IsInside(offsetPixelIndex))
          {
          offsetPixelIntensity = input->GetPixel(offsetPixelIndex);

          if ( offsetPixelIntensity < 0 || offsetPixelIntensity >= m_NumberOfBinsPerAxis )
            {
            continue; // don't put a pixel in the histogram if the value
                      // is out-of-bounds.
            }

          output->IncrementFrequency(centerPixelIntensity, offsetPixelIntensity);
          output->IncrementFrequency(offsetPixelIntensity, centerPixelIntensity);
          }
      }
    ++iterator;
  }
}

template< class TImageType >
void
ScalarImageToCooccurrenceMatrixFilter< TImageType >::FillCoocurrenceMatrixWithMask(const ImageType *maskImage)
{
  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  CoocurrenceMatrixType *output =
    static_cast< CoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  itk::ImageRegionConstIteratorWithIndex< ImageType > iterator(input, m_RegionOfInterest);
  PixelType centerPixelIntensity, offsetPixelIntensity;
  typename ImageType::IndexType centerPixelIndex, offsetPixelIndex;
  typename OffsetVector::ConstIterator off_it, off_it_begin = m_Offsets->Begin(), off_it_end = m_Offsets->End();

  iterator.GoToBegin();
  while(!iterator.IsAtEnd())
  {
    centerPixelIndex = iterator.GetIndex();
    if ( maskImage->GetPixel(centerPixelIndex) != m_InsidePixelValue )
      {
      continue; // Go to the next loop if we're not in the mask
      }

    centerPixelIntensity = input->GetPixel(centerPixelIndex);
    if ( centerPixelIntensity < 0 || centerPixelIntensity >= m_NumberOfBinsPerAxis )
      {
      continue; // don't put a pixel in the histogram if the value
                // is out-of-bounds.
      }

    for ( off_it = off_it_begin; off_it != off_it_end; ++off_it )
      {
        offsetPixelIndex = centerPixelIndex + off_it.Value();
        if(m_RegionOfInterest.IsInside(offsetPixelIndex))
          {
          if ( maskImage->GetPixel(offsetPixelIndex) != m_InsidePixelValue )
            {
            continue; // Go to the next loop if we're not in the mask
            }

          offsetPixelIntensity = input->GetPixel(offsetPixelIndex);

          if ( offsetPixelIntensity < 0 || offsetPixelIntensity >= m_NumberOfBinsPerAxis )
            {
            continue; // don't put a pixel in the histogram if the value
                      // is out-of-bounds.
            }

          output->IncrementFrequency(centerPixelIntensity, offsetPixelIntensity);
          output->IncrementFrequency(offsetPixelIntensity, centerPixelIntensity);
          }
      }
    ++iterator;
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
ScalarImageToCooccurrenceMatrixFilter< TImageType >::ComputeOffsetsMinRadius(void)
{
  unsigned int distance, i;
  typename OffsetVector::ConstIterator off_it;

  m_OffsetsMinRadius.Fill(0);

  for ( off_it = m_Offsets->Begin(); off_it != m_Offsets->End(); off_it++ )
  {
    for ( i = 0; i < ::itk::GetImageDimension< ImageType >::ImageDimension; i++ )
    {
      distance = vnl_math_abs(off_it.Value()[i]);
      if ( distance > m_OffsetsMinRadius[i] )
      {
        m_OffsetsMinRadius[i] = distance;
      }
    }
  }
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
