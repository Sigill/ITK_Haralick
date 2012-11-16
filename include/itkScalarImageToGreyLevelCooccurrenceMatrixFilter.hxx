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
#ifndef __itkScalarImageToGreyLevelCooccurrenceMatrixFilter_hxx
#define __itkScalarImageToGreyLevelCooccurrenceMatrixFilter_hxx

#include "itkScalarImageToGreyLevelCooccurrenceMatrixFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "vnl/vnl_math.h"

namespace itk
{
namespace Statistics
{
template< class TImageType >
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::ScalarImageToGreyLevelCooccurrenceMatrixFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(1);

  this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );

  //mask inside pixel value
  this->m_InsidePixelValue = NumericTraits< PixelType >::One;

  this->m_NumberOfBinsPerAxis = DefaultBinsPerAxis;
}

template< class TImageType >
void
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >
::SetOffset(const OffsetType offset)
{
  OffsetVectorPointer offsetVector = OffsetVector::New();

  offsetVector->push_back(offset);
  this->SetOffsets(offsetVector);
  this->ComputeOffsetsMinRadius();

  this->Modified();
}

template< class TImageType >
void
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >
::SetInput(const ImageType *image)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                    const_cast< ImageType * >( image ) );
}

template< class TImageType >
void
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >
::SetMaskImage(const ImageType *image)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 1,
                                    const_cast< ImageType * >( image ) );
}

template< class TImageType >
const TImageType *
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >
::GetInput() const
{
  return static_cast< const ImageType * >( this->GetPrimaryInput() );
}

template< class TImageType >
const TImageType *
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >
::GetMaskImage() const
{
  return static_cast< const ImageType * >( this->ProcessObject::GetInput(1) );
}

template< class TImageType >
const typename ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::GreyLevelCooccurrenceMatrixType *
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >
::GetOutput() const
{
  const GreyLevelCooccurrenceMatrixType *output =
    static_cast< const GreyLevelCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  return output;
}

template< class TImageType >
typename ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::DataObjectPointer
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
{
  typename GreyLevelCooccurrenceMatrixType::Pointer output = GreyLevelCooccurrenceMatrixType::New();
  return static_cast< DataObject * >( output.GetPointer() );
}

template< class TImageType >
void
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::GenerateData(void)
{
  itkDebugMacro( << "GenerateData called" );
  GreyLevelCooccurrenceMatrixType *output =
    static_cast< GreyLevelCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  output->SetToZero();

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
    this->FillGreyLevelCooccurrenceMatrixWithMask(maskImage);
    }
  else
    {
    this->FillGreyLevelCooccurrenceMatrix();
    }

  itkDebugMacro( << "Coocurrence Matrix generated" );
}

template< class TImageType >
void
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::FillGreyLevelCooccurrenceMatrix(void)
{
  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  GreyLevelCooccurrenceMatrixType *output =
    static_cast< GreyLevelCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  itkDebugMacro( << "Processing region: " << m_RegionOfInterest );

  itk::ImageRegionConstIteratorWithIndex< ImageType > iterator(input, m_RegionOfInterest);
  PixelType centerPixelIntensity, offsetPixelIntensity;
  typename ImageType::IndexType centerPixelIndex, offsetPixelIndex;
  typename OffsetVector::ConstIterator off_it, off_it_begin = m_Offsets->Begin(), off_it_end = m_Offsets->End();

  iterator.GoToBegin();
  while(!iterator.IsAtEnd())
  {
    centerPixelIndex = iterator.GetIndex();

    itkDebugMacro( << "Processing pixel at " << centerPixelIndex );

    centerPixelIntensity = input->GetPixel(centerPixelIndex);
    // don't put a pixel in the histogram if the value is out of bound
    if ( centerPixelIntensity >= 0 && centerPixelIntensity < m_NumberOfBinsPerAxis )
      {
      for ( off_it = off_it_begin; off_it != off_it_end; ++off_it )
        {
        offsetPixelIndex = centerPixelIndex + off_it.Value();

        itkDebugMacro( << "\tProcessing pixel at " << offsetPixelIndex );

        if(m_RegionOfInterest.IsInside(offsetPixelIndex))
          {
          offsetPixelIntensity = input->GetPixel(offsetPixelIndex);

          // don't put a pixel in the histogram if the value is out-of-bounds.
          if ( offsetPixelIntensity >= 0 && offsetPixelIntensity < m_NumberOfBinsPerAxis )
            {
            output->IncrementCounter(centerPixelIntensity, offsetPixelIntensity);
            output->IncrementCounter(offsetPixelIntensity, centerPixelIntensity);
            }
          else 
            {
              itkDebugMacro( << "\t\tValue (" << offsetPixelIntensity << ") out of bound [0; " << m_NumberOfBinsPerAxis << "[" );
            }
          }
        }
      }
      else 
      {
        itkDebugMacro( << "\tValue (" << centerPixelIntensity << ") out of bound [0; " << m_NumberOfBinsPerAxis << "[" );
      }

    ++iterator;
  }
}

template< class TImageType >
void
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::FillGreyLevelCooccurrenceMatrixWithMask(const ImageType *maskImage)
{
  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  GreyLevelCooccurrenceMatrixType *output =
    static_cast< GreyLevelCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

  itk::ImageRegionConstIteratorWithIndex< ImageType > iterator(input, m_RegionOfInterest);
  PixelType centerPixelIntensity, offsetPixelIntensity;
  typename ImageType::IndexType centerPixelIndex, offsetPixelIndex;
  typename OffsetVector::ConstIterator off_it, off_it_begin = m_Offsets->Begin(), off_it_end = m_Offsets->End();

  iterator.GoToBegin();
  while(!iterator.IsAtEnd())
  {
    centerPixelIndex = iterator.GetIndex();
    // Go to the next loop if we're not in the mask
    if ( maskImage->GetPixel(centerPixelIndex) == m_InsidePixelValue )
      {
      centerPixelIntensity = input->GetPixel(centerPixelIndex);

      // don't put a pixel in the histogram if the value is out-of-bounds.
      if ( centerPixelIntensity >= 0 && centerPixelIntensity < m_NumberOfBinsPerAxis )
        {
        for ( off_it = off_it_begin; off_it != off_it_end; ++off_it )
          {
          offsetPixelIndex = centerPixelIndex + off_it.Value();
          if(m_RegionOfInterest.IsInside(offsetPixelIndex))
            {
            // Go to the next loop if we're not in the mask
            if ( maskImage->GetPixel(offsetPixelIndex) == m_InsidePixelValue )
              {
              offsetPixelIntensity = input->GetPixel(offsetPixelIndex);

              // don't put a pixel in the histogram if the value is out-of-bounds.
              if ( offsetPixelIntensity >= 0 && offsetPixelIntensity < m_NumberOfBinsPerAxis )
                {
                output->IncrementCounter(centerPixelIntensity, offsetPixelIntensity);
                output->IncrementCounter(offsetPixelIntensity, centerPixelIntensity);
                }
              }
            }
          }
        }
      }
    ++iterator;
  }
}

template< class TImageType >
void
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::ComputeOffsetsMinRadius(void)
{
  unsigned int distance, i;
  typename OffsetVector::ConstIterator off_it;

  m_OffsetsMinRadius.Fill(0);

  for ( off_it = m_Offsets->Begin(); off_it != m_Offsets->End(); off_it++ )
  {
    for ( i = 0; i < ImageType::ImageDimension; i++ )
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
ScalarImageToGreyLevelCooccurrenceMatrixFilter< TImageType >::PrintSelf(std::ostream & os,
                                                                                 Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Offsets: " << this->GetOffsets() << std::endl;
  os << indent << "NumberOfBinsPerAxis: " << this->GetNumberOfBinsPerAxis() << std::endl;
  os << indent << "InsidePixelValue: " << this->GetInsidePixelValue() << std::endl;
}
} // end of namespace Statistics
} // end of namespace itk

#endif
