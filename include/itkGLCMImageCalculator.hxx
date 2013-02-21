#ifndef ITKGLCMIMAGECALCULATOR_HXX
#define ITKGLCMCIMAGEALCULATOR_HXX

#include "itkGLCMImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Statistics
{

template< typename TInputImage, typename TGLCMType >
GLCMImageCalculator< TInputImage, TGLCMType >
::GLCMImageCalculator():
  m_Image(TInputImage::New()),
  m_RegionSetByUser(false),
  m_CooccurrenceMatrix(GLCMType::New())
{}

template< typename TInputImage, typename TGLCMType >
void
GLCMImageCalculator< TInputImage, TGLCMType >
::Compute(void)
{
  if(m_MatrixSize != m_CooccurrenceMatrix->GetSize())
    {
    this->m_CooccurrenceMatrix->SetSize(m_MatrixSize);
    }

  if ( !m_RegionSetByUser )
    {
    m_Region = m_Image->GetRequestedRegion();
    }

  if(!m_Image->GetBufferedRegion().IsInside(m_Region))
    {
     itkExceptionMacro( << "The requested region is outside of the buffered region" );
    }

  itkDebugMacro( << "Processing region: " << m_Region );

  itk::ImageRegionConstIteratorWithIndex< ImageType > iterator(m_Image, m_Region);
  PixelType centerPixelIntensity, offsetPixelIntensity;
  typename ImageType::IndexType centerPixelIndex, offsetPixelIndex;
  typename OffsetVectorType::ConstIterator off_it, off_it_begin = m_Offsets->Begin(), off_it_end = m_Offsets->End();

  iterator.GoToBegin();
  while(!iterator.IsAtEnd())
  {
    centerPixelIndex = iterator.GetIndex();

    itkDebugMacro( << "Processing pixel at " << centerPixelIndex );

    centerPixelIntensity = m_Image->GetPixel(centerPixelIndex);
    // don't put a pixel in the histogram if the value is out of bound
    if ( centerPixelIntensity >= 0 && centerPixelIntensity < m_MatrixSize )
      {
      for ( off_it = off_it_begin; off_it != off_it_end; ++off_it )
        {
        offsetPixelIndex = centerPixelIndex + off_it.Value();

        itkDebugMacro( << "\tProcessing pixel at " << offsetPixelIndex );

        if(m_Region.IsInside(offsetPixelIndex))
          {
          offsetPixelIntensity = m_Image->GetPixel(offsetPixelIndex);

          // don't put a pixel in the histogram if the value is out-of-bounds.
          if ( offsetPixelIntensity >= 0 && offsetPixelIntensity < m_MatrixSize )
            {
            m_CooccurrenceMatrix->IncrementCounter(centerPixelIntensity, offsetPixelIntensity);
            m_CooccurrenceMatrix->IncrementCounter(offsetPixelIntensity, centerPixelIntensity);
            }
          else 
            {
              itkExceptionMacro( << "Value (" << offsetPixelIntensity << ") out of bound [0; " << m_MatrixSize << "[" );
            }
          }
        }
      }
      else 
      {
        itkExceptionMacro( << "Value (" << centerPixelIntensity << ") out of bound [0; " << m_MatrixSize << "[" );
      }

    ++iterator;
  }
}

template< typename TInputImage, typename TGLCMType >
void
GLCMImageCalculator< TInputImage, TGLCMType >
::SetRegion(const RegionType & region)
{
  m_Region = region;
  m_RegionSetByUser = true;

  this->Modified();
}

template< typename TInputImage, typename TGLCMType >
void
GLCMImageCalculator< TInputImage, TGLCMType >
::ResetMatrix()
{
  m_CooccurrenceMatrix->SetToZero();
}

template< typename TInputImage, typename TGLCMType >
void
GLCMImageCalculator< TInputImage, TGLCMType >
::SetOffsets(const OffsetVectorType * os)
{
  itkDebugMacro("setting offsets to " << os);
  this->m_Offsets = os;
  this->Modified();
}

template< typename TInputImage, typename TGLCMType >
void
GLCMImageCalculator< TInputImage, TGLCMType >
::SetOffset(const OffsetType offset)
{
  OffsetVectorPointer offsetVector = OffsetVectorType::New();

  offsetVector->push_back(offset);
  this->SetOffsets(offsetVector);

  this->Modified();
}

template< typename TInputImage, typename TGLCMType >
void
GLCMImageCalculator< TInputImage, TGLCMType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Image: " << std::endl;
  m_Image->Print( os, indent.GetNextIndent() );
  os << indent << "Region: " << std::endl;
  m_Region.Print( os, indent.GetNextIndent() );
  os << indent << "Region set by User: " << m_RegionSetByUser << std::endl;
}

} // end namespace Statistics
} // end namespace itk

#endif /* ITKGLCMIMAGECALCULATOR_HXX */
