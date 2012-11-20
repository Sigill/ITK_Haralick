#ifndef __itkScalarImageToHaralickTextureFeaturesImageFilter_h
#define __itkScalarImageToHaralickTextureFeaturesImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorImage.h"
#include "itkGreyLevelCooccurrenceMatrix.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGLCMImageCalculator.h"
#include "itkHaralickFeaturesGLCMCalculator.h"

namespace itk
{
namespace Statistics
{

/**
 * \class ScalarImageToHaralickTextureFeaturesImageFilter
 * \brief Filter used to compute local Haralick features from an image
 * and returns it as a VectorImage of the features.
 *
 * The class is templated over the type of input image (must be a fixed
 * point scalar image whose pixel range goes from 0 to an user defined
 * number of grey-levels to consider), and the type of features expected
 * (floating point type).
 */
template< typename TInputImageType, typename TFeatureType >
class ITK_EXPORT ScalarImageToHaralickTextureFeaturesImageFilter:
  public ImageToImageFilter< TInputImageType, VectorImage< TFeatureType, TInputImageType::ImageDimension> >
{
public:
  typedef TInputImageType                                            InputImageType;
  typedef TFeatureType                                               FeatureType;
  typedef VariableLengthVector< FeatureType >                        OutputPixelType;
  typedef VectorImage< TFeatureType, InputImageType::ImageDimension> OutputImageType;

  typedef ScalarImageToHaralickTextureFeaturesImageFilter            Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >      Superclass;
  typedef SmartPointer< Self >                                       Pointer;
  typedef SmartPointer< const Self >                                 ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  typedef Statistics::GreyLevelCooccurrenceMatrix< unsigned short >  GLCMType;
  typedef GLCMImageCalculator< InputImageType, GLCMType >            GLCMCalculatorType;
  typedef HaralickFeaturesGLCMCalculator < GLCMType, float >         FeaturesCalculatorType;

  typedef typename InputImageType::PixelType                         InputPixelType;
  typedef ::itk::Size< itkGetStaticConstMacro(ImageDimension) >      RadiusType;

  typedef typename GLCMCalculatorType::OffsetType                    OffsetType;
  typedef typename GLCMCalculatorType::OffsetVectorType              OffsetVectorType;

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( FeatureTypeIsFloatingPointCheck,
    ( Concept::IsFloatingPoint< FeatureType > ) );
#endif

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScalarImageToHaralickTextureFeaturesImageFilter, ImageToImageFilter)

  /** standard New() method support */
  itkNewMacro(Self)

  ScalarImageToHaralickTextureFeaturesImageFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);

    this->m_GLCMCalculator = GLCMCalculatorType::New();
    this->m_FeaturesCalculator = FeaturesCalculatorType::New();

	this->m_FeaturesCalculator->SetCooccurrenceMatrix(this->m_GLCMCalculator->GetCooccurrenceMatrix());
  }

  virtual void GenerateInputRequestedRegion() throw( InvalidRequestedRegionError )
  {
    typename InputImageType::Pointer inputPtr =
      const_cast< InputImageType * >( this->GetInput() );

    if ( !inputPtr )
      return;

    typename InputImageType::RegionType inputRequestedRegion;
    inputRequestedRegion = inputPtr->GetRequestedRegion();

    // pad the input requested region by the operator radius
    inputRequestedRegion.PadByRadius(m_WindowRadius);

    // crop the input requested region at the input's largest possible region
    if ( inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() ) )
      {
      inputPtr->SetRequestedRegion(inputRequestedRegion);
      return;
      }
    else
      {
      // Couldn't crop the region (requested region is outside the largest
      // possible region).  Throw an exception.

      // store what we tried to request (prior to trying to crop)
      inputPtr->SetRequestedRegion(inputRequestedRegion);

      // build an exception
      InvalidRequestedRegionError e(__FILE__, __LINE__);
      e.SetLocation(ITK_LOCATION);
      e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
      e.SetDataObject(inputPtr);
      throw e;
      }
  }

  void GenerateOutputInformation(void)
  {
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage
    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetNumberOfComponentsPerPixel( 8 );
  }

  void SetInput(const InputImageType *image)
  {
    this->SetNthInput( 0, const_cast< InputImageType * >( image ) );
  }

  itkSetMacro(WindowRadius, RadiusType);
  itkGetConstMacro(WindowRadius, RadiusType);

  inline void SetOffset(const OffsetType offset)
  {
    this->m_GLCMCalculator->SetOffset(offset);
  }

  inline void SetOffsets(const OffsetVectorType * offsets)
  {
    this->m_GLCMCalculator->SetOffsets(offsets);
  }

  inline const OffsetVectorType* GetOffsets()
  {
    return this->m_GLCMCalculator->GetOffsets();
  }

  inline void SetNumberOfBinsPerAxis(const unsigned char size)
  {
    this->m_GLCMCalculator->SetMatrixSize(size);
  }

  inline const unsigned char GetNumberOfBinsPerAxis()
  {
    return this->m_GLCMCalculator->GetMatrixSize();
  }

  void GenerateData(void)
  {
    itkDebugMacro( << "GenerateData() called");
    const InputImageType *input = this->GetInput();
    OutputImageType *output = this->GetOutput();

    this->m_GLCMCalculator->SetImage(input);

    typename InputImageType::RegionType windowRegion, requestedRegion = output->GetRequestedRegion();
    typename InputImageType::IndexType windowIndex;
    typename InputImageType::SizeType windowSize;

    for(unsigned int i = 0; i < itkGetStaticConstMacro(ImageDimension); ++i)
      {
      windowSize.SetElement(i, (m_WindowRadius.GetElement(i) << 1) + 1);
      }

    OutputPixelType features;

    output->SetBufferedRegion( output->GetRequestedRegion() );
    output->Allocate();

    ImageRegionConstIteratorWithIndex< InputImageType > imageIterator(input, requestedRegion);
    ImageRegionIteratorWithIndex< OutputImageType > outputIterator(output, requestedRegion);
    imageIterator.GoToBegin();
    outputIterator.GoToBegin();
    while(!imageIterator.IsAtEnd())
    {
      windowIndex = imageIterator.GetIndex();
      windowIndex -= m_WindowRadius;

      windowRegion.SetIndex(windowIndex);
      windowRegion.SetSize(windowSize);
      windowRegion.Crop(requestedRegion);

      itkDebugMacro( << "Processing Region: " << std::endl << windowRegion);
	  this->m_GLCMCalculator->ResetMatrix();
      this->m_GLCMCalculator->SetRegion(windowRegion);
	  this->m_GLCMCalculator->Compute();
      this->m_FeaturesCalculator->Compute();

      outputIterator.Set(this->m_FeaturesCalculator->GetFeatures());

      ++imageIterator;
      ++outputIterator;
    }
  }

private:
  RadiusType                               m_WindowRadius;
  typename GLCMCalculatorType::Pointer     m_GLCMCalculator;
  typename FeaturesCalculatorType::Pointer m_FeaturesCalculator;
};

} // End of namespace Statistics
} // End of namespace itk

#endif /* __itkScalarImageToHaralickTextureFeaturesImageFilter_h */

