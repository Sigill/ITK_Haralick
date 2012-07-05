#ifndef __itkScalarImageToHaralickTextureFeaturesImageFilter_h
#define __itkScalarImageToHaralickTextureFeaturesImageFilter_h

#include "itkScalarImageToHaralickTextureFeaturesFilter.h"
#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
namespace Statistics
{

template< class TInputImageType, class TOutputPixelType >
class ITK_EXPORT ScalarImageToHaralickTextureFeaturesImageFilter:
  public ImageToImageFilter< TInputImageType, VectorImage< TOutputPixelType, ::itk::GetImageDimension< TInputImageType >::ImageDimension > >
{
public:
  /** Standard typedefs */
  typedef ScalarImageToHaralickTextureFeaturesImageFilter   Self;
  typedef TInputImageType                                   InputImageType;
  typedef TOutputPixelType                                  OutputPixelType;

  itkStaticConstMacro(ImageDimension, unsigned int, ::itk::GetImageDimension< InputImageType >::ImageDimension);

  typedef VectorImage< OutputPixelType, itkGetStaticConstMacro(ImageDimension) >   OutputImageType;

  typedef ImageToImageFilter< InputImageType, OutputImageType >   Superclass;
  typedef SmartPointer< Self >                                    Pointer;
  typedef SmartPointer< const Self >                              ConstPointer;

  typedef typename InputImageType::PixelType    InputPixelType;
  typedef ::itk::Size< itkGetStaticConstMacro(ImageDimension) >    RadiusType;

  typedef ScalarImageToHaralickTextureFeaturesFilter< InputImageType, OutputPixelType > LocalHaralickComputer;

  typedef typename LocalHaralickComputer::OffsetType         OffsetType;
  typedef typename LocalHaralickComputer::OffsetVectorType   OffsetVectorType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScalarImageToHaralickTextureFeaturesImageFilter, ImageToImageFilter);

  /** standard New() method support */
  itkNewMacro(Self);

  ScalarImageToHaralickTextureFeaturesImageFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);

    this->m_LocalHaralickComputer = LocalHaralickComputer::New();

    //this->ProcessObject::SetOutput(this->MakeOutput
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
    this->m_LocalHaralickComputer->SetOffset(offset);
  }

  inline void SetOffsets(const OffsetVectorType * offsets)
  {
    this->m_LocalHaralickComputer->SetOffsets(offsets);
  }

  inline const OffsetVectorType* GetOffsets()
  {
    return this->m_LocalHaralickComputer->GetOffset();
  }


  inline void SetNumberOfBinsPerAxis(const unsigned char size)
  {
    this->m_LocalHaralickComputer->SetNumberOfBinsPerAxis(size);
  }

  inline const unsigned char GetNumberOfBinsPerAxis()
  {
    return this->m_LocalHaralickComputer->GetNumberOfBinsPerAxis();
  }

  void GenerateData(void)
  {
    itkDebugMacro( << "GenerateData() called");
    const InputImageType *input = this->GetInput();
    OutputImageType *output = this->GetOutput();

    this->m_LocalHaralickComputer->SetInput(input);

    typename InputImageType::RegionType windowRegion, requestedRegion = output->GetRequestedRegion();
    typename InputImageType::IndexType windowIndex;
    typename InputImageType::SizeType windowSize;

    for(unsigned int i = 0; i < itkGetStaticConstMacro(ImageDimension); ++i)
      {
      windowSize.SetElement(i, (m_WindowRadius.GetElement(i) << 1) + 1);
      }

    typedef itk::VariableLengthVector<double> VariableVectorType;
    VariableVectorType features;
    features.SetSize(8);

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

      this->m_LocalHaralickComputer->SetRegionOfInterest(windowRegion);
      this->m_LocalHaralickComputer->Update();

      features[0] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::Energy);
      features[1] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::Entropy);
      features[2] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::Correlation);
      features[3] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::InverseDifferenceMoment);
      features[4] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::Inertia);
      features[5] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::ClusterShade);
      features[6] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::ClusterProminence);
      features[7] = this->m_LocalHaralickComputer->GetFeature(LocalHaralickComputer::HaralickFeaturesComputer::HaralickCorrelation);

      outputIterator.Set(features);

      ++imageIterator;
      ++outputIterator;
    }
  }

private:
  typename LocalHaralickComputer::Pointer m_LocalHaralickComputer;
  RadiusType m_WindowRadius;
  
};

} // End of namespace Statistics
} // End of namespace itk

#endif /* __itkScalarImageToHaralickTextureFeaturesImageFilter_h */

