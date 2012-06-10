#ifndef __itkScalarImageToHaralickTextureFeaturesFilter_h
#define __itkScalarImageToHaralickTextureFeaturesFilter_h

#include "itkImageToImageFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkCooccurrenceMatrixNormalizerFilter.h"
#include "itkCooccurrenceMatrixToHaralickTextureFeaturesFilter.h"

namespace itk
{
namespace Statistics
{

template< class TInputImageType, class TOutputPixelType >
class ITK_EXPORT ScalarImageToHaralickTextureFeaturesFilter:public ImageToImageFilter< TInputImageType, TInputImageType>
{
public:
  /** Standard typedefs */
  typedef ScalarImageToHaralickTextureFeaturesFilter              Self;
  typedef ImageToImageFilter< TInputImageType, TInputImageType>   Superclass;
  typedef SmartPointer< Self >                                    Pointer;
  typedef SmartPointer< const Self >                              ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScalarImageToHaralickTextureFeaturesFilter, ImageToImageFilter);

  /** standard New() method support */
  itkNewMacro(Self);

  typedef TInputImageType                                   InputImageType;
  typedef typename InputImageType::Pointer                  InputImagePointer;
  typedef typename InputImageType::ConstPointer             InputImageConstPointer;
  typedef typename InputImageType::PixelType                PixelType;
  typedef typename InputImageType::RegionType               RegionType;
  typedef typename InputImageType::SizeType                 RadiusType;
  typedef typename InputImageType::OffsetType               OffsetType;

  typedef TOutputPixelType                                  OutputPixelType;

  typedef ScalarImageToCooccurrenceMatrixFilter< InputImageType >                                CooccurrenceMatrixComputer;
  typedef typename CooccurrenceMatrixComputer::CooccurrenceMatrixType                            CooccurrenceMatrixType;
  typedef typename CooccurrenceMatrixComputer::OffsetVector                                      OffsetVectorType;

  typedef CooccurrenceMatrixNormalizerFilter< CooccurrenceMatrixType, OutputPixelType >          CooccurrenceMatrixNormalizer;
  typedef typename CooccurrenceMatrixNormalizer::NormalizedCooccurrenceMatrixType                NormalizedCooccurrenceMatrixType;
  typedef CooccurrenceMatrixToHaralickTextureFeaturesFilter< NormalizedCooccurrenceMatrixType >  HaralickFeaturesComputer;

  ScalarImageToHaralickTextureFeaturesFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);

    this->ProcessObject::SetNthOutput( 0, static_cast< DataObject * >( InputImageType::New().GetPointer() ) );

    this->m_CooccurrenceMatrixComputer = CooccurrenceMatrixComputer::New();
    this->m_CooccurrenceMatrixNormalizer = CooccurrenceMatrixNormalizer::New();
    this->m_CooccurrenceMatrixNormalizer->SetInput(this->m_CooccurrenceMatrixComputer->GetOutput());
    this->m_HaralickFeaturesComputer = HaralickFeaturesComputer::New();
    this->m_HaralickFeaturesComputer->SetInput(this->m_CooccurrenceMatrixNormalizer->GetOutput());
  }

  using Superclass::SetInput;
  void SetInput(const InputImageType *image)
  {
    this->ProcessObject::SetNthInput( 0, const_cast< InputImageType * >( image ) );

    this->m_CooccurrenceMatrixComputer->SetInput(image);

    this->GraftOutput(const_cast< InputImageType * >( image ));
  }


  inline typename NormalizedCooccurrenceMatrixType::MeasurementType GetFeature(typename HaralickFeaturesComputer::TextureFeatureName feature) const
  {
    return this->m_HaralickFeaturesComputer->GetFeature(feature);
  }


  inline void SetOffset(const OffsetType offset)
  {
    this->m_CooccurrenMatrixComputer->SetOffset(offset);
  }

  inline void SetOffsets(const OffsetVectorType * offsets)
  {
    this->m_CooccurrenceMatrixComputer->SetOffsets(offsets);
    this->Modified();
  }

  inline const OffsetVectorType GetOffsets() const
  {
    return this->m_CooccurrenceMatrixComputer->GetOffset();
  }


  inline void SetRegionOfInterest(const RegionType roi)
  {
    this->m_CooccurrenceMatrixComputer->SetRegionOfInterest(roi);
    this->Modified();
  }

  inline const RegionType GetRegionOfInterest() const
  {
    return this->m_CooccurrenceMatrixComputer->GetRegionOfInterest();
  }


  inline void SetNumberOfBinsPerAxis(const unsigned char size)
  {
    this->m_CooccurrenceMatrixComputer->SetNumberOfBinsPerAxis(size);
    this->Modified();
  }

  inline const unsigned char GetNumberOfBinsPerAxis()
  {
    return this->m_CooccurrenceMatrixComputer->GetNumberOfBinsPerAxis();
  }

  
  void GenerateData(void)
  {
    itkDebugMacro(<< "GenerateData() called");
    this->m_HaralickFeaturesComputer->Update();
  }

private:
  typename CooccurrenceMatrixComputer::Pointer     m_CooccurrenceMatrixComputer;
  typename CooccurrenceMatrixNormalizer::Pointer   m_CooccurrenceMatrixNormalizer;
  typename HaralickFeaturesComputer::Pointer       m_HaralickFeaturesComputer;

}; 

} // End of namespace Statistics
} // End of namespace itk

#endif /* __itkScalarImageToHaralickTextureFeaturesFilter_h */
