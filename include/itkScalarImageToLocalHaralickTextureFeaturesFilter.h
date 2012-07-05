#ifndef __itkScalarImageToLocalHaralickTextureFeaturesFilter_h
#define __itkScalarImageToLocalHaralickTextureFeaturesFilter_h

#include "itkImageToImageFilter.h"
#include "itkScalarImageToGreyLevelCooccurrenceMatrixFilter.h"
#include "itkGreyLevelCooccurrenceMatrixNormalizerFilter.h"
#include "itkGreyLevelCooccurrenceMatrixToHaralickTextureFeaturesFilter.h"

namespace itk
{
namespace Statistics
{

template< class TInputImageType, class TOutputPixelType >
class ITK_EXPORT ScalarImageToLocalHaralickTextureFeaturesFilter:public ImageToImageFilter< TInputImageType, TInputImageType>
{
public:
  /** Standard typedefs */
  typedef ScalarImageToLocalHaralickTextureFeaturesFilter              Self;
  typedef ImageToImageFilter< TInputImageType, TInputImageType>   Superclass;
  typedef SmartPointer< Self >                                    Pointer;
  typedef SmartPointer< const Self >                              ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScalarImageToLocalHaralickTextureFeaturesFilter, ImageToImageFilter);

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

  typedef ScalarImageToGreyLevelCooccurrenceMatrixFilter< InputImageType >                                GreyLevelCooccurrenceMatrixComputer;
  typedef typename GreyLevelCooccurrenceMatrixComputer::GreyLevelCooccurrenceMatrixType                            GreyLevelCooccurrenceMatrixType;
  typedef typename GreyLevelCooccurrenceMatrixComputer::OffsetVector                                      OffsetVectorType;

  typedef GreyLevelCooccurrenceMatrixNormalizerFilter< GreyLevelCooccurrenceMatrixType, OutputPixelType >          GreyLevelCooccurrenceMatrixNormalizer;
  typedef typename GreyLevelCooccurrenceMatrixNormalizer::NormalizedGreyLevelCooccurrenceMatrixType                NormalizedGreyLevelCooccurrenceMatrixType;
  typedef GreyLevelCooccurrenceMatrixToHaralickTextureFeaturesFilter< NormalizedGreyLevelCooccurrenceMatrixType >  HaralickFeaturesComputer;

  ScalarImageToLocalHaralickTextureFeaturesFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);

    this->ProcessObject::SetNthOutput( 0, static_cast< DataObject * >( InputImageType::New().GetPointer() ) );

    this->m_GreyLevelCooccurrenceMatrixComputer = GreyLevelCooccurrenceMatrixComputer::New();
    this->m_GreyLevelCooccurrenceMatrixNormalizer = GreyLevelCooccurrenceMatrixNormalizer::New();
    this->m_GreyLevelCooccurrenceMatrixNormalizer->SetInput(this->m_GreyLevelCooccurrenceMatrixComputer->GetOutput());
    this->m_HaralickFeaturesComputer = HaralickFeaturesComputer::New();
    this->m_HaralickFeaturesComputer->SetInput(this->m_GreyLevelCooccurrenceMatrixNormalizer->GetOutput());
  }

  using Superclass::SetInput;
  void SetInput(const InputImageType *image)
  {
    this->ProcessObject::SetNthInput( 0, const_cast< InputImageType * >( image ) );

    this->m_GreyLevelCooccurrenceMatrixComputer->SetInput(image);

    this->GraftOutput(const_cast< InputImageType * >( image ));
  }


  inline typename NormalizedGreyLevelCooccurrenceMatrixType::MeasurementType GetFeature(typename HaralickFeaturesComputer::TextureFeatureName feature) const
  {
    return this->m_HaralickFeaturesComputer->GetFeature(feature);
  }


  inline void SetOffset(const OffsetType offset)
  {
    this->m_CooccurrenMatrixComputer->SetOffset(offset);
    this->Modified();
  }

  inline void SetOffsets(const OffsetVectorType * offsets)
  {
    this->m_GreyLevelCooccurrenceMatrixComputer->SetOffsets(offsets);
    this->Modified();
  }

  inline const OffsetVectorType GetOffsets() const
  {
    return this->m_GreyLevelCooccurrenceMatrixComputer->GetOffset();
  }


  inline void SetRegionOfInterest(const RegionType roi)
  {
    this->m_GreyLevelCooccurrenceMatrixComputer->SetRegionOfInterest(roi);
    this->Modified();
  }

  inline const RegionType GetRegionOfInterest() const
  {
    return this->m_GreyLevelCooccurrenceMatrixComputer->GetRegionOfInterest();
  }


  inline void SetNumberOfBinsPerAxis(const unsigned char size)
  {
    this->m_GreyLevelCooccurrenceMatrixComputer->SetNumberOfBinsPerAxis(size);
    this->Modified();
  }

  inline const unsigned char GetNumberOfBinsPerAxis()
  {
    return this->m_GreyLevelCooccurrenceMatrixComputer->GetNumberOfBinsPerAxis();
  }

  
  void GenerateData(void)
  {
    itkDebugMacro(<< "GenerateData() called");
    this->m_HaralickFeaturesComputer->Update();
  }

private:
  typename GreyLevelCooccurrenceMatrixComputer::Pointer     m_GreyLevelCooccurrenceMatrixComputer;
  typename GreyLevelCooccurrenceMatrixNormalizer::Pointer   m_GreyLevelCooccurrenceMatrixNormalizer;
  typename HaralickFeaturesComputer::Pointer       m_HaralickFeaturesComputer;

}; 

} // End of namespace Statistics
} // End of namespace itk

#endif /* __itkScalarImageToLocalHaralickTextureFeaturesFilter_h */
