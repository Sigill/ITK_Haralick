#ifndef ITKHARALICKFEATURESGLCMCALCULATOR_H
#define ITKHARALICKFEATURESGLCMCALCULATOR_H

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkVariableLengthVector.h"

namespace itk
{
template< typename TGLCMType, typename TFeatureType >
class ITK_EXPORT HaralickFeaturesGLCMCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef HaralickFeaturesGLCMCalculator               Self;
  typedef Object                                       Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  typedef TGLCMType                                    GLCMType;
  typedef typename GLCMType::Pointer                   GLCMPointer;
  typedef typename GLCMType::ConstPointer              GLCMConstPointer;
  typedef typename GLCMType::SizeType                  GLCMSizeType;
  typedef typename GLCMType::IndexType                 GLCMIndexType;
  typedef typename GLCMType::ConstIterator             GLCMIterator;

  typedef TFeatureType                                 FeatureType;

  typedef itk::VariableLengthVector< FeatureType >     FeaturesVectorType;

  itkNewMacro(Self)

  itkTypeMacro(HaralickFeaturesGLCMCalculator, Object)

  itkSetConstObjectMacro(CooccurrenceMatrix, GLCMType)

  itkGetConstReferenceMacro(Features, FeaturesVectorType)

  void Compute(void);

  /** Texture feature types */
  typedef enum {
    Energy,
    Entropy,
    Correlation,
    InverseDifferenceMoment,
    Inertia,
    ClusterShade,
    ClusterProminence,
    HaralickCorrelation,
    InvalidFeatureName
    } TextureFeatureName;

protected:
  HaralickFeaturesGLCMCalculator();
  virtual ~HaralickFeaturesGLCMCalculator() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  void ComputeMeansAndVariances(FeatureType & pixelMean, FeatureType & marginalMean,
                                FeatureType & marginalDevSquared, FeatureType & pixelVariance);

private:
  HaralickFeaturesGLCMCalculator(const Self &); //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  GLCMConstPointer m_CooccurrenceMatrix;

  FeaturesVectorType m_Features;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHaralickFeaturesGLCMCalculator.hxx"
#endif

#endif /* ITKHARALICKFEATURESGLCMCALCULATOR_H */
