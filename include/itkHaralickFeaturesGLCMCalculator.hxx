#ifndef ITKHARALICKFEATURESGLCMCALCULATOR_HXX
#define ITKHARALICKFEATURESGLCMCALCULATOR_HXX

#include "itkHaralickFeaturesGLCMCalculator.h"
#include "vnl/vnl_math.h"
#include "IncrementalWeightedMeanVarianceComputer.h"

namespace itk
{
namespace Statistics
{

template< typename TGLCMType, typename TFeatureType >
HaralickFeaturesGLCMCalculator< TGLCMType, TFeatureType >
::HaralickFeaturesGLCMCalculator()
{
  m_Features.SetSize(NumberOfFeatures);
}

template< typename TGLCMType, typename TFeatureType >
void
HaralickFeaturesGLCMCalculator< TGLCMType, TFeatureType >
::Compute(void)
{
  const typename GLCMType::SizeType matrixSize = m_CooccurrenceMatrix->GetSize();
  const typename GLCMType::SizeType maxGreylevel = matrixSize - 1;

  std::vector< FeatureType > sum_occurrences( 2 * matrixSize - 1, NumericTraits< FeatureType >::Zero),
                             diff_occurrences(matrixSize,         NumericTraits< FeatureType >::Zero);

  GLCMIterator cooc_begin = m_CooccurrenceMatrix->Begin(), cooc_end = m_CooccurrenceMatrix->End();
  GLCMIterator cooc_it = cooc_begin;
  GLCMIndexType i1, i2;

  IncrementalWeightedMeanVarianceComputer< FeatureType > 
    marginal_x_stats_computer, marginal_y_stats_computer,
    sum_stats_computer, diff_stats_computer;

  FeatureType angularSecondMoment     = NumericTraits< FeatureType >::Zero;
  FeatureType contrast                = NumericTraits< FeatureType >::Zero;
  FeatureType variance                = NumericTraits< FeatureType >::Zero;
  FeatureType inverseDifferenceMoment = NumericTraits< FeatureType >::Zero;
  FeatureType sumAverage              = NumericTraits< FeatureType >::Zero;
  FeatureType sumVariance             = NumericTraits< FeatureType >::Zero;
  FeatureType sumEntropy              = NumericTraits< FeatureType >::Zero;
  FeatureType entropy                 = NumericTraits< FeatureType >::Zero;
  FeatureType differenceVariance      = NumericTraits< FeatureType >::Zero;
  FeatureType differenceEntropy       = NumericTraits< FeatureType >::Zero;

  // How to handle a texture with a single color (no correlation) ?
  //FeatureType correlation             = NumericTraits< FeatureType >::Zero;
  // What is the range ?
  //FeatureType clusterShade            = NumericTraits< FeatureType >::Zero;
  // What is the range ?
  //FeatureType clusterProminence       = NumericTraits< FeatureType >::Zero;
  // This is the same thing as correlation, with the same problem.
  //FeatureType haralickCorrelation     = NumericTraits< FeatureType >::Zero;

  FeatureType cooccurrenceMatrixTotalCount = m_CooccurrenceMatrix->GetTotalCount();
  FeatureType frequency;
  const FeatureType log2 = vcl_log(2.0);

  while ( cooc_it != cooc_end )
    {
    frequency = (*cooc_it) / cooccurrenceMatrixTotalCount;

    if(frequency > 0)
      {
      m_CooccurrenceMatrix->GetIndexes( cooc_it - cooc_begin, &i1, &i2 );

      marginal_x_stats_computer.insert(i1, frequency);
      marginal_y_stats_computer.insert(i2, frequency);
      sum_occurrences[i1+i2] += frequency;
      diff_occurrences[abs(i1-i2)] += frequency;

      angularSecondMoment += frequency * frequency;
      entropy -= ( frequency > NumericTraits< FeatureType >::Zero ) ? frequency *vcl_log(frequency) / log2 : 0;
      inverseDifferenceMoment += frequency / (1 + (i1 - i2) * (i1 - i2));
      }
    ++cooc_it;
    }

  for(int k = 0; k < matrixSize; ++k)
    {
    frequency = diff_occurrences[k];

    if(frequency > 0)
      {
      diff_stats_computer.insert(k, frequency);

      contrast += k * k * frequency;
      differenceEntropy -= ( frequency > NumericTraits< FeatureType >::Zero ) ? frequency *vcl_log(frequency) / log2 : 0;
      }
    }

  differenceVariance = diff_stats_computer.getVariance();

  for(int k = 0; k < 2 * matrixSize - 1; ++k)
    {
    frequency = sum_occurrences[k];

    if(frequency > 0)
      {
      sum_stats_computer.insert(k, frequency);

      sumEntropy -= ( frequency > NumericTraits< FeatureType >::Zero ) ? frequency *vcl_log(frequency) / log2 : 0;
	  }
    }

  variance = marginal_x_stats_computer.getVariance();
  sumAverage = sum_stats_computer.getMean();
  sumVariance = sum_stats_computer.getVariance();


  m_Features[0] = angularSecondMoment;
  m_Features[1] = contrast / (maxGreylevel  * maxGreylevel);
  m_Features[2] = variance * 4.0 / (maxGreylevel  * maxGreylevel);
  m_Features[3] = inverseDifferenceMoment;
  m_Features[4] = sumAverage / (2 * maxGreylevel);
  m_Features[5] = sumVariance / (maxGreylevel * maxGreylevel);
  m_Features[6] = sumEntropy / (vcl_log(2 * matrixSize - 1) / log2);
  m_Features[7] = entropy / (2*vcl_log(matrixSize) / log2);
  m_Features[8] = differenceVariance * 4.0 / (maxGreylevel  * maxGreylevel);
  m_Features[9] = differenceEntropy / (vcl_log(matrixSize) / log2);
}

template< typename TGLCMType, typename TFeatureType >
void
HaralickFeaturesGLCMCalculator< TGLCMType, TFeatureType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  /*
  os << indent << "Image: " << std::endl;
  m_Image->Print( os, indent.GetNextIndent() );
  os << indent << "Region: " << std::endl;
  m_Region.Print( os, indent.GetNextIndent() );
  os << indent << "Region set by User: " << m_RegionSetByUser << std::endl;
  */
}

} // end namespace Statistics
} // end namespace itk

#endif /* ITKHARALICKFEATURESGLCMCALCULATOR_HXX */
