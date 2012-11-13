#ifndef ITKHARALICKFEATURESGLCMCALCULATOR_HXX
#define ITKHARALICKFEATURESGLCMCALCULATOR_HXX

#include "itkHaralickFeaturesGLCMCalculator.h"
#include "vnl/vnl_math.h"

namespace itk
{

template< typename TGLCMType, typename TFeatureType >
HaralickFeaturesGLCMCalculator< TGLCMType, TFeatureType >
::HaralickFeaturesGLCMCalculator()
{
	m_Features.SetSize(8);
}

template< typename TGLCMType, typename TFeatureType >
void
HaralickFeaturesGLCMCalculator< TGLCMType, TFeatureType >
::Compute(void)
{
  // Now get the various means and variances. This takes two passes
  // through the cooccurrence matrix.
  FeatureType pixelMean, marginalMean, marginalDevSquared, pixelDevSquared;

  this->ComputeMeansAndVariances(pixelMean, marginalMean, marginalDevSquared,
                                 pixelDevSquared);

  itkDebugMacro("PixelMean: " << pixelMean);
  itkDebugMacro("MarginalMean: " << marginalMean);
  itkDebugMacro("MarginalDevSquared: " << marginalDevSquared);
  itkDebugMacro("PixelDevSquared: " << pixelDevSquared);

  // Finally compute the texture features. Another one pass.
  FeatureType energy                  = NumericTraits< FeatureType >::Zero;
  FeatureType entropy                 = NumericTraits< FeatureType >::Zero;
  FeatureType correlation             = NumericTraits< FeatureType >::Zero;
  FeatureType inverseDifferenceMoment = NumericTraits< FeatureType >::Zero;
  FeatureType inertia                 = NumericTraits< FeatureType >::Zero;
  FeatureType clusterShade            = NumericTraits< FeatureType >::Zero;
  FeatureType clusterProminence       = NumericTraits< FeatureType >::Zero;
  FeatureType haralickCorrelation     = NumericTraits< FeatureType >::Zero;

  FeatureType log2 = vcl_log(2.0);

  FeatureType cooccurrenceMatrixTotalCount = m_CooccurrenceMatrix->GetTotalCount();

  GLCMIterator cooc_begin = m_CooccurrenceMatrix->Begin(), cooc_end = m_CooccurrenceMatrix->End();
  GLCMIterator cooc_it = cooc_begin;
  GLCMIndexType i1, i2;
  while ( cooc_it != cooc_end )
    {
    FeatureType frequency = (*cooc_it) / cooccurrenceMatrixTotalCount;
    if ( frequency != 0 )
      {
      m_CooccurrenceMatrix->GetIndexes( cooc_it - cooc_begin, &i1, &i2 );

      energy += frequency * frequency;
      entropy -= ( frequency > 0.0001 ) ? frequency *vcl_log(frequency) / log2:0;
      if(pixelDevSquared > 0)
        {
        correlation += ( ( i1 - pixelMean ) * ( i2 - pixelMean ) * frequency )
                       / pixelDevSquared;
        }
      inverseDifferenceMoment += frequency
                                 / ( 1.0 + ( i1 - i2 ) * ( i1 - i2 ) );
      inertia += ( i1 - i2 ) * ( i1 - i2 ) * frequency;
      clusterShade += vcl_pow( ( i1 - pixelMean ) + ( i2 - pixelMean ), 3 )
                      * frequency;
      clusterProminence += vcl_pow( ( i1 - pixelMean ) + ( i2 - pixelMean ), 4 )
                           * frequency;
      if(marginalDevSquared > 0)
        {
        haralickCorrelation += i1 * i2 * frequency;
        }
      }

    ++cooc_it;
    }

  if(marginalDevSquared > 0)
    {
    haralickCorrelation = ( haralickCorrelation - marginalMean * marginalMean )
                          / marginalDevSquared;
    }

  m_Features[0] = energy;
  m_Features[1] = entropy;
  m_Features[2] = correlation;
  m_Features[3] = inverseDifferenceMoment;
  m_Features[4] = inertia;
  m_Features[5] = clusterShade;
  m_Features[6] = clusterProminence;
  m_Features[7] = haralickCorrelation;

}

template< typename TGLCMType, typename TFeatureType >
void
HaralickFeaturesGLCMCalculator< TGLCMType, TFeatureType >
::ComputeMeansAndVariances(FeatureType & pixelMean, FeatureType & marginalMean,
                           FeatureType & marginalDevSquared, FeatureType & pixelDevSquared)
{
  // This function takes two passes through the cooccurrenceMatrix and two passes through
  // an array of the same length as a cooccurrenceMatrix axis. This could probably be
  // cleverly compressed to one pass, but it's not clear that that's necessary.

  // Initialize everything
  typename GLCMType::SizeType matrixSize = m_CooccurrenceMatrix->GetSize();
  FeatureType *marginalSums = new FeatureType[matrixSize];
  std::fill(marginalSums, marginalSums + matrixSize, 0.0);

  GLCMIterator cooc_begin = m_CooccurrenceMatrix->Begin(), cooc_end = m_CooccurrenceMatrix->End();
  GLCMIterator cooc_it = cooc_begin;
  GLCMIndexType i1, i2;

  /*  Now get the mean and deviaton of the marginal sums.
      Compute incremental mean and SD, a la Knuth, "The  Art of Computer
      Programming, Volume 2: Seminumerical Algorithms",  section 4.2.2.
      Compute mean and standard deviation using the recurrence relation:
      M(1) = x(1), M(k) = M(k-1) + (x(k) - M(k-1) ) / k
      S(1) = 0, S(k) = S(k-1) + (x(k) - M(k-1)) * (x(k) - M(k))
      for 2 <= k <= n, then
      sigma = vcl_sqrt(S(n) / n) (or divide by n-1 for sample SD instead of
      population SD).

      For weighted distribution, uses West1979, "Updating mean and variance 
      estimates: An improved method".
      SW(1) = w(1), SW(k) = SW(k-1) + SW(k)
      M(1) = x(1), M(k) = M(k-1) + (x(k) - M(k-1) ) * w(k) / SW(k)
      S(1) = 0, S(k) = S(k-1) + (x(k) - M(k-1)) * (x(k) - M(k)) * W(k)
      for 2 <= k <= n, then
      sigma = vcl_sqrt(S(n) / SW(n))
  */

  FeatureType SW_k = 0.0, W_k, M_k = 0.0, S_k = 0.0, delta;

  while ( cooc_it != cooc_end )
    {
    W_k = *cooc_it; // No need to devide by the total, because the algorithm is for weighted distributions
    if(W_k > 0)
      {
      m_CooccurrenceMatrix->GetIndexes( cooc_it - cooc_begin, &i1, &i2 );

      SW_k += W_k;

      delta = i1 - M_k;
      M_k += (W_k / SW_k) * delta;
      S_k += W_k * delta * (i1 - M_k);

      // Filling the marginal sums array
      marginalSums[i1] += W_k;
      }
    ++cooc_it;
    }


  pixelMean = M_k;
  pixelDevSquared = S_k / SW_k;

  
  W_k = marginalSums[0]; SW_k = W_k; M_k = 0.0; S_k = 0.0;

  for ( unsigned int X_k = 1; X_k < matrixSize; ++X_k )
    {
    W_k = marginalSums[X_k];
    if(W_k > 0)
      {
      SW_k += W_k;

      delta = X_k - M_k;
      M_k += (W_k / SW_k) * delta;
      S_k += W_k * delta * (X_k - M_k);
      }
    }

  marginalMean = M_k;
  marginalDevSquared = S_k / SW_k;

  delete[] marginalSums;
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

} // end namespace itk

#endif /* ITKHARALICKFEATURESGLCMCALCULATOR_HXX */
