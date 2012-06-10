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
#ifndef __itkCooccurrenceMatrixToHaralickTextureFeaturesFilter_hxx
#define __itkCooccurrenceMatrixToHaralickTextureFeaturesFilter_hxx

#include "itkCooccurrenceMatrixToHaralickTextureFeaturesFilter.h"

#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"

namespace itk
{
namespace Statistics
{
//constructor
template< class TCooccurrenceMatrix >
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::CooccurrenceMatrixToHaralickTextureFeaturesFilter(void)
{
  this->ProcessObject::SetNumberOfRequiredInputs(1);

  // allocate the data objects for the outputs which are
  // just decorators real types
  for ( int i = 0; i < 8; ++i )
    {
    this->ProcessObject::SetNthOutput( i, this->MakeOutput(i) );
    }
}

template< class TCooccurrenceMatrix >
void
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::SetInput(const CooccurrenceMatrixType *cooccurrenceMatrix)
{
  this->ProcessObject::SetNthInput( 0, const_cast< CooccurrenceMatrixType * >( cooccurrenceMatrix ) );
}

template< class TCooccurrenceMatrix >
const typename
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::CooccurrenceMatrixType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetInput() const
{
  return static_cast< const CooccurrenceMatrixType * >( this->GetPrimaryInput() );
}

template< class TCooccurrenceMatrix >
typename
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::DataObjectPointer
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
{
  return static_cast< DataObject * >( MeasurementObjectType::New().GetPointer() );
}

template< class TCooccurrenceMatrix >
void
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::GenerateData(void)
{
  itkDebugMacro(<< "GenerateData() called");
  typedef typename CooccurrenceMatrixType::ConstIterator CooccurrenceMatrixIterator;

  const CooccurrenceMatrixType *inputCooccurrenceMatrix = this->GetInput();

  // Now get the various means and variances. This is takes two passes
  // through the cooccurrenceMatrix.
  double pixelMean;
  double marginalMean;
  double marginalDevSquared;
  double pixelDevSquared;

  this->ComputeMeansAndVariances(pixelMean, marginalMean, marginalDevSquared,
                                 pixelDevSquared);

  // Finally compute the texture features. Another one pass.
  MeasurementType energy      = NumericTraits< MeasurementType >::Zero;
  MeasurementType entropy     = NumericTraits< MeasurementType >::Zero;
  MeasurementType correlation = NumericTraits< MeasurementType >::Zero;

  MeasurementType inverseDifferenceMoment      =
    NumericTraits< MeasurementType >::Zero;

  MeasurementType inertia             = NumericTraits< MeasurementType >::Zero;
  MeasurementType clusterShade        = NumericTraits< MeasurementType >::Zero;
  MeasurementType clusterProminence   = NumericTraits< MeasurementType >::Zero;
  MeasurementType haralickCorrelation = NumericTraits< MeasurementType >::Zero;

  double log2 = vcl_log(2.0);

  CooccurrenceMatrixIterator cooc_begin = inputCooccurrenceMatrix->Begin(), cooc_end = inputCooccurrenceMatrix->End();
  CooccurrenceMatrixIterator cooc_it = cooc_begin;
  IndexType i1, i2;
  while ( cooc_it != cooc_end )
    {
    MeasurementType frequency = *cooc_it;
    if ( frequency != 0 )
      {
      inputCooccurrenceMatrix->GetIndexes( cooc_it - cooc_begin, &i1, &i2 );

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

  MeasurementObjectType *energyOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(0) );
  energyOutputObject->Set(energy);

  MeasurementObjectType *entropyOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(1) );
  entropyOutputObject->Set(entropy);

  MeasurementObjectType *correlationOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(2) );
  correlationOutputObject->Set(correlation);

  MeasurementObjectType *inverseDifferenceMomentOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(3) );
  inverseDifferenceMomentOutputObject->Set(inverseDifferenceMoment);

  MeasurementObjectType *inertiaOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(4) );
  inertiaOutputObject->Set(inertia);

  MeasurementObjectType *clusterShadeOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(5) );
  clusterShadeOutputObject->Set(clusterShade);

  MeasurementObjectType *clusterProminenceOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(6) );
  clusterProminenceOutputObject->Set(clusterProminence);

  MeasurementObjectType *haralickCorrelationOutputObject =
    static_cast< MeasurementObjectType * >( this->ProcessObject::GetOutput(7) );
  haralickCorrelationOutputObject->Set(haralickCorrelation);
}

template< class TCooccurrenceMatrix >
void
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::ComputeMeansAndVariances(double & pixelMean,
                                                                                                 double & marginalMean,
                                                                                                 double & marginalDevSquared,
                                                                                                 double & pixelDevSquared)
{
  // This function takes two passes through the cooccurrenceMatrix and two passes through
  // an array of the same length as a cooccurrenceMatrix axis. This could probably be
  // cleverly compressed to one pass, but it's not clear that that's necessary.

  typedef typename CooccurrenceMatrixType::ConstIterator CooccurrenceMatrixIterator;

  const CooccurrenceMatrixType *inputCooccurrenceMatrix =  this->GetInput();

  // Initialize everything
  typename CooccurrenceMatrixType::SizeType matrixSize = inputCooccurrenceMatrix->GetSize();
  double *marginalSums = new double[matrixSize];
  std::fill(marginalSums, marginalSums + matrixSize, 0.0);

  CooccurrenceMatrixIterator cooc_begin = inputCooccurrenceMatrix->Begin(), cooc_end = inputCooccurrenceMatrix->End();
  CooccurrenceMatrixIterator cooc_it = cooc_begin;
  IndexType i1, i2;

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

  double SW_k = 0.0, W_k, M_k = 0.0, S_k = 0.0, delta;

  while ( cooc_it != cooc_end )
    {
    W_k = *cooc_it;
    if(W_k > 0)
      {
      inputCooccurrenceMatrix->GetIndexes( cooc_it - cooc_begin, &i1, &i2 );

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

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetEnergyOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(0) );
}

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetEntropyOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(1) );
}

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetCorrelationOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(2) );
}

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetInverseDifferenceMomentOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(3) );
}

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetInertiaOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(4) );
}

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetClusterShadeOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(5) );
}

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetClusterProminenceOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(6) );
}

template< class TCooccurrenceMatrix >
const
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementObjectType *
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetHaralickCorrelationOutput() const
{
  return static_cast< const MeasurementObjectType * >( this->ProcessObject::GetOutput(7) );
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetEnergy() const
{
  return this->GetEnergyOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetEntropy() const
{
  return this->GetEntropyOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetCorrelation() const
{
  return this->GetCorrelationOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetInverseDifferenceMoment() const
{
  return this->GetInverseDifferenceMomentOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetInertia() const
{
  return this->GetInertiaOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetClusterShade() const
{
  return this->GetClusterShadeOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetClusterProminence() const
{
  return this->GetClusterProminenceOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetHaralickCorrelation() const
{
  return this->GetHaralickCorrelationOutput()->Get();
}

template< class TCooccurrenceMatrix >
typename CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >::MeasurementType
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::GetFeature(TextureFeatureName feature)
{
  switch ( feature )
    {
    case Energy :
      return this->GetEnergy();
    case Entropy:
      return this->GetEntropy();
    case Correlation:
      return this->GetCorrelation();
    case InverseDifferenceMoment:
      return this->GetInverseDifferenceMoment();
    case Inertia:
      return this->GetInertia();
    case ClusterShade:
      return this->GetClusterShade();
    case ClusterProminence:
      return this->GetClusterProminence();
    case HaralickCorrelation:
      return this->GetHaralickCorrelation();
    default:
      return 0;
    }
}

template< class TCooccurrenceMatrix >
void
CooccurrenceMatrixToHaralickTextureFeaturesFilter< TCooccurrenceMatrix >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end of namespace Statistics
} // end of namespace itk

#endif
