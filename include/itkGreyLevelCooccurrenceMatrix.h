#ifndef ITKGREYLEVELCOOCCURRENCEMATRIX_H
#define ITKGREYLEVELCOOCCURRENCEMATRIX_H

#include <vector>

#include "itkMacro.h"
#include "itkObjectFactory.h"
#include "itkDataObject.h"
#include "itkSmartPointer.h"
#include "itkNumericTraits.h"
#include "itkConceptChecking.h"

namespace itk
{
namespace Statistics
{

/**
 * \class GreyLevelCooccurrenceMatrix
 * \brief This class is used to count coocurrence of two grey levels
 * in a region of an image.
 *
 * The class is templated over the type of measurements.
 *
 * A cooccurrence matrix can be seen as a 2-dimensional histogram where each
 * dimension represent the occurrence of grey level.
 *
 * Grey levels are represented by their value, from 0 to N-1, where N is the
 * size of the matrix.
 **/

template< class TMeasurementType = unsigned int >
class GreyLevelCooccurrenceMatrix:private std::vector<TMeasurementType>, public DataObject
{
public:

  /** Standard class typedefs */
  typedef GreyLevelCooccurrenceMatrix Self;
  typedef DataObject                  Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(GreyLevelCooccurrenceMatrix, DataObject);

  typedef TMeasurementType                          MeasurementType;
  typedef std::vector<MeasurementType>              CounterContainer;
  typedef unsigned int                              SizeType;
  typedef SizeType                                  IndexType;
  typedef typename CounterContainer::iterator       Iterator;
  typedef typename CounterContainer::const_iterator ConstIterator;

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( MeasurementHasNumericTraitsCheck,
    ( Concept::HasNumericTraits< MeasurementType > ) );
#endif

  GreyLevelCooccurrenceMatrix()
    :CounterContainer(), m_Size(itk::NumericTraits<SizeType>::Zero), m_TotalCount(itk::NumericTraits<MeasurementType>::Zero)
  {
    itkDebugMacro(<< "Constructing an empty GreyLevelCooccurrenceMatrix.");
  }

  itkNewMacro(Self);

  /** Resizes the matrix to the specified size (only if necessary)
   * and set every measures to 0. */
  void SetSize(const SizeType size);

  inline unsigned int GetSize(void) const
  {
    return m_Size;
  }

  inline MeasurementType GetTotalCount(void) const
  {
    return m_TotalCount;
  }
  inline void SetTotalCount(MeasurementType v)
  {
    m_TotalCount = v;
  }

  inline void SetToZero()
  {
    std::fill( CounterContainer::begin(), CounterContainer::end(), NumericTraits< MeasurementType >::Zero );
    m_TotalCount = 0;
  }

  inline void IncrementCounter(const unsigned int v1, const unsigned int v2)
  {
    CounterContainer::operator[](ComputeOffset(v1, v2)) += 1;
    ++m_TotalCount;
  }

  inline ConstIterator Begin(void) const
  {
    return CounterContainer::begin();
  }

  inline ConstIterator End(void) const
  {
    return CounterContainer::end();
  }

  inline Iterator Begin(void)
  {
    return CounterContainer::begin();
  }

  inline Iterator End(void)
  {
    return CounterContainer::end();
  }

  inline void GetIndexes(const IndexType index, IndexType *v1, IndexType *v2) const
  {
    *v1 = index % m_Size;
    *v2 = index / m_Size;
  }

  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  IndexType m_Size;
  MeasurementType m_TotalCount;

  inline IndexType ComputeOffset(const IndexType v1, const IndexType v2) const
  {
    return v2 * m_Size + v1;
  }

};
} // End namespace Statistics
} // End namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGreyLevelCooccurrenceMatrix.hxx"
#endif

#endif /* ITKGREYLEVELCOOCCURRENCEMATRIX_H */

