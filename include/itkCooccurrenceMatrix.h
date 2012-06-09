#ifndef __itkCooccurrenceMatrix_h
#define __itkCooccurrenceMatrix_h

#include <vector>

#include "itkMacro.h"
#include "itkObjectFactory.h"
#include "itkDataObject.h"
#include "itkSmartPointer.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Statistics
{

template< class TMeasurementType = unsigned int > 
class CooccurrenceMatrix:private std::vector<TMeasurementType>, public DataObject
{
public:

  /** Standard class typedefs */
  typedef CooccurrenceMatrix          Self;
  typedef DataObject                 Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(CooccurrenceMatrix, DataObject);

  typedef TMeasurementType MeasurementType;
  typedef std::vector<MeasurementType> CounterContainer;
  typedef unsigned int SizeType;
  typedef SizeType IndexType;
  typedef typename CounterContainer::iterator Iterator;
  typedef typename CounterContainer::const_iterator ConstIterator;

  CooccurrenceMatrix()
    :CounterContainer(), m_Size(0)
  {}

  itkNewMacro(Self);

  void SetSize(const SizeType size)
  {
    if(m_Size != size)
      {
      m_Size = size;
      CounterContainer::resize(m_Size * m_Size);
      }
  }

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
#include "itkCooccurrenceMatrix.hxx"
#endif

#endif /* __itkCooccurrenceMatrix_h */

