#ifndef __itkCoocurrenceMatrix_h
#define __itkCoocurrenceMatrix_h

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

template< class TCounterType = unsigned int > 
class CoocurrenceMatrix:private std::vector<TCounterType>, public DataObject
{
public:

  /** Standard class typedefs */
  typedef CoocurrenceMatrix          Self;
  typedef DataObject                 Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(CoocurrenceMatrix, DataObject);

  typedef TCounterType CounterType;
  typedef std::vector<CounterType> CounterContainer;
  typedef typename CounterContainer::iterator Iterator;
  typedef typename CounterContainer::const_iterator ConstIterator;

  CoocurrenceMatrix()
    :CounterContainer(), m_Size(0)
  {}

  itkNewMacro(Self);

  void SetSize(const unsigned int size)
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

  inline CounterType GetTotalCount(void) const
  {
    return m_TotalCount;
  }
  inline void SetTotalCount(CounterType v)
  {
    m_TotalCount = v;
  }

  inline void SetToZero()
  {
    std::fill( CounterContainer::begin(), CounterContainer::end(), NumericTraits< CounterType >::Zero );
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

private:
  unsigned int m_Size;
  CounterType m_TotalCount;

  inline unsigned int ComputeOffset(const unsigned int v1, const unsigned int v2) const
  {
    return v2 * m_Size + v1;
  }

};
} // End namespace Statistics
} // End namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCoocurrenceMatrix.hxx"
#endif

#endif /* __itkCoocurrenceMatrix_h */

