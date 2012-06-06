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
  typedef typename CounterContainer::const_iterator ConstIterator;

  CoocurrenceMatrix()
    :CounterContainer()
  {}

  itkNewMacro(Self);


  void SetSize(const unsigned int size)
  {
    m_Size = size;
    CounterContainer::resize(m_Size * m_Size);
    SetToZero();
  }

  inline unsigned int GetSize(void) const
  {
    return m_Size;
  }

  inline unsigned int GetTotalCount(void) const
  {
    return m_TotalCount;
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

  /*
  void Normalize() {
    typename CounterContainer::iterator it = CounterContainer::begin(),
             end = CounterContainer::end();
    while(it < end)
    {
      *it /= m_TotalCount;
      ++it;
    }
    m_TotalCount = 1;
  }
  */

  inline ConstIterator Begin(void) const
  {
    return CounterContainer::begin();
  }

  inline ConstIterator End(void) const
  {
    return CounterContainer::end();
  }

private:
  unsigned int m_Size;
  unsigned int m_TotalCount;

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

