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

template< class TValueType = unsigned int, class TFrequencyType = float > 
class CoocurrenceMatrix:private std::vector<TFrequencyType>, public DataObject
{
public:

  /** Standard class typedefs */
  typedef CoocurrenceMatrix          Self;
  typedef DataObject                 Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(CoocurrenceMatrix, DataObject);

  typedef TValueType ValueType;
  typedef TFrequencyType FrequencyType;
  typedef std::vector<FrequencyType> FrequencyContainer;
  typedef typename FrequencyContainer::const_iterator ConstIterator;

  CoocurrenceMatrix()
    :FrequencyContainer()
  {}

  itkNewMacro(Self);


  void SetSize(const unsigned int size)
  {
    m_Size = size;
    FrequencyContainer::resize(m_Size * m_Size);
    Reset();
  }

  inline unsigned int GetSize(void) const
  {
    return m_Size;
  }

  inline void Reset()
  {
    std::fill( FrequencyContainer::begin(), FrequencyContainer::end(), NumericTraits< FrequencyType >::Zero );
    m_TotalFrequency = 0;
  }

  inline void IncrementFrequency(const ValueType v1, const ValueType v2)
  {
    FrequencyContainer::operator[](ComputeOffset(v1, v2)) += 1;
    ++m_TotalFrequency;
  }

  void Normalize() {
    typename FrequencyContainer::iterator it = FrequencyContainer::begin(),
             end = FrequencyContainer::end();
    while(it < end)
    {
      *it /= m_TotalFrequency;
      ++it;
    }
    m_TotalFrequency = 1;
  }

  inline ConstIterator Begin(void) const
  {
    return FrequencyContainer::begin();
  }

  inline ConstIterator End(void) const
  {
    return FrequencyContainer::end();
  }

private:
  unsigned int m_Size;
  unsigned int m_TotalFrequency;

  inline unsigned int ComputeOffset(const ValueType v1, const ValueType v2)
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

