#ifndef __itkCoocurrenceMatrix_h
#define __itkCoocurrenceMatrix_h

#include <vector>

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkMacro.h>
#include <itkSmartPointer.h>
#include <itkNumericTraits.h>

namespace itk
{
namespace Statistics
{

template< class TValueType = unsigned int, class TFrequencyType = float > 
class CoocurrenceMatrix:private std::vector<TFrequencyType>
{
public:

  typedef TValueType ValueType;
  typedef TFrequencyType FrequencyType;
  typedef std::vector<FrequencyType> FrequencyContainer;
  typedef typename FrequencyContainer::const_iterator ConstIterator;

  CoocurrenceMatrix(const unsigned int size = 0)
    :FrequencyContainer()
  {
    SetSize(size);
  }

  inline void SetSize(const unsigned int size)
  {
    m_Size = size;
    if(size > 0)
    {
      FrequencyContainer::resize(m_Size);
      Reset();
    }
  }

  inline unsigned int GetSize(void)
  {
    return m_Size;
  }

  void Initialize(const unsigned int size);

  inline void Reset()
  {
    std::fill( FrequencyContainer::begin(), FrequencyContainer::end(), NumericTraits< FrequencyType >::Zero );
    m_TotalFrequency = 0;
  }

  inline void IncrementFrequency(const ValueType v1, const ValueType v2)
  {
    FrequencyContainer::operator[](ComputeOffset(v1, v2)) += 1;
  }

  using FrequencyContainer::begin;
  using FrequencyContainer::end;

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

