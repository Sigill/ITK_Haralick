#ifndef ITKGREYLEVELCOOCCURRENCEMATRIX_HXX
#define ITKGREYLEVELCOOCCURRENCEMATRIX_HXX

#include <iomanip>

#include "itkGreyLevelCooccurrenceMatrix.h"

namespace itk
{
namespace Statistics
{

template< typename TMeasurementType >
void
GreyLevelCooccurrenceMatrix< TMeasurementType >
::SetSize(const SizeType size)
{
  itkDebugMacro(<< "Resizing the GreyLevelCooccurrenceMatrix to " << size << ".");
  if(m_Size != size)
    {
    m_Size = size;
    CounterContainer::resize(m_Size * m_Size);
    }
  this->SetToZero();
}

template< typename TMeasurementType >
void
GreyLevelCooccurrenceMatrix< TMeasurementType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Size: " << this->GetSize() << std::endl;
  os << indent << "TotalCount: " << this->GetTotalCount() << std::endl;
  os << indent << "Values: " << std::endl;

  ConstIterator it = this->Begin(), begin = this->Begin(), end = this->End();

  IndexType i1, i2;
  while(it != end)
  {
    GetIndexes(it - begin, &i1, &i2);
    os << indent.GetNextIndent() << "(" << std::setw(4) << setiosflags(std::ios::left) << i1 << "; " 
      << std::setw(4) << setiosflags(std::ios::left) << i2 << "): " << (*it) << std::endl;
    ++it;
  }
}

} // End namespace Statistics
} // End namespace itk

#endif /* ITKGREYLEVELCOOCCURRENCEMATRIX_HXX */

