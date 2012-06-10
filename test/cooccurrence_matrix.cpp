#include "itkCooccurrenceMatrix.h"
#include "itkNumericTraits.h"

typedef itk::Statistics::CooccurrenceMatrix< unsigned int > CooccurrenceMatrix;

bool assertSizeIs(const CooccurrenceMatrix * const m, const CooccurrenceMatrix::SizeType s)
{
  if(m->GetSize() != s)
  {
    std::cout << "The size of the matrix is " << m->GetSize() << ", it should be " << s << "." << std::endl;
    return false;
  }
  return true;
}

bool assertTotalCountIs(const CooccurrenceMatrix * const m, const CooccurrenceMatrix::MeasurementType s)
{
  if(m->GetTotalCount() != s)
  {
    std::cout << "The total count of the matrix is " << m->GetTotalCount() << ", it should be " << s << "." << std::endl;
    return false;
  }
  return true;
}

bool assertIsZeroEverywhere(const CooccurrenceMatrix * const m)
{
  typename CooccurrenceMatrix::ConstIterator it = m->Begin(), end = m->End();
  typename CooccurrenceMatrix::MeasurementType zero = itk::NumericTraits< CooccurrenceMatrix::MeasurementType >::Zero;

  while(it != end)
  {
    if((*it) != zero)
    {
      std::cout << "The matrix does not contains only zeros." << std::endl;
      return false;
    }
    ++it;
  }
  return true;
}

int main(void) {
  CooccurrenceMatrix::Pointer matrix = CooccurrenceMatrix::New();
  if(!assertSizeIs(matrix, 0))
    exit(-1);

  if(!assertTotalCountIs(matrix, 0))
    exit(-1);

  matrix->SetSize(16);
  if(!assertSizeIs(matrix, 16))
    exit(-1);
  if(!assertTotalCountIs(matrix, 0))
    exit(-1);
  if(!assertIsZeroEverywhere(matrix))
    exit(-1);

  matrix->IncrementCounter(3, 3);
  if(assertIsZeroEverywhere(matrix))
    exit(-1);
  if(!assertTotalCountIs(matrix, 1))
    exit(-1);

  typename CooccurrenceMatrix::ConstIterator it = matrix->Begin(), begin = matrix->Begin(), end = matrix->End();
  typename CooccurrenceMatrix::IndexType i1, i2;
  while(it != end)
  {
    matrix->GetIndexes(it - begin, &i1, &i2);
    if(i1 == 3 && i2 == 3)
    {
      if(*it != 1)
      {
        std::cout << "Cooccurrence (3, 3) has a value of " << (*it) << ", shoule have a value of 1." << std::endl;
        exit(-1);
      }
      break;
    }
    ++it;
  }

  matrix->SetToZero();
  if(!assertTotalCountIs(matrix, 0))
    exit(-1);
  if(!assertIsZeroEverywhere(matrix))
    exit(-1);

  std::cout << std::endl << "Test passed." << std::endl;
}
