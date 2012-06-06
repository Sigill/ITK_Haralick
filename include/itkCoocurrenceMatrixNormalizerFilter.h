#include "itkCoocurrenceMatrix.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Statistics
{

template< class TCoocurrenceMatrix = CoocurrenceMatrix< unsigned int >, class TFrequencyType = float >
class ITK_EXPORT CoocurrenceMatrixNormalizerFilter:public ProcessObject
{ 
public:
  typedef CoocurrenceMatrixNormalizerFilter Self;
  typedef ProcessObject                    Superclass;
  typedef SmartPointer< Self >             Pointer;
  typedef SmartPointer< const Self >       ConstPointer;

  typedef TCoocurrenceMatrix CoocurrenceMatrixType;
  typedef TFrequencyType FrequencyType;
  typedef CoocurrenceMatrix< FrequencyType > NormalizedCoocurrenceMatrixType;

  itkTypeMacro(CoocurrenceMatrixNormalizerFilter, ProcessObject);

  itkNewMacro(Self);

  CoocurrenceMatrixNormalizerFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);

    this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );
  }

  const CoocurrenceMatrixType *GetInput() const
  { 
    return static_cast< const CoocurrenceMatrixType * >( this->GetPrimaryInput() );
  }

  void SetInput(const CoocurrenceMatrixType * coocurrenceMatrix)
  {
    this->ProcessObject::SetNthInput( 0, const_cast< CoocurrenceMatrixType * >(coocurrenceMatrix) );
  }

  const NormalizedCoocurrenceMatrixType * GetOutput() const
  {
    const NormalizedCoocurrenceMatrixType *output =
      static_cast< const NormalizedCoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

    return output;
  }

  DataObject::Pointer MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
  { 
    typename NormalizedCoocurrenceMatrixType::Pointer output = NormalizedCoocurrenceMatrixType::New();
      return static_cast< DataObject * >( output.GetPointer() );
  }

  void GenerateData(void)
  {
    const CoocurrenceMatrixType * const input = this->GetInput();
    typename CoocurrenceMatrixType::CounterType totalCount = input->GetTotalCount();
    NormalizedCoocurrenceMatrixType *output =
      static_cast< NormalizedCoocurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

    output->SetSize(input->GetSize());

    if(totalCount == 0)
      {
      output->SetToZero();
      }
    else
      {
      typename CoocurrenceMatrixType::ConstIterator iit = input->Begin(), iit_end = input->End();
      typename NormalizedCoocurrenceMatrixType::Iterator oit = output->Begin();

      typename NormalizedCoocurrenceMatrixType::CounterType nTotalCount = totalCount;

      while(iit != iit_end)
        {
          *oit = *iit / nTotalCount;
          ++iit;
          ++oit;
        }

      output->SetTotalCount(NumericTraits< typename NormalizedCoocurrenceMatrixType::CounterType >::One);
      }

  }

};
} // end of namespace Statistics
} // end of namespace itk

