#include "itkCooccurrenceMatrix.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Statistics
{

template< class TCooccurrenceMatrix = CooccurrenceMatrix< unsigned int >, class TFrequencyType = float >
class ITK_EXPORT CooccurrenceMatrixNormalizerFilter:public ProcessObject
{ 
public:
  typedef CooccurrenceMatrixNormalizerFilter Self;
  typedef ProcessObject                    Superclass;
  typedef SmartPointer< Self >             Pointer;
  typedef SmartPointer< const Self >       ConstPointer;

  typedef TCooccurrenceMatrix CooccurrenceMatrixType;
  typedef TFrequencyType FrequencyType;
  typedef CooccurrenceMatrix< FrequencyType > NormalizedCooccurrenceMatrixType;

  itkTypeMacro(CooccurrenceMatrixNormalizerFilter, ProcessObject);

  itkNewMacro(Self);

  CooccurrenceMatrixNormalizerFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);

    this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );
  }

  const CooccurrenceMatrixType *GetInput() const
  { 
    return static_cast< const CooccurrenceMatrixType * >( this->GetPrimaryInput() );
  }

  void SetInput(const CooccurrenceMatrixType * cooccurrenceMatrix)
  {
    this->ProcessObject::SetNthInput( 0, const_cast< CooccurrenceMatrixType * >(cooccurrenceMatrix) );
  }

  const NormalizedCooccurrenceMatrixType * GetOutput() const
  {
    const NormalizedCooccurrenceMatrixType *output =
      static_cast< const NormalizedCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

    return output;
  }

  DataObject::Pointer MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
  { 
    typename NormalizedCooccurrenceMatrixType::Pointer output = NormalizedCooccurrenceMatrixType::New();
      return static_cast< DataObject * >( output.GetPointer() );
  }

  void GenerateData(void)
  {
    const CooccurrenceMatrixType * const input = this->GetInput();
    typename CooccurrenceMatrixType::MeasurementType totalCount = input->GetTotalCount();
    NormalizedCooccurrenceMatrixType *output =
      static_cast< NormalizedCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

    output->SetSize(input->GetSize());

    if(totalCount == 0)
      {
      output->SetToZero();
      }
    else
      {
      typename CooccurrenceMatrixType::ConstIterator iit = input->Begin(), iit_end = input->End();
      typename NormalizedCooccurrenceMatrixType::Iterator oit = output->Begin();

      typename NormalizedCooccurrenceMatrixType::MeasurementType nTotalCount = totalCount;

      while(iit != iit_end)
        {
          *oit = *iit / nTotalCount;
          ++iit;
          ++oit;
        }

      output->SetTotalCount(NumericTraits< typename NormalizedCooccurrenceMatrixType::MeasurementType >::One);
      }

  }

};
} // end of namespace Statistics
} // end of namespace itk

