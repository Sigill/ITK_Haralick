#include "itkGreyLevelCooccurrenceMatrix.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Statistics
{

template< class TGreyLevelCooccurrenceMatrix = GreyLevelCooccurrenceMatrix< unsigned int >, class TFrequencyType = float >
class ITK_EXPORT GreyLevelCooccurrenceMatrixNormalizerFilter:public ProcessObject
{ 
public:
  typedef GreyLevelCooccurrenceMatrixNormalizerFilter Self;
  typedef ProcessObject                    Superclass;
  typedef SmartPointer< Self >             Pointer;
  typedef SmartPointer< const Self >       ConstPointer;

  typedef TGreyLevelCooccurrenceMatrix GreyLevelCooccurrenceMatrixType;
  typedef TFrequencyType FrequencyType;
  typedef GreyLevelCooccurrenceMatrix< FrequencyType > NormalizedGreyLevelCooccurrenceMatrixType;

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( OutputIsFloatingPointCheck,
      ( Concept::IsFloatingPoint< FrequencyType > ) );
#endif

  itkTypeMacro(GreyLevelCooccurrenceMatrixNormalizerFilter, ProcessObject);

  itkNewMacro(Self);

  GreyLevelCooccurrenceMatrixNormalizerFilter()
  {
    itkDebugMacro("Constructor called");

    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);

    this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );
  }

  const GreyLevelCooccurrenceMatrixType *GetInput() const
  { 
    return static_cast< const GreyLevelCooccurrenceMatrixType * >( this->GetPrimaryInput() );
  }

  void SetInput(const GreyLevelCooccurrenceMatrixType * cooccurrenceMatrix)
  {
    this->ProcessObject::SetNthInput( 0, const_cast< GreyLevelCooccurrenceMatrixType * >(cooccurrenceMatrix) );
  }

  const NormalizedGreyLevelCooccurrenceMatrixType * GetOutput() const
  {
    const NormalizedGreyLevelCooccurrenceMatrixType *output =
      static_cast< const NormalizedGreyLevelCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

    return output;
  }

  DataObject::Pointer MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
  { 
    typename NormalizedGreyLevelCooccurrenceMatrixType::Pointer output = NormalizedGreyLevelCooccurrenceMatrixType::New();
      return static_cast< DataObject * >( output.GetPointer() );
  }

  void GenerateData(void)
  {
    itkDebugMacro("Starting normalization");
    const GreyLevelCooccurrenceMatrixType * const input = this->GetInput();
    typename GreyLevelCooccurrenceMatrixType::MeasurementType totalCount = input->GetTotalCount();
    NormalizedGreyLevelCooccurrenceMatrixType *output =
      static_cast< NormalizedGreyLevelCooccurrenceMatrixType * >( this->ProcessObject::GetOutput(0) );

    output->SetSize(input->GetSize());

    if(totalCount == 0)
      {
      itkDebugMacro("Input GreyLevelCooccurrenceMatrix is empty (TotalCounter = 0)");
      output->SetToZero();
      }
    else
      {
      typename GreyLevelCooccurrenceMatrixType::ConstIterator iit = input->Begin(), iit_end = input->End();
      typename NormalizedGreyLevelCooccurrenceMatrixType::Iterator oit = output->Begin();

      typename NormalizedGreyLevelCooccurrenceMatrixType::MeasurementType nTotalCount = totalCount;

      while(iit != iit_end)
        {
          *oit = *iit / nTotalCount;
          ++iit;
          ++oit;
        }

      output->SetTotalCount(NumericTraits< typename NormalizedGreyLevelCooccurrenceMatrixType::MeasurementType >::One);
      }

    itkDebugMacro("Normalization done");
  }

};
} // end of namespace Statistics
} // end of namespace itk

