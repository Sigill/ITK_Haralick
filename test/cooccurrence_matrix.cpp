#include <itkCooccurrenceMatrix.h>

typedef itk::Statistics::CooccurrenceMatrix< unsigned int > CooccurrenceMatrixType;

int main(void) {
  CooccurrenceMatrixType matrix;
  matrix.SetSize(16); 

}
