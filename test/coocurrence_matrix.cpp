#include <itkCoocurrenceMatrix.h>

typedef itk::Statistics::CoocurrenceMatrix< unsigned int, float > CoocurrenceMatrixType;

int main(void) {
  CoocurrenceMatrixType matrix;
  matrix.SetSize(16); 

}
