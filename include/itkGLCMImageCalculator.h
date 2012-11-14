#ifndef ITKGLCMIMAGECALGULATOR_H
#define ITKGLCMIMAGECALGULATOR_H

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkVectorContainer.h"

namespace itk
{
template< typename TInputImage, typename TGLCMType >
class ITK_EXPORT GLCMImageCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef GLCMImageCalculator                          Self;
  typedef Object                                       Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  typedef TInputImage                                  ImageType;
  typedef typename ImageType::OffsetType               OffsetType;
  typedef typename ImageType::Pointer                  ImagePointer;
  typedef typename ImageType::ConstPointer             ImageConstPointer;
  typedef typename ImageType::PixelType                PixelType;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::RegionType               RegionType;
  typedef typename ImageType::SizeType                 RadiusType;

  typedef TGLCMType                                    GLCMType;
  typedef typename GLCMType::Pointer                   GLCMPointer;
  typedef typename GLCMType::ConstPointer              GLCMConstPointer;
  typedef typename GLCMType::SizeType                  GLCMSizeType;

  typedef VectorContainer< unsigned char, OffsetType > OffsetVector;
  typedef typename OffsetVector::Pointer               OffsetVectorPointer;
  typedef typename OffsetVector::ConstPointer          OffsetVectorConstPointer;

  itkNewMacro(Self)

  itkTypeMacro(GLCMImageCalculator, Object)

  itkSetConstObjectMacro(Image, ImageType)

  itkSetMacro(MatrixSize, GLCMSizeType)
  itkGetConstMacro(MatrixSize, GLCMSizeType)

  itkGetConstObjectMacro(CooccurrenceMatrix, GLCMType)

  void Compute(void);

  void SetRegion(const RegionType & region);

  void ResetMatrix(void);

  itkGetConstObjectMacro(Offsets, OffsetVector)
  void SetOffsets(const OffsetVector * os);
  void SetOffset(const OffsetType o);

protected:
  GLCMImageCalculator();
  virtual ~GLCMImageCalculator() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  GLCMImageCalculator(const Self &); //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  void ComputeOffsetsMinRadius(void);

  ImageConstPointer m_Image;

  GLCMPointer m_CooccurrenceMatrix;
  GLCMSizeType m_MatrixSize;

  OffsetVectorConstPointer m_Offsets;
  RadiusType m_OffsetsMinRadius;

  RegionType m_Region;
  bool       m_RegionSetByUser;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGLCMImageCalculator.hxx"
#endif

#endif /* ITKGLCMIMAGECALGULATOR_H */
