#ifndef ITKGLCMIMAGECALGULATOR_H
#define ITKGLCMIMAGECALGULATOR_H

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkVectorContainer.h"

namespace itk
{
namespace Statistics
{

/**
 * \class GLCMImageCalculator
 * \brief Used to compute a grey-level cooccurrence matrix
 * on a region of an image.
 *
 * GLCMImageCalculator is templated over the type of image
 * and the type of cooccurrence matrix used internally.
 *
 * The algorithm iterates over a region specified by the user
 * and fills the cooccurrence matrix based on the cooccurrences
 * of grey-levels found in this region at the specified offsets.
 */

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

  typedef TGLCMType                                    GLCMType;
  typedef typename GLCMType::Pointer                   GLCMPointer;
  typedef typename GLCMType::ConstPointer              GLCMConstPointer;
  typedef typename GLCMType::SizeType                  GLCMSizeType;

  typedef VectorContainer< unsigned char, OffsetType > OffsetVectorType;
  typedef typename OffsetVectorType::Pointer           OffsetVectorPointer;
  typedef typename OffsetVectorType::ConstPointer      OffsetVectorConstPointer;

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( InputHasNumericTraitsCheck, 
    ( Concept::HasNumericTraits< typename ImageType::PixelType > ) );
#endif

  itkNewMacro(Self)

  itkTypeMacro(GLCMImageCalculator, Object)

  itkSetConstObjectMacro(Image, ImageType)

  itkSetMacro(MatrixSize, GLCMSizeType)
  itkGetConstMacro(MatrixSize, GLCMSizeType)

  itkGetConstObjectMacro(CooccurrenceMatrix, GLCMType)

  void Compute(void);

  void SetRegion(const RegionType & region);

  void ResetMatrix(void);

  itkGetConstObjectMacro(Offsets, OffsetVectorType)
  void SetOffsets(const OffsetVectorType * os);
  void SetOffset(const OffsetType o);

protected:
  GLCMImageCalculator();
  virtual ~GLCMImageCalculator() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  GLCMImageCalculator(const Self &); //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  ImageConstPointer        m_Image;

  GLCMPointer              m_CooccurrenceMatrix;
  GLCMSizeType             m_MatrixSize;

  OffsetVectorConstPointer m_Offsets;

  RegionType               m_Region;
  bool                     m_RegionSetByUser;
};
} // end namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGLCMImageCalculator.hxx"
#endif

#endif /* ITKGLCMIMAGECALGULATOR_H */
