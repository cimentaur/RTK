#include "rtkJosephForwardProjectionImageFilter.h"
#ifdef RTK_USE_CUDA
#include "rtkCudaForwardProjectionImageFilter.h"
#endif
#include "rtkRayCastInterpolatorForwardProjectionImageFilter.h"
#include "rtkJosephBackProjectionImageFilter.h"
#include "rtkNormalizedJosephBackProjectionImageFilter.h"
#ifdef RTK_USE_CUDA
#include "rtkCudaBackProjectionImageFilter.h"
#include "rtkCudaRayCastBackProjectionImageFilter.h"
#endif

// #include "rtkPhaseReader.h"
#include "rtkThreeDCircularProjectionGeometryXMLFile.h"

#include <itkTimeProbe.h>
#include <itkImageFileWriter.h>

// typedef itk::Image<float,3> volumeType projType;
typedef itk::Image< float, 3 > OutputImageType;
typedef OutputImageType::Pointer volType;
//typedef itk::Image< float, 3 > OutputImageType;
typedef rtk::ThreeDCircularProjectionGeometry::Pointer geomType;

typedef rtk::ForwardProjectionImageFilter<OutputImageType,OutputImageType>::Pointer fType;
typedef rtk::BackProjectionImageFilter<OutputImageType,OutputImageType>::Pointer bType;


																				 
struct paramType
{
  int nIter;
  int nFac;
  int nSplit;
  int nProj;
  bool accelerate;
  volType volOld;
  volType y;  // y is the measurements
};

struct ctSystemType
{
  fType forProj;
  bType backProj;
  geomType geom;
};