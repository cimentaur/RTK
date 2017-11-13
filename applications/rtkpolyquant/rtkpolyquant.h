#include "rtkJosephForwardProjectionImageFilter.h"
#ifdef RTK_USE_CUDA
#include "rtkCudaForwardProjectionImageFilter.h"
#endif
#include "rtkRayCastInterpolatorForwardProjectionImageFilter.h"
#include "rtkJosephBackProjectionImageFilter.h"
#include "rtkNormalizedJosephBackProjectionImageFilter.h"
#ifdef RTK_USE_CUDA
#  include "rtkCudaBackProjectionImageFilter.h"
#  include "rtkCudaRayCastBackProjectionImageFilter.h"
#endif
//#include "rtkThreeDCircularProjectionGeometryXMLFile.h"

// typedef itk::Image<float,3> volumeType projType;
typedef itk::Image< float, 3 > OutputImageType;
//typedef itk::Image< float, 3 > OutputImageType;
typedef rtk::ThreeDCircularProjectionGeometry::Pointer geomType;

typedef rtk::ForwardProjectionImageFilter<OutputImageType,OutputImageType>::Pointer fType;
typedef rtk::BackProjectionImageFilter<OutputImageType,OutputImageType>::Pointer bType;
																				 
struct paramType
{
  int nIter;
  int nFac;
  int nSplit;
  bool accelerate;
};

struct ctSystemType
{
  fType forProj;
  bType backProj;
  geomType geom;
};