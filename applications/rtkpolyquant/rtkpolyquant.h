#include "rtkRayCastInterpolatorForwardProjectionImageFilter.h"
#include "rtkJosephBackProjectionImageFilter.h"
#include "rtkThreeDCircularProjectionGeometryXMLFile.h"

typedef itk::Image<float, 3> volumeType projType;
typedef rtk::ThreeDCircularProjectionGeometryXMLFileReader geomType;
typedef rtk::JosephForwardProjectionImageFilter<OutputImageType,
																								OutputImageType> forProjType;
typedef rtk::JosephBackProjectionImageFilter<OutputImageType,
																						 OutputImageType> backProjType;

struct paramType
{
  int nIter nFac nSplit;
  bool accelerate;
}

struct ctSystemType
{
  forProjType forProj;
  backProjType backProj;
  geomType geom;
}