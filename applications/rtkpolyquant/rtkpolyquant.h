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
#include "rtkDisplacedDetectorImageFilter.h"

#include <itkTimeProbe.h>
#include <itkImageFileWriter.h>
#include <itkSubtractImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkExpImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkMaskImageFilter.h>

#include<fstream>

//#include <itkMatrix>

// typedef itk::Image<float,3> volumeType projType;
typedef itk::Image< float, 3 > OutputImageType;
typedef OutputImageType::Pointer volType;
//typedef itk::Image< float, 3 > OutputImageType;
typedef rtk::ThreeDCircularProjectionGeometry::Pointer geomType;

typedef rtk::DisplacedDetectorImageFilter<OutputImageType> ddType;

//typedef rtk::DisplacedDetectorImageFilter<OutputImageType>::Pointer ddFilterType;
typedef rtk::ForwardProjectionImageFilter<OutputImageType,OutputImageType>::Pointer fType;
typedef rtk::BackProjectionImageFilter<OutputImageType,OutputImageType>::Pointer bType;
typedef itk::SubtractImageFilter<OutputImageType,OutputImageType> subtractType;
typedef itk::AddImageFilter<OutputImageType,OutputImageType> addType;
typedef itk::MultiplyImageFilter<OutputImageType,OutputImageType> multiplyType;
typedef itk::ThresholdImageFilter<OutputImageType> thresholdType;
typedef itk::BinaryThresholdImageFilter<OutputImageType,OutputImageType> maskType;
typedef itk::ExpImageFilter<OutputImageType,OutputImageType> expType;
typedef itk::DivideImageFilter<OutputImageType,OutputImageType,OutputImageType> divType;
typedef itk::MaskImageFilter<OutputImageType,OutputImageType> maskingType;
																				 
struct paramType
{
  int nIter;
  int nFac;
  int nSplit;
  int nProj;
  float stepSize;
  float up;  // up is the maximum allowable value
  bool accelerate;
  volType volOld;
  volType recon;
  volType y;  // y is the measurements
  std::vector<float> spectrum;
  std::vector<float> facKnee;
  std::vector< std::vector<float> > knee;
  //itk::Matrix kneeData;
};

struct ctSystemType
{
  fType forProj;
  bType backProj;
  geomType geom;
};
