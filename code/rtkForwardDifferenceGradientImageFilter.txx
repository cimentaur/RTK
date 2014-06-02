/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef __rtkForwardDifferenceGradientImageFilter_txx
#define __rtkForwardDifferenceGradientImageFilter_txx
#include "rtkForwardDifferenceGradientImageFilter.h"

#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include <itkImageRegionIterator.h>
#include <itkForwardDifferenceOperator.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkOffset.h>
#include <itkProgressReporter.h>

namespace rtk
{
  using namespace itk;
//
// Constructor
//
template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::ForwardDifferenceGradientImageFilter()
{
  // default boundary condition
  m_BoundaryCondition = new ZeroFluxNeumannBoundaryCondition<TInputImage>();

  // default behaviour is to take into account both spacing and direction
  this->m_UseImageSpacing   = true;
  this->m_UseImageDirection = true;

  // default behaviour is to process all dimensions
  for (int dim = 0; dim < TInputImage::ImageDimension; dim++)
    {
    m_DimensionsProcessed[dim] = true;
    }
}

//
// Destructor
//
template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::~ForwardDifferenceGradientImageFilter()
{
  delete m_BoundaryCondition;
}

// This should be handled by an itkMacro, but it doesn't seem to work with pointer types
template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
void
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::SetDimensionsProcessed(bool* DimensionsProcessed)
{
  bool Modified=false;
  for (int dim=0; dim<TInputImage::ImageDimension; dim++)
    {
    if (m_DimensionsProcessed[dim] != DimensionsProcessed[dim])
      {
      m_DimensionsProcessed[dim] = DimensionsProcessed[dim];
      Modified = true;
      }
    }
  if(Modified) this->Modified();
}

template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
void
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::OverrideBoundaryCondition(ImageBoundaryCondition< TInputImage >* boundaryCondition)
{
  m_BoundaryCondition = boundaryCondition;
}

template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
void
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer inputPtr =
    const_cast< InputImageType * >( this->GetInput() );
  OutputImagePointer outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // Build an operator so that we can determine the kernel size
  ForwardDifferenceOperator< OperatorValueType, InputImageDimension > oper;
  oper.SetDirection(0);
  oper.CreateDirectional();
  SizeValueType radius = oper.GetRadius()[0];

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius(radius);

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() ) )
    {
    inputPtr->SetRequestedRegion(inputRequestedRegion);
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion(inputRequestedRegion);

    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}

template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
void
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::BeforeThreadedGenerateData()
{
  Superclass::BeforeThreadedGenerateData();
  this->AllocateOutputs();
  CovariantVectorType zero;
  for (int dim=0; dim<TInputImage::ImageDimension; dim++)  zero[dim] = 0;
  this->GetOutput()->FillBuffer(zero);
}

template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
void
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  unsigned int    i;
  CovariantVectorType gradient;

  ConstNeighborhoodIterator< InputImageType > nit;
  ImageRegionIterator< OutputImageType >      it;

  NeighborhoodInnerProduct< InputImageType, OperatorValueType,
                            OutputValueType > SIP;

  // Get the input and output
  OutputImageType *     outputImage = this->GetOutput();
  const InputImageType *inputImage  = this->GetInput();

  // Generate a list of indices of the dimensions to process
  std::vector<int> dimsToProcess;
  for (int dim = 0; dim < TInputImage::ImageDimension; dim++)
    {
    if(m_DimensionsProcessed[dim]) dimsToProcess.push_back(dim);
    }

  // Set up operators
  ForwardDifferenceOperator< OperatorValueType, InputImageDimension > op[InputImageDimension];

  for ( i = 0; i < InputImageDimension; i++ )
    {
    op[i].SetDirection(0);
    op[i].CreateDirectional();

    // Take into account the pixel spacing if necessary
    if ( m_UseImageSpacing == true )
      {
      if ( this->GetInput()->GetSpacing()[i] == 0.0 )
        {
        itkExceptionMacro(<< "Image spacing cannot be zero.");
        }
      else
        {
        op[i].ScaleCoefficients(1.0 / this->GetInput()->GetSpacing()[i]);
        }
      }
    }

  // Calculate iterator radius
  Size< InputImageDimension > radius;
  for ( i = 0; i < InputImageDimension; ++i )
    {
    radius[i]  = op[0].GetRadius()[0];
    }

  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType > bC;
  faceList = bC(inputImage, outputRegionForThread, radius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >::FaceListType::iterator fit;
  fit = faceList.begin();

  // support progress methods/callbacks
  ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  // Initialize the x_slice array
  nit = ConstNeighborhoodIterator< InputImageType >(radius, inputImage, *fit);

  std::slice          x_slice[InputImageDimension];
  const SizeValueType center = nit.Size() / 2;
  for ( i = 0; i < InputImageDimension; ++i )
    {
    x_slice[i] = std::slice( center - nit.GetStride(i) * radius[i],
                             op[i].GetSize()[0], nit.GetStride(i) );
    }

  // Process non-boundary face and then each of the boundary faces.
  // These are N-d regions which border the edge of the buffer.
  for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    nit = ConstNeighborhoodIterator< InputImageType >(radius,
                                                      inputImage, *fit);
    it = ImageRegionIterator< OutputImageType >(outputImage, *fit);
    nit.OverrideBoundaryCondition(m_BoundaryCondition);
    nit.GoToBegin();

    while ( !nit.IsAtEnd() )
      {
      for ( i = 0; i < dimsToProcess.size(); ++i )
        {
        gradient[dimsToProcess[i]] = SIP(x_slice[dimsToProcess[i]], nit, op[dimsToProcess[i]]);
        }

      // This method optionally performs a tansform for Physical
      // coordiantes and potential conversion to a different output
      // pixel type.
      this->SetOutputPixel( it, gradient );

      ++nit;
      ++it;
      progress.CompletedPixel();
      }
    }
}

template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
void
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::GenerateOutputInformation()
{
  // this methods is overloaded so that if the output image is a
  // VectorImage then the correct number of components are set.

  Superclass::GenerateOutputInformation();
  OutputImageType* output = this->GetOutput();

  if ( !output )
    {
    return;
    }
  if ( output->GetNumberOfComponentsPerPixel() != InputImageDimension )
    {
    output->SetNumberOfComponentsPerPixel( InputImageDimension );
    }
}


/**
 * Standard "PrintSelf" method
 */
template< typename TInputImage, typename TOperatorValueType, typename TOuputValue , typename TOuputImage >
void
ForwardDifferenceGradientImageFilter< TInputImage, TOperatorValueType, TOuputValue, TOuputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "UseImageSpacing: "
     << ( this->m_UseImageSpacing ? "On" : "Off" ) << std::endl;
  os << indent << "UseImageDirection = "
     << ( this->m_UseImageDirection ? "On" : "Off" ) << std::endl;
}
} // end namespace itk

#endif