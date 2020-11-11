#include "rtkTestConfiguration.h"

#include "itkImageFileWriter.h"

#ifdef USE_CUDA
#  include "rtkCudaFFTRampImageFilter.h"
#else
#  include "rtkFFTRampImageFilter.h"
#endif

/**
 * \file rtkrampfiltertest.cxx
 *
 * \brief Functional test for the ramp filter of the FDK reconstruction.
 *
 * \author Julien JOmier
 */

int
main(int, char **)
{
  constexpr unsigned int Dimension = 3;
  using PixelType = float;
#ifdef USE_CUDA
  using ImageType = itk::CudaImage<PixelType, Dimension>;
  using RampFilterType = rtk::CudaFFTRampImageFilter;
#else
  using ImageType = itk::Image<PixelType, Dimension>;
  using RampFilterType = rtk::FFTRampImageFilter<ImageType, ImageType>;
#endif

  ImageType::Pointer    image = ImageType::New();
  ImageType::RegionType region;
  ImageType::SizeType   size;
  size[0] = 64;
  size[1] = 64;
  size[2] = 64;
  region.SetSize(size);
  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(10);

  RampFilterType::Pointer rampFilter = RampFilterType::New();
  rampFilter->SetInput(image);

  try
  {
    rampFilter->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }


  // Check the results
  ImageType::IndexType index;
  index[0] = 3;
  index[1] = 21;
  index[2] = 26;

  float value = 0.132652;
  if (fabs(rampFilter->GetOutput()->GetPixel(index) - value) > 0.000001)
  {
    std::cout << "Output value #0 should be " << value << " found " << rampFilter->GetOutput()->GetPixel(index)
              << " instead." << std::endl;
    return EXIT_FAILURE;
  }

  // Testing the HannCutFrequency
  rampFilter->SetHannCutFrequency(0.8);
  rampFilter->Update();
  value = 0.149724;
  if (fabs(rampFilter->GetOutput()->GetPixel(index) - value) > 0.000001)
  {
    std::cout << "Output value #1 should be " << value << " found " << rampFilter->GetOutput()->GetPixel(index)
              << " instead." << std::endl;
    return EXIT_FAILURE;
  }

  // Testing the HanncutFrequencyY
  rampFilter->SetHannCutFrequencyY(0.1);
  rampFilter->Update();
  value = 0.150181;
  if (fabs(rampFilter->GetOutput()->GetPixel(index) - value) > 0.000001)
  {
    std::cout << "Output value #2 should be " << value << " found " << rampFilter->GetOutput()->GetPixel(index)
              << " instead." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Test PASSED! " << std::endl;
  return EXIT_SUCCESS;
}
