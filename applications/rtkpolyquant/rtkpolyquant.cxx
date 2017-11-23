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

#include "rtkpolyquant_ggo.h"
#include "rtkGgoFunctions.h"

#include "rtkpolyquant.h"
#include "alg_polyquant.h"

//#include "rtkThreeDCircularProjectionGeometryXMLFile.h"

#ifdef RTK_USE_CUDA
  #include "itkCudaImage.h"
#endif



int main(int argc, char * argv[])
{
  GGO(rtkpolyquant, args_info);

  typedef float OutputPixelType;
  const unsigned int Dimension = 3;

//#ifdef RTK_USE_CUDA
//  typedef itk::CudaImage< OutputPixelType, Dimension > OutputImageType;
//#else
//  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
//#endif

  // Projections reader
  typedef rtk::ProjectionsReader< OutputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  rtk::SetProjectionsReaderFromGgo<ReaderType, args_info_rtkpolyquant>(reader, args_info);

  // Geometry
  if(args_info.verbose_flag)
    std::cout << "Reading geometry information from "
              << args_info.geometry_arg
              << "..."
              << std::endl;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( geometryReader->GenerateOutputInformation() )

  // Create input: either an existing volume read from a file or a blank image
  itk::ImageSource< OutputImageType >::Pointer inputFilter;
  if(args_info.input_given)
    {
    // Read an existing image to initialize the volume
    typedef itk::ImageFileReader<  OutputImageType > InputReaderType;
    InputReaderType::Pointer inputReader = InputReaderType::New();
    inputReader->SetFileName( args_info.input_arg );
    inputFilter = inputReader;
    }
  else
    {
    // Create new empty volume
    typedef rtk::ConstantImageSource< OutputImageType > ConstantImageSourceType;
    ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
    rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_rtkpolyquant>(constantImageSource, args_info);
    inputFilter = constantImageSource;
    }

  // CT system definition
  ctSystemType ctSystem;
  paramType param;

  // Set the forward and back projection filters
  switch(args_info.fp_arg)
  {
  case(fp_arg_Joseph):
    ctSystem.forProj = rtk::JosephForwardProjectionImageFilter<OutputImageType, OutputImageType>::New();
    break;
  case(fp_arg_RayCastInterpolator):
    ctSystem.forProj = rtk::RayCastInterpolatorForwardProjectionImageFilter<OutputImageType, OutputImageType>::New();
    break;
  case(fp_arg_CudaRayCast):
#ifdef RTK_USE_CUDA
    ctSystem.forProj = rtk::CudaForwardProjectionImageFilter<OutputImageType,OutputImageType>::New();
    dynamic_cast<rtk::CudaForwardProjectionImageFilter<OutputImageType,OutputImageType>*>(ctSystem.forProj.GetPointer())->SetStepSize(args_info.step_arg);
#else
    std::cerr << "The program has not been compiled with cuda option" << std::endl;
    return EXIT_FAILURE;
#endif
    break;
  default:
    std::cerr << "Unhandled --method value." << std::endl;
    return EXIT_FAILURE;
  }
  switch(args_info.bp_arg)
  {
    case(bp_arg_VoxelBasedBackProjection):
      ctSystem.backProj = rtk::BackProjectionImageFilter<OutputImageType,OutputImageType>::New();
      break;
    case(bp_arg_Joseph):
      ctSystem.backProj = rtk::JosephBackProjectionImageFilter<OutputImageType,OutputImageType>::New();
      break;
    case(bp_arg_NormalizedJoseph):
      ctSystem.backProj = rtk::NormalizedJosephBackProjectionImageFilter<OutputImageType,OutputImageType>::New();
      break;
    case(bp_arg_CudaBackProjection):
#ifdef RTK_USE_CUDA
      ctSystem.backProj = rtk::CudaBackProjectionImageFilter::New();
#else
      std::cerr << "The program has not been compiled with cuda option" << std::endl;
      return EXIT_FAILURE;
#endif
      break;
    case(bp_arg_CudaRayCast):
#ifdef RTK_USE_CUDA
      ctSystem.backProj = rtk::CudaRayCastBackProjectionImageFilter::New();
#else
      std::cerr << "The program has not been compiled with cuda option" << std::endl;
      return EXIT_FAILURE;
#endif
      break;
    default:
    std::cerr << "Unhandled --method value." << std::endl;
    return EXIT_FAILURE;
  }

  ctSystem.geom = geometryReader->GetOutputObject();
  param.y = reader->GetOutput();
  param.volOld = inputFilter->GetOutput();
  param.nIter = args_info.niterations_arg;
  param.nSplit = args_info.nsplit_arg;
  param.stepSize = args_info.lambda_arg*(float)param.nSplit;
  param.gamma = args_info.gamma_arg;
  param.nProj = ctSystem.geom->GetGantryAngles().size();
  param.accelerate = (bool)args_info.accelerate_flag;
  param.up = 10;
  std::cout << "Found spectrum and it is " << args_info.spectrumfile_arg << std::endl;
  std::ifstream spectrumFile(args_info.spectrumfile_arg);
  // test file open   
  if (spectrumFile)
  {        
	  float value;
	  // read the elements in the file into a vector  
	  while (spectrumFile >> value)
	  {
	    param.spectrum.push_back(value);
	  }
	  spectrumFile.close();
	  //std::cout << "Found spectrum and it is " << param.spectrum << std::endl;
  }
  else
  {
  	std::cout << "I cannot find the spectrum file .o0 FAILS!" << std::endl;
  	return EXIT_FAILURE;
  }
  std::cout << "Found knee and it is " << args_info.kneefile_arg << std::endl;
  std::cout << "Spectrum length is " << param.spectrum.size() << std::endl;
  std::ifstream kneeFile(args_info.kneefile_arg);
  // test file open   
  if (kneeFile)
  {
    float value;
    // insert the knee factors into the vector
    for (int k = 0; k<2; k++)
    {
      kneeFile >> value;
      param.facKnee.push_back(value);
    }      
	  // read the elements in the file into a vector  
	  for (int k = 0; k<param.spectrum.size(); k++)
	  {
	    std::vector<float> row; 
      for (int j = 0; j<5; j++)
      {
        kneeFile >> value;
        row.push_back(value);
      }
      param.knee.push_back(row);
	  }
	  //std::cout << "Found knee and it is " << param.knee << std::endl;
	  kneeFile.close();
  }
  else
  {
  	std::cout << "I cannot find the spectrum file .o0 FAILS!" << std::endl;
  	return EXIT_FAILURE;
  }  

  itk::TimeProbe totalTimeProbe;

  std::cout << "Recording elapsed time... " << std::endl << std::flush;
  totalTimeProbe.Start();

  // Perform the update
  os_polyquant(param,ctSystem);

  
  //TRY_AND_EXIT_ON_ITK_EXCEPTION( polyquant->Update() )

  //polyquant->PrintTiming(std::cout);
  totalTimeProbe.Stop();
  std::cout << "It took...  " << totalTimeProbe.GetMean() << ' ' << totalTimeProbe.GetUnit() << std::endl;

  // Write
  if(args_info.verbose_flag)
    std::cout << "Writing... " << std::endl;
  itk::TimeProbe writeProbe;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput(param.recon);
  //writer->GetInput(param.rec);
  writeProbe.Start();
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() )
  writeProbe.Stop();
  if(args_info.verbose_flag)
    std::cout << " done in "
              << writeProbe.GetMean() << ' ' << writeProbe.GetUnit()
              << '.' << std::endl;

  return EXIT_SUCCESS;
}
