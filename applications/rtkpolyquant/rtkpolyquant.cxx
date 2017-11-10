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

#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkPolyquantConeBeamReconstructionFilter.h"
#include "rtkNormalizedJosephBackProjectionImageFilter.h"

#ifdef RTK_USE_CUDA
  #include "itkCudaImage.h"
#endif

#include <itkImageFileWriter.h>

int main(int argc, char * argv[])
{
  GGO(rtkpolyquant, args_info);

  typedef float OutputPixelType;
  const unsigned int Dimension = 3;

#ifdef RTK_USE_CUDA
  typedef itk::CudaImage< OutputPixelType, Dimension > OutputImageType;
#else
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
#endif

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

  // SART reconstruction filter
  rtk::PolyquantConeBeamReconstructionFilter< OutputImageType >::Pointer polyquant =
      rtk::PolyquantConeBeamReconstructionFilter< OutputImageType >::New();

  // Set the forward and back projection filters
  polyquant->SetForwardProjectionFilter(args_info.fp_arg);
  polyquant->SetBackProjectionFilter(args_info.bp_arg);
  polyquant->SetInput( inputFilter->GetOutput() );
  if (args_info.signal_given)
    {
    polyquant->SetInput(1, phaseGating->GetOutput());
    polyquant->SetGeometry( phaseGating->GetOutputGeometry() );
    polyquant->SetGatingWeights( phaseGating->GetGatingWeightsOnSelectedProjections() );
    }
  else
    {
    polyquant->SetInput(1, reader->GetOutput());
    polyquant->SetGeometry( geometryReader->GetOutputObject() );
    }
  polyquant->SetNumberOfIterations( args_info.niterations_arg );
  polyquant->SetNumberOfProjectionsPerSubset( args_info.nprojpersubset_arg );
  polyquant->SetLambda( args_info.lambda_arg );
  polyquant->SetDisableDisplacedDetectorFilter(args_info.nodisplaced_flag);

  itk::TimeProbe totalTimeProbe;
  if(args_info.time_flag)
    {
    std::cout << "Recording elapsed time... " << std::endl << std::flush;
    totalTimeProbe.Start();
    }
  if(args_info.positivity_flag)
    {
    polyquant->SetEnforcePositivity(true);
    }

  TRY_AND_EXIT_ON_ITK_EXCEPTION( polyquant->Update() )

  if(args_info.time_flag)
    {
    polyquant->PrintTiming(std::cout);
    totalTimeProbe.Stop();
    std::cout << "It took...  " << totalTimeProbe.GetMean() << ' ' << totalTimeProbe.GetUnit() << std::endl;
    }

  // Write
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( polyquant->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() )

  return EXIT_SUCCESS;
}
