#include "rtkpolyquant.h"
#include "alg_polyquant.h"
#include "calc_polyquant.h"
#include "bit_reversal.h"

void os_polyquant(paramType &param,ctSystemType ctSystem)
{
	itk::TimeProbe subSetProbe;
	std::cout << "Number of iterations = " << param.nIter << std::endl;
	std::cout << "Number of projections = " << param.nProj << std::endl;
	volType ySub;
	int ind; // = sizeof(indArray)/sizeof(*indArray);
	float t = 1; float t1;
	ctSystemType subSetSystem = ctSystem;
	paramType subSetParam = param;
	volType volNow = param.volOld;//->GetOutput();
	subSetParam.recon = param.volOld;
	std::cout << "Running Polyquant reconstruction..." << std::endl;
	param.y->Update();
	std::vector<int> indArray;

  for (int j = 0; j<param.nProj/param.nSplit; j++)
  {
  	indArray.push_back(j*param.nSplit);
  }
  volType emptyProj = calc_subset_proj(param,indArray);
  
  multiplyType::Pointer multFilter = multiplyType::New();
  multFilter->SetInput1(emptyProj);
  multFilter->SetConstant2(0);
  multFilter->Update();
  subSetSystem.forProj->SetInput(multFilter->GetOutput());
  //subSetSystem.forProj->ReleaseDataFlagOn();
	subSetSystem.backProj->SetInput(param.volOld);
  volType grad;
  subtractType::Pointer derivUpdate = subtractType::New();
  multiplyType::Pointer stepSizeFilter = multiplyType::New();
  stepSizeFilter->SetConstant2(param.stepSize);
  OutputImageType::RegionType largestRegion = param.volOld->GetLargestPossibleRegion();
  largestRegion = param.volOld->GetLargestPossibleRegion();
  
  if (param.nSplit == 1)
  {
    subSetSystem.forProj->SetGeometry(ctSystem.geom);
    subSetSystem.backProj->SetGeometry(ctSystem.geom);
  }
  
  multiplyType::Pointer fistaMult = multiplyType::New();
  addType::Pointer fistaAdd = addType::New();
  subtractType::Pointer fistaSub = subtractType::New();
	for (int k = 0; k < param.nIter; k++)
  {
  	derivUpdate->SetInput1(subSetParam.recon);
  	//subSetProbe.Reset();
  	subSetProbe.Start();
  	if (param.nSplit > 1)
  	{
  	  indArray.clear();
    	ind = (k%param.nSplit);
    	ind = calc_bit_reversal(ind,param.nSplit);
    	// calculate a subset
      for (int j = 0; j<param.nProj/param.nSplit; j++)
      {
      	indArray.push_back(ind+j*param.nSplit);
      }
      subSetParam.y = calc_subset_proj(param,indArray);
      subSetSystem.geom = calc_subset_geom(ctSystem,indArray);
      subSetSystem.forProj->SetGeometry(subSetSystem.geom);
      subSetSystem.backProj->SetGeometry(subSetSystem.geom);
    }
    grad = grad_polyquant(subSetParam,subSetSystem);
    grad->DisconnectPipeline();
    stepSizeFilter->SetInput1(grad);
    derivUpdate->SetInput2(stepSizeFilter->GetOutput());
    derivUpdate->GetOutput()->SetRequestedRegion( largestRegion );
    volNow = proj_nn(derivUpdate->GetOutput(),param);
    volNow->SetRegions(largestRegion);
    volNow->SetOrigin(param.volOld->GetOrigin());
    volNow->SetSpacing(param.volOld->GetSpacing());
    subSetProbe.Stop();
    print_stats(k,param,subSetProbe);
  	//std::cout << "\tCompleted in " << subSetProbe.GetTotal()
  	//															 << subSetProbe.GetUnit() << std::endl;
  	if (param.accelerate)
  	{
  	  t1 = 0.5*(1+sqrt(1+4*t*t));
  	  fistaSub->SetInput1(volNow);
  	  fistaSub->SetInput2(subSetParam.volOld);
  	  fistaSub->Update();
  	  fistaMult->SetInput1(fistaSub->GetOutput());
  	  fistaMult->SetConstant2((t-1)/t1);
  	  fistaMult->Update();
  	  fistaAdd->SetInput1(volNow);
  	  fistaAdd->SetInput2(fistaMult->GetOutput());
  	  fistaAdd->Update();
      subSetParam.recon = fistaAdd->GetOutput();
      subSetParam.volOld = volNow;
      t = t1;
  	}
  	else
  	{
  	  subSetParam.recon = volNow;
  	}			 				
  }
  std::cout << std::endl;
  // Output gradient for testing
  param.recon = subSetParam.recon;
}

geomType calc_subset_geom(ctSystemType &ctSystem,std::vector<int> indArray)
{
  // This function selects a subset of the full geometry file
  geomType newGeom = rtk::ThreeDCircularProjectionGeometry::New();
  for (int k = 0; k<indArray.size(); k++)
  {
    newGeom->AddProjection(ctSystem.geom->GetMatrices()[indArray[k]]);
  }
  return newGeom;
}

volType calc_subset_proj(paramType &param,std::vector<int> indArray)
{
	typedef rtk::ConstantImageSource< OutputImageType > SourceType;
  SourceType::Pointer source = SourceType::New();
  source->SetInformationFromImage(param.y);
  OutputImageType::SizeType outputSize = param.y->GetLargestPossibleRegion().GetSize();
  outputSize[2] = indArray.size();
  source->SetSize(outputSize);
  source->SetConstant(0);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(source->Update())

  // Fill in the outputGeometry and the output projections
  typedef itk::PasteImageFilter<OutputImageType> PasteType;
  PasteType::Pointer paste = PasteType::New();
  paste->SetSourceImage(param.y);
  paste->SetDestinationImage(source->GetOutput());
  
  OutputImageType::RegionType sourceRegion;
  OutputImageType::IndexType destinationIndex;
  for (unsigned int i=0; i<indArray.size(); i++)
  {
    // If it is not the first projection, we need to use the output of
    // the paste filter as input
    if(i)
    {
      OutputImageType::Pointer pimg = paste->GetOutput();
      pimg->DisconnectPipeline();
      paste->SetDestinationImage(pimg);
    }
    sourceRegion = param.y->GetLargestPossibleRegion();
    sourceRegion.SetIndex(2, indArray[i]);
    sourceRegion.SetSize(2, 1);
    paste->SetSourceRegion(sourceRegion);

    destinationIndex = param.y->GetLargestPossibleRegion().GetIndex();
    destinationIndex[2] = i;
    paste->SetDestinationIndex(destinationIndex);

    TRY_AND_EXIT_ON_ITK_EXCEPTION( paste->Update() )
  }
  return paste->GetOutput();
}

volType proj_nn(volType in,paramType &param)
{
  thresholdType::Pointer threshFilt = thresholdType::New();
  if (param.gamma > 0)
  {
    tvType::Pointer tvFilt = tvType::New();
    tvFilt->SetInput(in);
    tvFilt->SetNumberOfIterations(3);
    tvFilt->SetGamma(param.gamma);
    tvFilt->Update();
    in = tvFilt->GetOutput();
    in->DisconnectPipeline();
  }
  threshFilt->SetOutsideValue(0);
  threshFilt->ThresholdOutside(0,param.up);
  threshFilt->SetInput(in);
  threshFilt->Update();
  return threshFilt->GetOutput();
}

void print_stats(int iter,paramType &param,itk::TimeProbe &subSetProbe)
{
  int barWidth = 40;
  float progress = (float)(iter+1)/(float)param.nIter;
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  float remain = subSetProbe.GetMean()*(float)(param.nIter-iter+1);
  std::cout << std::setprecision(3) << "] (" << iter+1 << "/" << param.nIter
            << ") avg_t=" << subSetProbe.GetMean()
            << "s tot_t=" << subSetProbe.GetTotal()
            << "s rem_t=" << remain << "s \r";
  std::cout.flush();
}
