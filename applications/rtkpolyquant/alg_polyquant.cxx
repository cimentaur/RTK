#include "rtkpolyquant.h"
#include "alg_polyquant.h"
#include "calc_polyquant.h"

void os_polyquant(paramType &param,ctSystemType ctSystem)
{
	itk::TimeProbe subSetProbe;
	std::cout << "Number of iterations = " << param.nIter << std::endl;
	std::cout << "Number of projections = " << param.nProj << std::endl;
	volType ySub;
	int ind; // = sizeof(indArray)/sizeof(*indArray);
	ctSystemType subSetSystem = ctSystem;
	paramType subSetParam = param;
	volType volNow = param.volOld;
	
	std::vector<int> indArray;
  for (int j = 0; j<param.nProj/param.nSplit; j++)
  {
  	indArray.push_back(j*param.nSplit);
  }
  volType emptyProj = calc_subset_proj(param,indArray);
  typedef itk::MultiplyImageFilter<OutputImageType,OutputImageType> multiplyType;
  multiplyType::Pointer multFilter = multiplyType::New();
  multFilter->SetInput1(emptyProj);
  multFilter->SetConstant2(0);
  multFilter->Update();
  emptyProj = multFilter->GetOutput();
  std::cout << "Initialising measurements..." << std::endl;
  param.y->Update();
  subSetSystem.forProj->SetInput(emptyProj);
  subSetSystem.backProj->SetInput(param.volOld);
	for (int k = 0; k < param.nIter; k++)
  {
    std::cout << "Running iteration " << k+1 << " of " << param.nIter << std::endl;
  	indArray.clear();
  	subSetProbe.Reset();
  	subSetProbe.Start();
  	ind = (k%param.nSplit);
  	// calculate a subset
    for (int j = 0; j<param.nProj/param.nSplit; j++)
    {
    	indArray.push_back(ind+j*param.nSplit);
    }
    subSetParam.y = calc_subset_proj(param,indArray);
    subSetSystem.geom = calc_subset_geom(ctSystem,indArray);
    subSetSystem.forProj->SetGeometry(subSetSystem.geom);
    subSetSystem.backProj->SetGeometry(subSetSystem.geom);
    param.y = grad_polyquant(subSetParam,subSetSystem);	
    subSetProbe.Stop();
  	std::cout << "\tCompleted in " << subSetProbe.GetTotal()
  																 << subSetProbe.GetUnit() << std::endl;
  								 				
  }
}

int bit_reversal(int index, int max)
{
	// TODO: apply bit reversal strategy
	return index;
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