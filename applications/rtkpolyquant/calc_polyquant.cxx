#include "rtkpolyquant.h"
#include "calc_polyquant.h"

volType grad_polyquant(paramType &param,ctSystemType &ctSystem)
{
  // TODO: tidy and generalise -- store all temporary variable in 4D array
  // Calculate the forward projections
  thresholdType::Pointer threshFilt = thresholdType::New();
  threshFilt->SetOutsideValue(0);
  threshFilt->ThresholdOutside(0,1);
  threshFilt->SetInput(param.volOld);
  threshFilt->Update();
  ctSystem.forProj->SetInput(1,threshFilt->GetOutput());
  ctSystem.forProj->Update();
  volType projA = ctSystem.forProj->GetOutput();
  threshFilt->ThresholdOutside(1,3);
  threshFilt->Update();
  ctSystem.forProj->SetInput(1,threshFilt->GetOutput());
  ctSystem.forProj->Update();
  volType projB = ctSystem.forProj->GetOutput();
  threshFilt->ThresholdOutside(3,10);
  threshFilt->Update();
  ctSystem.forProj->SetInput(1,threshFilt->GetOutput());
  ctSystem.forProj->Update();
  volType projC = ctSystem.forProj->GetOutput();
  maskType::Pointer maskFilt = maskType::New();
  maskFilt->SetLowerThreshold(1);
  maskFilt->SetUpperThreshold(3);
  maskFilt->SetInsideValue(1);
  maskFilt->SetOutsideValue(0);
  maskFilt->SetInput(param.volOld);
  maskFilt->Update();
  ctSystem.forProj->SetInput(1,maskFilt->GetOutput());
  ctSystem.forProj->Update();
  volType constB = ctSystem.forProj->GetOutput();
  maskFilt->SetLowerThreshold(3);
  maskFilt->SetUpperThreshold(10);
  ctSystem.forProj->SetInput(1,maskFilt->GetOutput());
  ctSystem.forProj->Update();
  volType constC = ctSystem.forProj->GetOutput();
  
  // Assemble the estimate and derivative terms
  expType::Pointer expFilt = expType::New();
  expFilt->SetInput(calc_poly_projection(projA,projB,projC,constB,constC,param,0));
  expFilt->Update();
  
  std::cout << "Generated the forward projections OK!" << std::endl;
  projA->DisconnectPipeline();
  projB->DisconnectPipeline();
  projC->DisconnectPipeline();
  constB->DisconnectPipeline();
  constC->DisconnectPipeline();
  multiplyType::Pointer multFilt = multiplyType::New();
  multFilt->SetInput1(expFilt->GetOutput());
  multFilt->SetConstant2(param.spectrum[0]);
  std::cout << "Spectrum: " << param.spectrum[0] << std::endl;
  multFilt->Update();
  
  volType tmpProj = multFilt->GetOutput();
  tmpProj->DisconnectPipeline();
  addType::Pointer addSpecProj = addType::New();
  
  // First derivative value
  multFilt->SetInput1(tmpProj);
  //multFilt->ReleaseDataFlagOn();
  multFilt->SetConstant2(-param.knee[0]);
  multFilt->Update();
  volType addDerivA = multFilt->GetOutput();
  addDerivA->DisconnectPipeline();
  // Second derivative value
  multFilt->SetConstant2(-param.knee[1]);
  multFilt->Update();
  volType addDerivB = multFilt->GetOutput();
  addDerivB->DisconnectPipeline();
  // Third derivative value
  multFilt->SetConstant2(-param.knee[3]);
  multFilt->Update();
  volType addDerivC = multFilt->GetOutput();
  addDerivC->DisconnectPipeline();
  std::cout << "Initialised all variables OK!" << std::endl;
  
  volType tmpProjIter;
  for (int k = 1; k<2; k++)
  {
  	expFilt->SetInput(calc_poly_projection(projA,projB,projC,constB,constC,param,k));
  	expFilt->Update();
  	multFilt->SetInput1(expFilt->GetOutput());
  	multFilt->SetConstant2(param.spectrum[k]);
  	multFilt->Update();
  	
  	tmpProjIter = multFilt->GetOutput();
  	tmpProjIter->DisconnectPipeline();
  	addSpecProj->SetInput1(tmpProj);
  	addSpecProj->SetInput2(tmpProjIter);
  	addSpecProj->Update();
  	tmpProj = addSpecProj->GetOutput();
  	tmpProj->DisconnectPipeline();
  	// First derivative value
  	multFilt->SetInput1(tmpProjIter);
  	multFilt->SetConstant2(-param.knee[0]);
  	multFilt->Update();
  	addSpecProj->SetInput1(multFilt->GetOutput());
  	addSpecProj->SetInput2(addDerivA);
  	addSpecProj->Update();
  	addDerivA = addSpecProj->GetOutput();
  	addDerivA->DisconnectPipeline();
  	// Second derivative value
  	multFilt->SetConstant2(-param.knee[1]);
  	multFilt->Update();
  	addSpecProj->SetInput1(multFilt->GetOutput());
  	addSpecProj->SetInput2(addDerivA);
  	addSpecProj->Update();
  	addDerivB = addSpecProj->GetOutput();
  	addDerivB->DisconnectPipeline();
  	// Third derivative value
  	multFilt->SetConstant2(-param.knee[3]);
  	multFilt->Update();
  	addSpecProj->SetInput1(multFilt->GetOutput());
  	addSpecProj->SetInput2(addDerivA);
  	addSpecProj->Update();
  	addDerivC = addSpecProj->GetOutput();
  	addDerivC->DisconnectPipeline();
  }
  
  std::cout << "Finished all the precorrection... on to derivative!" << std::endl;
  	
  // Calculate the derivative factors
  
  divType::Pointer divFilt = divType::New();
  divFilt->SetInput1(param.y);
  divFilt->SetInput2(tmpProj);
  divFilt->Update();
  subtractType::Pointer subFilt = subtractType::New();
  subFilt->SetInput1(divFilt->GetOutput());
  subFilt->SetConstant2(1);
	subFilt->Update();
  
	multiplyType::Pointer multOutFilt = multiplyType::New();
  multOutFilt->SetInput1(subFilt->GetOutput());
  multOutFilt->SetInput2(addDerivA);
  multOutFilt->Update();
  
  std::cout << "Calculated the forward derivative OK!" << std::endl;
  ctSystem.backProj->SetInput(1,multOutFilt->GetOutput());
  ctSystem.backProj->Update();
  return ctSystem.backProj->GetOutput();
  /* Linearised gradient calculation
  ctSystem.forProj->SetInput(1,param.volOld);
  ctSystem.forProj->Update();
  
  subtractType::Pointer subFilter = subtractType::New();
  //subFilter->ReleaseDataFlagOn();
  subFilter->SetInput1(ctSystem.forProj->GetOutput());
  subFilter->SetInput2(param.y);
  subFilter->Update();
  
  ctSystem.backProj->SetInput(1,subFilter->GetOutput());
  ctSystem.backProj->Update();
  //ctSystem.backProj->GetOutput()->UpdateOutputInformation();
  //ctSystem.backProj->GetOutput()->PropagateRequestedRegion();
  return ctSystem.backProj->GetOutput();*/
}

volType calc_poly_projection(volType &projA,volType &projB,volType &projC,
            								 volType &constB,volType &constC,paramType &param,int enInd)
{
	multiplyType::Pointer multFilt = multiplyType::New();
	addType::Pointer addFilt = addType::New();
	addFilt->InPlaceOn();
	multFilt->SetInput1(projA);
	multFilt->SetConstant2(param.knee[0]);
	multFilt->Update();
	addFilt->SetInput2(multFilt->GetOutput());
	multFilt->SetInput1(projB);
	multFilt->SetConstant2(param.knee[1]);
	multFilt->Update();
	addFilt->SetInput1(multFilt->GetOutput());
	addFilt->Update();
	//addFilt->SetInput1(this->GetOutput());
	multFilt->SetInput1(constB);
	multFilt->SetConstant2(param.knee[2]);
	multFilt->Update();
	addFilt->SetInput1(multFilt->GetOutput());
	addFilt->Update();
	multFilt->SetInput1(projC);
	multFilt->SetConstant2(param.knee[3]);
	multFilt->Update();
	addFilt->SetInput1(multFilt->GetOutput());
	addFilt->Update();
	multFilt->SetInput1(constC);
	multFilt->SetConstant2(param.knee[4]);
	multFilt->Update();
	addFilt->SetInput1(multFilt->GetOutput());
	addFilt->Update();
	return addFilt->GetOutput();
}