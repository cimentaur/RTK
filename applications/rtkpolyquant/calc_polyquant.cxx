#include "rtkpolyquant.h"
#include "calc_polyquant.h"

volType grad_polyquant(paramType &param,ctSystemType &ctSystem)
{
  // TODO: tidy and generalise -- store all temporary variable in 4D array
  // Calculate the forward projections
  thresholdType::Pointer threshFilt = thresholdType::New();
  threshFilt->SetOutsideValue(0);
  threshFilt->ThresholdOutside(0,param.facKnee[0]);
  threshFilt->SetInput(param.recon);
  threshFilt->Update();
  ctSystem.forProj->SetInput(1,threshFilt->GetOutput());
  ctSystem.forProj->Update();
  volType projA = ctSystem.forProj->GetOutput();
  projA->DisconnectPipeline();
  threshFilt->ThresholdOutside(param.facKnee[0],param.facKnee[1]);
  threshFilt->Update();
  ctSystem.forProj->SetInput(1,threshFilt->GetOutput());
  ctSystem.forProj->Update();
  volType projB = ctSystem.forProj->GetOutput();
  projB->DisconnectPipeline();
  threshFilt->ThresholdOutside(param.facKnee[1],10);
  threshFilt->Update();
  ctSystem.forProj->SetInput(1,threshFilt->GetOutput());
  ctSystem.forProj->Update();
  volType projC = ctSystem.forProj->GetOutput();
  maskType::Pointer maskFilt = maskType::New();
  projC->DisconnectPipeline();
  maskFilt->SetLowerThreshold(param.facKnee[0]);
  maskFilt->SetUpperThreshold(param.facKnee[1]);
  maskFilt->SetInsideValue(1);
  maskFilt->SetOutsideValue(0);
  maskFilt->SetInput(param.recon);
  maskFilt->Update();
  ctSystem.forProj->SetInput(1,maskFilt->GetOutput());
  ctSystem.forProj->Update();
  volType constB = ctSystem.forProj->GetOutput();
  constB->DisconnectPipeline();
  maskFilt->SetLowerThreshold(param.facKnee[1]);
  maskFilt->SetUpperThreshold(10);
  ctSystem.forProj->SetInput(1,maskFilt->GetOutput());
  ctSystem.forProj->Update();
  volType constC = ctSystem.forProj->GetOutput();
  constC->DisconnectPipeline();
  
  // Assemble the estimate and derivative terms
  expType::Pointer expFilt = expType::New();
  expFilt->SetInput(calc_poly_projection(projA,projB,projC,constB,constC,param,0));
  expFilt->Update();

  multiplyType::Pointer multFilt = multiplyType::New();
  multFilt->SetInput1(expFilt->GetOutput());
  multFilt->SetConstant2(param.spectrum[0]);
  multFilt->Update();
  
  volType tmpProj = multFilt->GetOutput();
  tmpProj->DisconnectPipeline();
  addType::Pointer addSpecProj = addType::New();
  
  // First derivative value
  multFilt->SetInput1(tmpProj);
  multFilt->SetConstant2(param.knee[0][0]);
  multFilt->Update();
  volType addDerivA = multFilt->GetOutput();
  addDerivA->DisconnectPipeline();
  // Second derivative value
  multFilt->SetConstant2(param.knee[0][1]);
  multFilt->Update();
  volType addDerivB = multFilt->GetOutput();
  addDerivB->DisconnectPipeline();
  // Third derivative value
  multFilt->SetConstant2(param.knee[0][3]);
  multFilt->Update();
  volType addDerivC = multFilt->GetOutput();
  addDerivC->DisconnectPipeline();
  
  volType tmpProjIter;
  for (int k = 1; k<param.spectrum.size(); k++)
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
  	multFilt->SetConstant2(param.knee[k][0]);
  	multFilt->Update();
  	addSpecProj->SetInput1(multFilt->GetOutput());
  	addSpecProj->SetInput2(addDerivA);
  	addSpecProj->Update();
  	addDerivA = addSpecProj->GetOutput();
  	addDerivA->DisconnectPipeline();
  	// Second derivative value
  	multFilt->SetConstant2(param.knee[k][1]);
  	multFilt->Update();
  	addSpecProj->SetInput1(multFilt->GetOutput());
  	addSpecProj->SetInput2(addDerivB);
  	addSpecProj->Update();
  	addDerivB = addSpecProj->GetOutput();
  	addDerivB->DisconnectPipeline();
  	// Third derivative value
  	multFilt->SetConstant2(param.knee[k][3]);
  	multFilt->Update();
  	addSpecProj->SetInput1(multFilt->GetOutput());
  	addSpecProj->SetInput2(addDerivC);
  	addSpecProj->Update();
  	addDerivC = addSpecProj->GetOutput();
  	addDerivC->DisconnectPipeline();
  }
  
  //std::cout << "Finished all the precorrection... on to derivative!" << std::endl;
  	
  // Calculate the derivative factors
  
  divType::Pointer divFilt = divType::New();
  divFilt->SetInput1(param.y);
  divFilt->SetInput2(tmpProj);
  divFilt->Update();
  subtractType::Pointer subFilt = subtractType::New();
  subFilt->SetInput1(divFilt->GetOutput());
  subFilt->SetConstant2(1);
	subFilt->Update();
  
  OutputImageType::RegionType largestRegion = param.recon->GetLargestPossibleRegion();
  // First output term
	multiplyType::Pointer multOutFilt = multiplyType::New();
  multOutFilt->SetInput1(subFilt->GetOutput());
  multOutFilt->SetInput2(addDerivA);
  multOutFilt->Update();
  //ddFilt->SetOffsets(319.614,319.614);
  //ddFilt->SetPadOnTruncatedSide(false);
	ddType::Pointer ddFilt = ddType::New();
  ddFilt->SetInput(multOutFilt->GetOutput());
  ddFilt->SetGeometry(ctSystem.geom);
  ctSystem.backProj->SetInput(1,ddFilt->GetOutput());
  ctSystem.backProj->Update();
  maskingType::Pointer maskOutFilt = maskingType::New();
  maskFilt->SetLowerThreshold(0);
  maskFilt->SetUpperThreshold(param.facKnee[0]);
  maskFilt->Update();
  maskOutFilt->SetInput(ctSystem.backProj->GetOutput());
  maskOutFilt->SetMaskImage(maskFilt->GetOutput());
  maskOutFilt->Update();
  volType tmpOut = maskOutFilt->GetOutput();
  tmpOut->DisconnectPipeline();
  
  // Second output term
  multOutFilt->SetInput2(addDerivB);
  multOutFilt->Update();
  ddFilt->SetInput(multOutFilt->GetOutput());
  ctSystem.backProj->SetInput(1,ddFilt->GetOutput());
  ctSystem.backProj->Update();
  maskFilt->SetLowerThreshold(param.facKnee[0]);
  maskFilt->SetUpperThreshold(param.facKnee[1]);
  maskFilt->Update();
  maskOutFilt->SetInput(ctSystem.backProj->GetOutput());
  maskOutFilt->SetMaskImage(maskFilt->GetOutput());
  maskOutFilt->Update();
  addType::Pointer addOutFilt = addType::New();
  addOutFilt->SetInput1(tmpOut);
  addOutFilt->SetInput2(maskOutFilt->GetOutput());
  addOutFilt->Update();
  tmpOut = addOutFilt->GetOutput();
  tmpOut->DisconnectPipeline();
  
  // Third output term
  multOutFilt->SetInput2(addDerivC);
  multOutFilt->Update();
  ddFilt->SetInput(multOutFilt->GetOutput());
  ctSystem.backProj->SetInput(1,ddFilt->GetOutput());
  ctSystem.backProj->Update();
  maskFilt->SetLowerThreshold(param.facKnee[1]);
  maskFilt->SetUpperThreshold(10);
  maskFilt->Update();
  maskOutFilt->SetInput(ctSystem.backProj->GetOutput());
  maskOutFilt->SetMaskImage(maskOutFilt->GetOutput());
  maskOutFilt->Update();
  addOutFilt->SetInput1(tmpOut);
  addOutFilt->SetInput2(maskOutFilt->GetOutput());
  addOutFilt->Update();
  
  return addOutFilt->GetOutput();
}

volType calc_poly_projection(volType &projA,volType &projB,volType &projC,
            								 volType &constB,volType &constC,paramType &param,int enInd)
{
	multiplyType::Pointer multFilt = multiplyType::New();
	addType::Pointer addFilt = addType::New();
	multFilt->SetInput1(projA);
	multFilt->SetConstant2(-param.knee[enInd][0]);
	multFilt->Update();
	volType tmpProj = multFilt->GetOutput();
	tmpProj->DisconnectPipeline();
	multFilt->SetInput1(projB);
	multFilt->SetConstant2(-param.knee[enInd][1]);
	multFilt->Update();
	addFilt->SetInput1(tmpProj);
	addFilt->SetInput2(multFilt->GetOutput());
	addFilt->Update();
	tmpProj = addFilt->GetOutput();
	tmpProj->DisconnectPipeline();
	multFilt->SetInput1(constB);
	multFilt->SetConstant2(-param.knee[enInd][2]);
	multFilt->Update();
	addFilt->SetInput1(tmpProj);
	addFilt->SetInput2(multFilt->GetOutput());
	addFilt->Update();
	tmpProj = addFilt->GetOutput();
	tmpProj->DisconnectPipeline();
	multFilt->SetInput1(projC);
	multFilt->SetConstant2(-param.knee[enInd][3]);
	multFilt->Update();
	addFilt->SetInput1(tmpProj);
	addFilt->SetInput2(multFilt->GetOutput());
	addFilt->Update();
	tmpProj = addFilt->GetOutput();
	tmpProj->DisconnectPipeline();
	multFilt->SetInput1(constC);
	multFilt->SetConstant2(-param.knee[enInd][4]);
	multFilt->Update();
	addFilt->SetInput1(tmpProj);
	addFilt->SetInput2(multFilt->GetOutput());
	addFilt->Update();
	tmpProj = addFilt->GetOutput();
	tmpProj->DisconnectPipeline();
	return tmpProj;
}
