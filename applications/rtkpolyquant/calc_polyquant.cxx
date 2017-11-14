#include "rtkpolyquant.h"
#include "calc_polyquant.h"

volType grad_polyquant(paramType &param,ctSystemType &ctSystem)
{
  // TODO: calculate full Polyquant gradient term
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
  return ctSystem.backProj->GetOutput();
  /*specProb = specRat(1);

  tmp1 = specProb.*exp(-(specData.boneMa(1)*projA));
  tmpA = tmp1*specData.boneMa(1);

  for k = 2:length(specData.spectrum)
  {
    tmp = specProb.*exp(-(specData.boneMa(k)*projA));
    tmp1 = tmp1+tmp;
    tmpA = tmpA+tmp*specData.boneMa(k);
  }
	tmp1 = I0.*tmp1;
	tmpA = I0.*tmpA;

	deriFac = (y./(tmp1+s)-1);
	outA = w(tmpA.*deriFac);

	out = At(outA,ind);*/
}