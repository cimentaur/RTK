#include "other.h"
#include <iostream>

void calc_polyquant(const &in,const &param,const ctSystemType &ctSystem);
{
  // TODO: calculate full Polyquant gradient term
  projA = ctSystem.projFor;

  specProb = specRat(1);

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

	out = At(outA,ind);
}