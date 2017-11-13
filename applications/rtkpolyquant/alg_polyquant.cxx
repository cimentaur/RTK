#include "rtkpolyquant.h"
#include "alg_polyquant.h"
#include <iostream>

void os_polyquant(paramType param,ctSystemType ctSystem)
{
	std::cout << "Number of iterations = " << param.nIter << std::endl;
  /*x0 = in;
  x1 = x0;
  for (int k = 0; k < param.nIter; k++)
  {
  	subSet = bit_reversal(k,param.nIter);
		grad = grad_polyquant(x0,param,ctSystem);
		
		xNew = prox(x1-param.stepSize*grad);
		
		if (param.accelerate == true)
		{
		  // FISTA type acceleration
			t1 = 0.5*(1+sqrt(1+4*t^2));
      x1 = xNew+(t-1)/t1*(xNew-x0);
      x0 = xNew;
      t = t1;
		}
		else
		{
		  x0 = x1;
      x1 = xNew;
		}
		
  }*/
}

int bit_reversal(int index, int max)
{
	// TODO: apply bit reversal strategy
	return index;
}