#include "bit_reversal.h"

int calc_bit_reversal(int k,int nProj)
{
  std::vector<int> factorArray = prime_factors(nProj);
  std::vector<int> unit;
  unit.push_back(1);
  int len = factorArray.size();
  if (factorArray.size() > 1)
  {
    for (int i = 1; i<len; i++)
    {
      int tmpUnit = 1;
      for (int j = 0; j<i; j++)
      {
        tmpUnit *= factorArray[j];
      }
      unit.push_back(tmpUnit);
    }
  }
  int remain = k;
  std::vector<int> vec;
  for (int i = 0; i<len; i++)
  {
    vec.push_back(0);
  }
  int ind;

  while (remain>0)
  {
    ind = 0;
    while (unit[ind]*factorArray[ind]-1 < remain)
    {
      ind += 1;
    }
    vec[ind] += 1;
    remain = remain-unit[ind];
  }
  int out = vec.back();
  ind = 1;
  while (ind < len)
  {
    int tmpVal = 1;
    for (int i = len-ind; i<len; i++)
    {
      tmpVal *= factorArray[i];
    }
    out += tmpVal*vec[len-ind-1];
    ind += 1;
  }
  return out;
}

std::vector<int> prime_factors(int x)
{
  std::vector<int> out;
  while (x%2 == 0)
  {
    out.push_back(2);
    x = x/2;
  }
  // n must be odd at this point.  So we can skip 
  // one element (Note i = i +2)
  for (int i = 3; i <= sqrt(x); i = i+2)
  {
      // While i divides n, print i and divide n
      while (x%i == 0)
      {
          out.push_back(i);
          x = x/i;
      }
  }
  // This condition is to handle the case when n 
  // is a prime number greater than 2
  if (x > 2)
  {
    out.push_back(x);
  }
  return out;
}
