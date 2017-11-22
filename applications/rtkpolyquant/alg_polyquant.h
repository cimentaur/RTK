#include <iostream>
#include <vector>
#include <itkPasteImageFilter.h>
#include "rtkConstantImageSource.h"
#include <unistd.h>
#include <iomanip>

void os_polyquant(paramType &param,ctSystemType ctSystem);
int bit_reversal(int index,int max);
geomType calc_subset_geom(ctSystemType &ctSystem,std::vector<int> indArray);
volType calc_subset_proj(paramType &param,std::vector<int> indArray);
volType proj_nn(volType in,float upValue);
void print_stats(int iter,paramType &param,itk::TimeProbe &subSetProbe);
