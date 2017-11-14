#include <iostream>
#include <vector>
#include <itkPasteImageFilter.h>
#include "rtkConstantImageSource.h"
//#include "rtkSubSelectFromListImageFilter.h"

void os_polyquant(paramType &param,ctSystemType ctSystem);
int bit_reversal(int index,int max);
geomType calc_subset_geom(ctSystemType &ctSystem,std::vector<int> indArray);
volType calc_subset_proj(paramType &param,std::vector<int> indArray);