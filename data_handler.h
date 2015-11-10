/*
  This header file contains the prototypes for
  data_handler.c. 
  See accompanying .c file for more information.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"
#include "galaxy_structs.h"

#ifndef HANDLER_H

galaxy* readData(int* numGals);
cartesianGalaxy* convertToCartesian(galaxy* gals, int numGals);
galaxy* convertFromCartesian(cartesianGalaxy* cartGals, galaxy* origGals,int numGals);
galaxy *trimGalaxyList(galaxy *gals, int *numGals);


#endif
