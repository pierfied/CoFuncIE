/*
  This header file contains the structure used to contain 
  the map, or the histogram of the logarithmic density field.
 */

#ifndef MAPSTRUCT_H
#define MAPSTRUCT_H

/* TODO: Define the structure for the map 
   I wasn't sure if you wanted to do a 1D array or a 3D array.
   The 1D array will be of length numdivs^3, whereas
   the 3D array will have each dimension be of length numdivs.
   Just delete whichever one you don't use.
*/
typedef struct{
  int numdivs; //number of divisions on a side
  double* option1;
  double*** option2;
} map;

#endif
