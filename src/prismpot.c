/*
Program to calculate potential of a rectangular prism model on a set of points.
*/


#include "grav_prism.h"
#include "prismg_main.h"


/** Main */
int main(int argc, char **argv)
{
    return run_prismg_main(argc, argv, "prismpot", &prism_pot);
}
