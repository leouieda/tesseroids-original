/* *****************************************************************************
Copyright 2011 Leonardo Uieda

Tesseroids is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Tesseroids is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tesseroids.  If not, see <http://www.gnu.org/licenses/>.
***************************************************************************** */

/** \file
Set of misc utilities and data structures.

Defines the TESSEROID, SPHERE and PRISM structures.

\todo Make functions that calculate the center of mass

@author Leonardo Uieda
@date 25 Jan 2011
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "logger.h"
#include "utils.h"


/* Convert a tesseroid to a rectangular prism of equal volume. */
void tess2prism(TESSEROID tess, PRISM *prism)
{
    double deg2rad = PI/180., r0, dx, dy;

    r0 = 0.5*(tess.r1 + tess.r2);

    dx = r0*deg2rad*(tess.n - tess.s);
    dy = r0*cos(deg2rad*0.5*(tess.n + tess.s))*deg2rad*(tess.e - tess.w);

    prism->x1 = -0.5*dx;
    prism->x2 = 0.5*dx;
    prism->y1 = -0.5*dy;
    prism->y2 = 0.5*dy;
    /* r1 is not z1 because r1 is the bottom face */
    prism->z1 = MEAN_EARTH_RADIUS - tess.r2;
    prism->z2 = MEAN_EARTH_RADIUS - tess.r1;
    prism->density = tess.density;
}


/* Convert a tesseroid to a sphere of equal volume. */
void tess2sphere(TESSEROID tess, SPHERE *sphere)
{
    sphere->density = tess.density;
    /** \todo Put sphere in center of mass, not geometrical center */
    sphere->lonc = 0.5*(tess.e + tess.w);
    sphere->latc = 0.5*(tess.n + tess.s);
    sphere->rc = 0.5*(tess.r1 + tess.r2);
    sphere->r = cbrt(3*tess_volume(tess)/(4.*PI));
}


/* Convert a rectangular prism into a sphere of equal volume. */
void prism2sphere(PRISM prism, double lonc, double latc, double rc,
                  SPHERE *sphere)
{
    sphere->density = prism.density;
    sphere->lonc = lonc;
    sphere->latc = latc;
    sphere->rc = rc;
    sphere->r = cbrt(3*prism_volume(prism)/(4.*PI));
}


/* Calculate the volume of a tesseroid */
double tess_volume(TESSEROID tess)
{
    double d2r = PI/180., vol;

    vol = d2r*(tess.e - tess.w)*(pow(tess.r2, 3) - pow(tess.r1, 3))*
          (sin(d2r*tess.n) - sin(d2r*tess.s))/3.;

    return vol;
}


/* Calculate the volume of a sphere */
double sphere_volume(SPHERE sphere)
{
    return 4.*PI*pow(sphere.r, 3)/3.;
}


/* Calculate the volume of a prism */
double prism_volume(PRISM prism)
{
    return (prism.x2 - prism.x1)*(prism.y2 - prism.y1)*(prism.z2 - prism.z1);
}


/* Strip trailing spaces and newlines from the end of a string */
void strstrip(char *str)
{
    int i;
    for(i = strlen(str) - 1; i >= 0; i--)
    {
        if(str[i] != ' ' && str[i] != '\n' && str[i] != '\r' && str[i] != '\0')
            break;
    }
    str[i + 1] = '\0';
}


/* Read tesseroids from an open file and store them in an array */
TESSEROID * read_tess_model(FILE *modelfile, int *size)
{
    TESSEROID *model;
    int buffsize = 100;

    /* Start with a single buffer allocation and expand later if necessary */
    model = (TESSEROID *)malloc(buffsize*sizeof(TESSEROID));
    if(model == NULL)
    {
        log_error("problem allocating initial memory to load tesseroid model");
        return NULL;
    }

    int nread, nchars, line, badinput = 0;
    char sbuff[10000];
    double w, e, s, n, top, bot, dens;
    TESSEROID *tmp;

    *size = 0;
    
    for(line = 1; !feof(modelfile); line++)
    {
        if(fgets(sbuff, 10000, modelfile) != NULL)
        {
            /* Check for comments and blank lines */
            if(sbuff[0] == '#' || sbuff[0] == '\r' || sbuff[0] == '\n')
            {
                continue;
            }

            /* Remove any trailing spaces or newlines */
            strstrip(sbuff);

            nread = sscanf(sbuff, "%lf %lf %lf %lf %lf %lf %lf%n", &w, &e, &s,
                           &n, &top, &bot, &dens, &nchars);

            if(nread != 7 || sbuff[nchars] != '\0')
            {
                log_warning("bad/invalid tesseroid at line %d", line);
                badinput = 1;
            }

            if(*size == buffsize)
            {
                buffsize += buffsize;
                tmp = (TESSEROID *)realloc(model, buffsize*sizeof(TESSEROID));
                if(tmp == NULL)
                {
                    /* Need to free because realloc leaves unchanged in case of
                       error */
                    free(model);
                    log_error("problem expanding memory for tesseroid model");
                    return NULL;
                }
                model = tmp;
            }

            model[*size].w = w;
            model[*size].e = e;
            model[*size].s = s;
            model[*size].n = n;
            model[*size].r1 = MEAN_EARTH_RADIUS - bot;
            model[*size].r2 = MEAN_EARTH_RADIUS - top;
            model[*size].density = dens;

            *size += 1;
        }
    }

    if(badinput)
    {
        free(model);
        return NULL;
    }

    /* Adjust the size of the model */
    tmp = (TESSEROID *)realloc(model, (*size)*sizeof(TESSEROID));
    if(tmp == NULL)
    {
        /* Need to free because realloc leaves unchanged in case of
            error */
        free(model);
        log_error("problem freeing extra memory for tesseroid model");
        return NULL;
    }
    model = tmp;
    
    return model;
}