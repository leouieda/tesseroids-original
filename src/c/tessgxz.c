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
Program to calculate gxz of a tesseroid model on a set of points.

@author Leonardo Uieda
@date 02 Feb 2011
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "logger.h"
#include "version.h"
#include "grav_tess.h"
#include "glq.h"
#include "constants.h"
#include "utils.h"
#include "cmd.h"


/** Main */
int main(int argc, char **argv)
{
    log_init(LOG_INFO);
    char progname[] = "tessgxz";
    TESSG_ARGS args;

    int rc = parse_tessg_args(argc, argv, progname, &args);

    if(rc == 2)
    {
        return 0;
    }
    if(rc == 1)
    {
        log_warning("Terminating due to bad input");
        log_warning("Try '%s -h' for instructions", progname);

        return 1;
    }

    /* Set the appropriate logging level and log to file if necessary */
    if(!args.verbose) { log_init(LOG_WARNING); }

    FILE *logfile;
    if(args.logtofile)
    {
        logfile = fopen(args.logfname, "w");
        if(logfile == NULL)
        {
            log_error("unable to create log file %s\n", args.logfname);
            log_warning("Terminating due to bad input");
            log_warning("Try '%s -h' for instructions", progname);
            return 1;
        }
        log_tofile(logfile, LOG_INFO);
    }

    /* Print standard verbose */
    log_info("%s (Tesseroids project) %s", progname, tesseroids_version);
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    log_info("(local time) %s", asctime(timeinfo));

    /* Make the necessary GLQ structures */
    log_info("Using GLQ orders: %d lon / %d lat / %d r", args.lon_order,
             args.lat_order, args.r_order);
    GLQ *glq_lon, *glq_lat, *glq_r;
    glq_lon = glq_new(args.lon_order, -1, 1);
    glq_lat = glq_new(args.lat_order, -1, 1);
    glq_r = glq_new(args.r_order, -1, 1);
    if(glq_lon == NULL || glq_lat == NULL || glq_r == NULL)
    {
        log_error("failed to create required GLQ structures\n");
        log_warning("Terminating due to bad input");
        log_warning("Try '%s -h' for instructions", progname);
        return 1;
    }

    /* Read the tesseroid model file */
    log_info("Reading tesseroid model from file %s", args.modelfname);
    FILE *modelfile = fopen(args.modelfname, "r");
    if(modelfile == NULL)
    {
        log_error("failed to open model file %s\n", args.modelfname);
        log_warning("Terminating due to bad input");
        log_warning("Try '%s -h' for instructions", progname);
        return 1;
    }
    TESSEROID *model;
    int modelsize;
    model = read_tess_model(modelfile, &modelsize);
    fclose(modelfile);
    if(modelsize == 0)
    {
        log_error("failed to read model from file %s\n", args.modelfname);
        log_warning("Terminating due to bad input");
        log_warning("Try '%s -h' for instructions", progname);
        return 1;
    }
    log_info("Total of %d tesseroid(s) read", modelsize);

    /* Read each computation point from stdin and calculate */
    log_info("Calculating %s component (this may take a while)...", progname+4);
    clock_t tstart = clock();
    int line, points = 0;
    char buff[10000];
    double lon, lat, height, res;

    for(line = 1; !feof(stdin); line++)
    {
        if(fgets(buff, 10000, stdin) != NULL)
        {
            /* Check for comments and blank lines */
            if(buff[0] == '#' || buff[0] == '\r' || buff[0] == '\n')
            {
                printf("%s", buff);
                continue;
            }
            if(sscanf(buff, "%lf %lf %lf", &lon, &lat, &height) != 3)
            {
                log_warning("bad/invalid computation point at line %d", line);
                log_warning("skipping this line and continuing");
                continue;
            }

            /* Need to remove \n and \r from end of buff first to print the
               result in the end */
            strstrip(buff);

            res = calc_tess_model(model, modelsize, lon, lat,
                                  height + MEAN_EARTH_RADIUS, glq_lon, glq_lat,
                                  glq_r, tess_gxz);

            printf("%s %g\n", buff, res);
            points++;
        }
    }

    log_info("Calculated on %d points in %.5g seconds", points,
             (double)(clock() - tstart)/CLOCKS_PER_SEC);

    /* Clean up */
    if(args.logtofile)
    {
        log_info("Closing log file");
        fclose(logfile);
    }
    free(model);
    glq_free(glq_lon);
    glq_free(glq_lat);
    glq_free(glq_r);

    log_info("Done");

    return 0;
}