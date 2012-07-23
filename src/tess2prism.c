/*
Convert a tesseroid model into a prism model in spherical coordinates
*/


#include <stdio.h>
#include <time.h>
#include "version.h"
#include "parsers.h"
#include "logger.h"
#include "geometry.h"


/** Print the help message */
void print_help()
{
    printf("Usage: tess2prim TESSFILE [OPTIONS]\n\n");
    printf("Convert a tesseroid model into a rectangular prism model\n");
    printf("(for use with the prism* programs).\n\n");
    printf("The converted prisms have the same mass as the tesseroids.\n\n");
    printf("Along with each prism are given the spherical coordinates of the\n");
    printf("center of the top face of the tesseroid (used as the origin of\n");
    printf("the prisms coordinate system). This is needed to compute the.\n");
    printf("effect of the prisms in spherical coordinates.\n\n");
    printf("If option --flatten is used, the tesseroids are converted by\n");
    printf("approximating 1 degree by 111.11km and no spherical coordinates\n");
    printf("are given. Use this option when you want to calculate in\n");
    printf("Cartesian coordinates.\n\n");
    printf("In both cases, the density of the prism is adjusted so that it\n");
    printf("has the same mass as the tesseroid.\n\n");
    printf("All units either SI or degrees!\n\n");
    printf("Input:\n");
    printf("  If TESSFILE is omited, will read from standard input (stdin)\n");
    printf("  TESSFILE: File containing the tesseroid model\n");
    printf("  * Each tesseroid is specified by the values of its borders\n");
    printf("    and density\n");
    printf("  * The file should contain one tesseroid per line\n");
    printf("  * Each line should have the following column format:\n");
    printf("      West East South North Top Bottom Density\n");
    printf("  * Top and Bottom should be read as 'height to top' and \n");
    printf("    'height to bottom' from the mean Earth radius. Use negative\n");
    printf("    values if bellow the surface, for example when modeling\n");
    printf("    deep structures, and positive if above the surface, for\n");
    printf("    example when modeling topography.\n");
    printf("  * If a line starts with # it will be considered a comment\n");
    printf("    and will be ignored\n\n");
    printf("Output:\n");
    printf("  Printed to standard output (stdout) one prism per line in the\n");
    printf("  format:\n");
    printf("    dx dy dz density lon lat r\n");
    printf("  lon, lat, r are the spherical coordinates of the center of\n");
    printf("  top face of the prism. This is used as the origin of the\n");
    printf("  local coordinate system of the prism.\n");
    printf("  If options --flatten is used, the output format is:\n");
    printf("    x1 x2 y1 y2 z1 z2 density\n");
    printf("  Comments about the provenance of the data are inserted into\n");
    printf("  the top of the output.\n\n");
    printf("Options:\n");
    printf("  --flatten    Convert the tesseroids by approximating 1 degree\n");
    printf("               by 111.11 km (for compatibility with prism*\n");
    printf("               programs).\n");
    printf("  -h           Print instructions.\n");
    printf("  --version    Print version and license information.\n");
    printf("  -v           Enable verbose printing to stderr.\n");
    printf("  -lFILENAME   Print log messages to file FILENAME.\n");
    printf("\nPart of the Tesseroids package.\n");
    printf("Project site: <http://fatiando.org/software/tesseroids>\n");
    printf("Report bugs at: ");
    printf("<http://code.google.com/p/tesseroids/issues/list>\n");
}


/** Main */
int main(int argc, char **argv)
{
    char *progname = "tess2prism";
    TESS2PRISM_ARGS args;
    int rc, line, converted = 0, error_exit = 0, bad_input = 0;
    char buff[10000];
    TESSEROID tess;
    PRISM prism;
    FILE *logfile = NULL, *modelfile = NULL;
    time_t rawtime;
    struct tm * timeinfo;

    log_init(LOG_INFO);
    rc = parse_tess2prism_args(argc, argv, progname, &args, &print_help);
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
    if(!args.verbose)
    {
        log_init(LOG_WARNING);
    }
    if(args.logtofile)
    {
        logfile = fopen(args.logfname, "w");
        if(logfile == NULL)
        {
            log_error("unable to create log file %s", args.logfname);
            log_warning("Terminating due to bad input");
            log_warning("Try '%s -h' for instructions", progname);
            return 1;
        }
        log_tofile(logfile, LOG_INFO);
    }

    /* Print standard verbose */
    log_info("%s (Tesseroids project) %s", progname, tesseroids_version);
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    log_info("(local time) %s", asctime(timeinfo));

    /* If an input file is not given, read from stdin. Else open the file */
    if(rc == 3)
    {
        log_info("Reading tesseroids from stdin");
        modelfile = stdin;
    }
    else
    {
        log_info("Reading tesseroids from file %s", args.inputfname);
        modelfile = fopen(args.inputfname, "r");
        if(modelfile == NULL)
        {
            log_error("failed to open file %s", args.inputfname);
            log_warning("Terminating due to bad input");
            log_warning("Try '%s -h' for instructions", progname);
            if(args.logtofile)
                fclose(logfile);
            return 1;
        }
    }

    /* Print provenance data to stdout */
    printf("# Prisms converted from tesseroid model with %s %s\n", progname,
           tesseroids_version);
    printf("#   local time: %s", asctime(timeinfo));
    printf("#   tesseroids file: %s\n", rc == 3 ? "stdin" : args.inputfname);
    printf("#   conversion type: %s\n",
        args.flatten ? "equal mass|flatten" :
                       "equal mass|spherical coordinates");
    if(args.flatten)
    {
        printf("#   format: x1 x2 y1 y2 z1 z2 density\n");
    }
    else
    {
        printf("#   format: dx dy dz density lon lat r\n");
    }

    /* Read the tesseroids, convert and print to stdout */
    for(line = 1; !feof(modelfile); line++)
    {
        if(fgets(buff, 10000, modelfile) == NULL)
        {
            if(ferror(stdin))
            {
                log_error("problem encountered reading line %d", line);
                error_exit = 1;
                break;
            }
        }
        else
        {
            /* Check for comments and blank lines */
            if(buff[0] == '#' || buff[0] == '\r' || buff[0] == '\n')
            {
                printf("%s", buff);
                continue;
            }
            /* Remove any trailing spaces or newlines */
            strstrip(buff);
            if(gets_tess(buff, &tess))
            {
                log_warning("bad/invalid tesseroid at line %d", line);
                bad_input++;
                continue;
            }
            if(args.flatten)
            {
                tess2prism_flatten(tess, &prism);
                printf("%.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
                       prism.x1, prism.x2, prism.y1, prism.y2, prism.z1,
                       prism.z2, prism.density);
            }
            else
            {
                tess2prism(tess, &prism);
                printf("%.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
                       prism.x2 - prism.x1, prism.y2 - prism.y1,
                       prism.z2 - prism.z1, prism.density,
                       prism.lon, prism.lat, prism.r);
            }
            converted++;
        }
    }
    if(bad_input)
    {
        log_warning("Encountered %d bad input line(s) which were skipped",
                    bad_input);
    }
    if(error_exit)
    {
        log_warning("Terminating due to error in input");
        log_warning("Try '%s -h' for instructions", progname);
    }
    else
    {
        log_info("Converted %d tesseroids", converted);
    }
    /* Clean up */
    fclose(modelfile);
    if(args.logtofile)
        fclose(logfile);
    return 0;
}
