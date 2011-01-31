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
 **************************************************************************** */

/*
Functions that calculate the gravitational potential and its first and second
derivatives for the tesseroid.

@author Leonardo Uieda
@date 27 Jan 2011
*/

/** \todo Possible speed up, use pointers for weights and nodes */


#include "grav_tess.h"


/* Calculates gx caused by a tesseroid. */
double tess_gx(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
               GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, kphi, coslatp, coslatc, sinlatp, sinlatc,
           coslon, rc, kappa, res;
    register int i, j, k;
    
    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;
    
    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));

                l_sqr = rp*rp + rc*rc - 2*rp*rc*(sinlatp*sinlatc +
                                                 coslatp*coslatc*coslon);

                kphi = coslatp*sinlatc - sinlatp*coslatc*coslon;

                kappa = rc*rc*coslatc;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(rc*kphi)/pow(l_sqr, 1.5);
            }
        }
    }

    res *= SI2MGAL*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gy caused by a tesseroid. */
double tess_gy(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
               GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, coslatp, coslatc, sinlatp, sinlatc,
           coslon, sinlon, rc, kappa, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));
                sinlon = sin(d2r*(glq_lon.nodes[k] - lonp));

                l_sqr = rp*rp + rc*rc - 2*rp*rc*(sinlatp*sinlatc +
                                                 coslatp*coslatc*coslon);

                kappa = rc*rc*coslatc;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(rc*coslatc*sinlon)/pow(l_sqr, 1.5);
            }
        }
    }

    res *= SI2MGAL*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gz caused by a tesseroid. */
double tess_gz(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
               GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, coslatp, coslatc, sinlatp, sinlatc,
           coslon, cospsi, rc, kappa, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));

                cospsi = sinlatp*sinlatc + coslatp*coslatc*coslon;

                l_sqr = rp*rp + rc*rc - 2*rp*rc*cospsi;

                kappa = rc*rc*coslatc;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(rc*cospsi - rp)/pow(l_sqr, 1.5);
            }
        }
    }

    res *= SI2MGAL*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gxx caused by a tesseroid. */
double tess_gxx(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
                GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, kphi, coslatp, coslatc, sinlatp, sinlatc,
           coslon, rc, kappa, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));

                l_sqr = rp*rp + rc*rc - 2*rp*rc*(sinlatp*sinlatc +
                                                 coslatp*coslatc*coslon);

                kphi = coslatp*sinlatc - sinlatp*coslatc*coslon;

                kappa = rc*rc*coslatc;
                
                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(3*rc*kphi*rc*kphi - l_sqr)/pow(l_sqr, 2.5);
            }
        }
    }

    res *= SI2EOTVOS*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gxy caused by a tesseroid. */
double tess_gxy(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
                GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, kphi, coslatp, coslatc, sinlatp, sinlatc,
           coslon, sinlon, rc, kappa, deltax, deltay, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));
                sinlon = sin(d2r*(glq_lon.nodes[k] - lonp));

                l_sqr = rp*rp + rc*rc - 2*rp*rc*(sinlatp*sinlatc +
                                                 coslatp*coslatc*coslon);

                kphi = coslatp*sinlatc - sinlatp*coslatc*coslon;

                kappa = rc*rc*coslatc;

                deltax = rc*kphi;

                deltay = rc*coslatc*sinlon;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(3*deltax*deltay)/pow(l_sqr, 2.5);
            }
        }
    }

    res *= SI2EOTVOS*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gxz caused by a tesseroid. */
double tess_gxz(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
                GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, kphi, coslatp, coslatc, sinlatp, sinlatc,
           coslon, cospsi, rc, kappa, deltax, deltaz, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));

                cospsi = sinlatp*sinlatc + coslatp*coslatc*coslon;

                l_sqr = rp*rp + rc*rc - 2*rp*rc*cospsi;

                kphi = coslatp*sinlatc - sinlatp*coslatc*coslon;

                kappa = rc*rc*coslatc;

                deltax = rc*kphi;

                deltaz = rc*cospsi - rp;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(3*deltax*deltaz)/pow(l_sqr, 2.5);
            }
        }
    }

    res *= SI2EOTVOS*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gyy caused by a tesseroid. */
double tess_gyy(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
                GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, coslatp, coslatc, sinlatp, sinlatc,
           coslon, sinlon, rc, kappa, deltay, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));
                sinlon = sin(d2r*(glq_lon.nodes[k] - lonp));

                l_sqr = rp*rp + rc*rc - 2*rp*rc*(sinlatp*sinlatc +
                                                 coslatp*coslatc*coslon);

                kappa = rc*rc*coslatc;

                deltay = rc*coslatc*sinlon;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(3*deltay*deltay - l_sqr)/pow(l_sqr, 2.5);
            }
        }
    }

    res *= SI2EOTVOS*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gyz caused by a tesseroid. */
double tess_gyz(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
                GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, coslatp, coslatc, sinlatp, sinlatc,
           coslon, sinlon, cospsi, rc, kappa, deltay, deltaz, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));
                sinlon = sin(d2r*(glq_lon.nodes[k] - lonp));

                cospsi = sinlatp*sinlatc + coslatp*coslatc*coslon;

                l_sqr = rp*rp + rc*rc - 2*rp*rc*cospsi;

                kappa = rc*rc*coslatc;

                deltay = rc*coslatc*sinlon;

                deltaz = rc*cospsi - rp;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(3*deltay*deltaz)/pow(l_sqr, 2.5);
            }
        }
    }

    res *= SI2EOTVOS*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}


/* Calculates gzz caused by a tesseroid. */
double tess_gzz(TESSEROID tess, double lonp, double latp, double rp, GLQ glq_lon,
                GLQ glq_lat, GLQ glq_r)
{
    double d2r = PI/180., l_sqr, coslatp, coslatc, sinlatp, sinlatc,
           coslon, cospsi, rc, kappa, deltaz, res;
    register int i, j, k;

    coslatp = cos(d2r*latp);
    sinlatp = sin(d2r*latp);

    res = 0;

    for(k = 0; k < glq_lon.order; k++)
    {
        for(j = 0; j < glq_lat.order; j++)
        {
            for(i = 0; i < glq_r.order; i++)
            {
                rc = glq_r.nodes[i];
                sinlatc = sin(d2r*glq_lat.nodes[j]);
                coslatc = cos(d2r*glq_lat.nodes[j]);
                coslon = cos(d2r*(lonp - glq_lon.nodes[k]));

                cospsi = sinlatp*sinlatc + coslatp*coslatc*coslon;

                l_sqr = rp*rp + rc*rc - 2*rp*rc*cospsi;

                kappa = rc*rc*coslatc;

                deltaz = rc*cospsi - rp;

                res += glq_lon.weights[k]*glq_lat.weights[j]*glq_r.weights[i]*
                       kappa*(3*deltaz*deltaz - l_sqr)/pow(l_sqr, 2.5);
            }
        }
    }

    res *= SI2EOTVOS*G*tess.density*d2r*(tess.e - tess.w)*d2r*(tess.n - tess.s)*
           (tess.r2 - tess.r1)*0.125;

    return res;
}