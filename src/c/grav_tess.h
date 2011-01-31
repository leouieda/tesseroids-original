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

/** \file
Functions that calculate the gravitational potential and its first and second
derivatives for the tesseroid.

The gravity gradients can be calculated using the general formula:

\f[
g_{\alpha\beta}(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    I_{\alpha\beta}\ d r' d \phi' d \lambda'
    \ \ \alpha,\beta \in \{1,2,3\}
\f]

\f[
I_{\alpha\beta} = \left(\frac{3\Delta x_i \Delta x_j}{\ell^5} -
    \frac{\delta_{ij}}{\ell^3} \right)\kappa\
\f]

and solved using the Gauss-Legendre Quadrature rule:

\f[
g_{\alpha\beta}(r_p,\phi_p,\lambda_p) \approx G \rho \frac{(\lambda_2 - \lambda_1)
    (\phi_2 - \phi_1)(r_2 - r_1)}{8} \displaystyle\sum_{k=0}^{N^{\lambda} - 1}
    \displaystyle\sum_{j=0}^{N^{\phi} - 1} \displaystyle\sum_{i=0}^{N^r - 1}
    W^r_i W^{\phi}_j W^{\lambda}_k
    I_{\alpha\beta}({r'}_i, {\phi'}_j, {\lambda'}_k )\kappa\ \ \alpha,\beta \in \{1,2,3\}
\f]

where \f$ \rho \f$ is density, the subscripts 1, 2, and 3 should be
interpreted as the x, y, and z axis and

\f{eqnarray*}{
\Delta x_1 &=& r' K_{\phi} \\
\Delta x_2 &=& r' \cos \phi' \sin(\lambda' - \lambda_p) \\
\Delta x_3 &=& r' \cos \psi - r_p\\
\ell &=& \sqrt{r'^2 + r_p^2 - 2 r' r_p \cos \psi} \\
\cos\psi &=& \sin\phi_p\sin\phi' + \cos\phi_p\cos\phi'
             \cos(\lambda' - \lambda_p) \\
K_{\phi} &=& \cos\phi_p\sin\phi' - \sin\phi_p\cos\phi'
             \cos(\lambda' - \lambda_p)\\
\kappa &=& {r'}^2 \cos \phi'
\f}

\f$ \phi \f$ is latitude, \f$ \lambda \f$ is longitude, \f$ r \f$ is radius. The
subscript \f$ p \f$ is for the computation point.

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>. So it would be normal for a tesseroid of
positive density to have negative gz.

<b>Example</b>:

To calculate the gzz component due to a tesseroid on a regular grid.

\verbatim
#include <stdio.h>
#include "glq.h"
#include "constants.h"
#include "grav_tess.h"

int main()
{
    TESSEROID tess = {1000, 44, 46, -1, 1, MEAN_EARTH_RADIUS - 100000,
                      MEAN_EARTH_RADIUS};
    GLQ *glqlon, *glqlat, *glqr;
    double lon, lat, r = MEAN_EARTH_RADIUS + 1500000, res;
    int order = 8;
    
    glqlon = glq_new(order, tess.w, tess.e);
    glqlat = glq_new(order, tess.s, tess.n);
    glqr = glq_new(order, tess.r1, tess.r2);

    for(lat = 20; lat <= 70; lat += 0.5)
    {
        for(lon = -25; lon <= 25; lon += 0.5)
        {
            res = tess_gzz(tess, lon, lat, r, *glqlon, *glqlat, *glqr);
            printf("%g %g %g\n", lon, lat, res);
        }
    }

    glq_free(glqlon);
    glq_free(glqlat);
    glq_free(glqr);

    return 0;
}
\endverbatim 

@author Leonardo Uieda
@date 27 Jan 2011
*/

#ifndef _GRAV_TESS_H_
#define _GRAV_TESS_H_


#include <math.h>
#include "utils.h"
#include "glq.h"
#include "constants.h"


/** Calculates gx caused by a tesseroid.

\f[
g_x(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{r'K_{\phi}}{\ell^3}\kappa \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in mGal!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gx(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);
                    
/** Calculates gy caused by a tesseroid.

\f[
g_y(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{r'\cos\phi'\sin(\lambda'-\lambda)}{\ell^3}\kappa
    \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in mGal!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gy(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

/** Calculates gz caused by a tesseroid.

\f[
g_z(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{r'\cos\psi - r_p}{\ell^3}\kappa \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in mGal!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gz(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

/** Calculates gxx caused by a tesseroid.

\f[
g_{xx}(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{3(r' K_{\phi})^2 - \ell^2}{\ell^5}\kappa \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in Eotvos!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gxx(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

/** Calculates gxy caused by a tesseroid.

\f[
g_{xy}(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{3{r'}^2 K_{\phi}\cos\phi'\sin(\lambda' - \lambda_p)}{\ell^5}
    \kappa \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in Eotvos!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gxy(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

/** Calculates gxz caused by a tesseroid.

\f[
g_{xz}(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{3 r' K_{\phi}(r' \cos\psi - r_p)}{\ell^5}\kappa
    \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in Eotvos!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gxz(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

/** Calculates gyy caused by a tesseroid.

\f[
g_{yy}(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{3(r'\cos\phi'\sin(\lambda' - \lambda_p))^2 - \ell^2}{\ell^5}
    \kappa \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in Eotvos!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gyy(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

/** Calculates gyz caused by a tesseroid.

\f[
g_{yz}(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{3 r' \cos\phi' \sin(\lambda' - \lambda_p)(r'\cos\psi - r_p)}{\ell^5}
    \kappa \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in Eotvos!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gyz(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

/** Calculates gzz caused by a tesseroid.

\f[
g_{zz}(r_p,\phi_p,\lambda_p) = G \rho \displaystyle\int_{\lambda_1}^{\lambda_2}
    \displaystyle\int_{\phi_1}^{\phi_2} \displaystyle\int_{r_1}^{r_2}
    \frac{3(r'\cos\psi-r_p)^2 - \ell^2}{\ell^5}\kappa \ d r' d \phi' d \lambda'
\f]

The derivatives of the potential are made with respect to the local coordinate
system <b>x->North, y->East, z->out</b>

<b>Input values in SI units and <b>degrees</b> and returns values in Eotvos!</b>

Use function glq_new() to create the GLQ parameters required. The integration
limits should be set to:
    - glq_lon: lower = tess.w and upper = tess.e  (in degrees)
    - glq_lat: lower = tess.s and upper = tess.n  (in degrees)
    - glq_r: lower = tess.r1 and upper = tess.r2

@param tess data structure describing the tesseroid
@param lonp longitude of the computation point P
@param latp latitude of the computation point P
@param rp radial coordinate of the computation point P
@param glq_lon GLQ structure with the nodes, weights and integration limits set
    for the longitudinal integration
@param glq_lat GLQ structure with the nodes, weights and integration limits set
    for the latitudinal integration
@param glq_r GLQ structure with the nodes, weights and integration limits set
    for the radial integration

@return field calculated at P
*/
extern double tess_gzz(TESSEROID tess, double lonp, double latp, double rp,
                      GLQ glq_lon, GLQ glq_lat, GLQ glq_r);

#endif
