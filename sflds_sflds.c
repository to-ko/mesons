/******************************************************************************
 * sflds_sflds.c
 *
 * Copyright (C) 2013, 2014 Tomasz Korzec
 *
 * Based on openQCD
 * Copyright (C) 2012 Martin Luescher and Stefan Schaefer
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************
 *
 *
 * functions that should go into
 * modules/sflds/sflds.c
 *
 ******************************************************************************/


#define SLFDS_SFLDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "sflds.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

void random_Z2_sd(int vol,spinor_dble *sd)
{
   double r[12];
   spinor_dble *sm;

   sm=sd+vol;

   for (;sd<sm;sd++)
   {
      ranlxd(r,12);

      (*sd).c1.c1.re=1.0-2.0*(r[ 0]>0.5);
      (*sd).c1.c1.im=0;
      (*sd).c1.c2.re=1.0-2.0*(r[ 1]>0.5);
      (*sd).c1.c2.im=0;
      (*sd).c1.c3.re=1.0-2.0*(r[ 2]>0.5);
      (*sd).c1.c3.im=0;

      (*sd).c2.c1.re=1.0-2.0*(r[ 3]>0.5);
      (*sd).c2.c1.im=0;
      (*sd).c2.c2.re=1.0-2.0*(r[ 4]>0.5);
      (*sd).c2.c2.im=0;
      (*sd).c2.c3.re=1.0-2.0*(r[ 5]>0.5);
      (*sd).c2.c3.im=0;

      (*sd).c3.c1.re=1.0-2.0*(r[ 6]>0.5);
      (*sd).c3.c1.im=0;
      (*sd).c3.c2.re=1.0-2.0*(r[ 7]>0.5);
      (*sd).c3.c2.im=0;
      (*sd).c3.c3.re=1.0-2.0*(r[ 8]>0.5);
      (*sd).c3.c3.im=0;

      (*sd).c4.c1.re=1.0-2.0*(r[ 9]>0.5);
      (*sd).c4.c1.im=0;
      (*sd).c4.c2.re=1.0-2.0*(r[10]>0.5);
      (*sd).c4.c2.im=0;
      (*sd).c4.c3.re=1.0-2.0*(r[11]>0.5);
      (*sd).c4.c3.im=0;
   }
}

void random_U1_sd(int vol,spinor_dble *sd)
{
   double r[12];
   spinor_dble *sm;

   sm=sd+vol;

   for (;sd<sm;sd++)
   {
      ranlxd(r,12);

      (*sd).c1.c1.re=cos(2*M_PI*r[0]);
      (*sd).c1.c1.im=sin(2*M_PI*r[0]);
      (*sd).c1.c2.re=cos(2*M_PI*r[1]);
      (*sd).c1.c2.im=sin(2*M_PI*r[1]);
      (*sd).c1.c3.re=cos(2*M_PI*r[2]);
      (*sd).c1.c3.im=sin(2*M_PI*r[2]);

      (*sd).c2.c1.re=cos(2*M_PI*r[3]);
      (*sd).c2.c1.im=sin(2*M_PI*r[3]);
      (*sd).c2.c2.re=cos(2*M_PI*r[4]);
      (*sd).c2.c2.im=sin(2*M_PI*r[4]);
      (*sd).c2.c3.re=cos(2*M_PI*r[5]);
      (*sd).c2.c3.im=sin(2*M_PI*r[5]);

      (*sd).c3.c1.re=cos(2*M_PI*r[6]);
      (*sd).c3.c1.im=sin(2*M_PI*r[6]);
      (*sd).c3.c2.re=cos(2*M_PI*r[7]);
      (*sd).c3.c2.im=sin(2*M_PI*r[7]);
      (*sd).c3.c3.re=cos(2*M_PI*r[8]);
      (*sd).c3.c3.im=sin(2*M_PI*r[8]);

      (*sd).c4.c1.re=cos(2*M_PI*r[9]);
      (*sd).c4.c1.im=sin(2*M_PI*r[9]);
      (*sd).c4.c2.re=cos(2*M_PI*r[10]);
      (*sd).c4.c2.im=sin(2*M_PI*r[10]);
      (*sd).c4.c3.re=cos(2*M_PI*r[11]);
      (*sd).c4.c3.im=sin(2*M_PI*r[11]);
   }
}

/* assign minus sd to rd */
void assign_msd2sd(int vol,spinor_dble *sd,spinor_dble *rd)
{
   spinor_dble *sm;

   sm=sd+vol;

   for (;sd<sm;sd++)
   {
      (*rd).c1.c1.re=-(*sd).c1.c1.re; (*rd).c1.c1.im=-(*sd).c1.c1.im;
      (*rd).c1.c2.re=-(*sd).c1.c2.re; (*rd).c1.c2.im=-(*sd).c1.c2.im;
      (*rd).c1.c3.re=-(*sd).c1.c3.re; (*rd).c1.c3.im=-(*sd).c1.c3.im;
      (*rd).c2.c1.re=-(*sd).c2.c1.re; (*rd).c2.c1.im=-(*sd).c2.c1.im;
      (*rd).c2.c2.re=-(*sd).c2.c2.re; (*rd).c2.c2.im=-(*sd).c2.c2.im;
      (*rd).c2.c3.re=-(*sd).c2.c3.re; (*rd).c2.c3.im=-(*sd).c2.c3.im;
      (*rd).c3.c1.re=-(*sd).c3.c1.re; (*rd).c3.c1.im=-(*sd).c3.c1.im;
      (*rd).c3.c2.re=-(*sd).c3.c2.re; (*rd).c3.c2.im=-(*sd).c3.c2.im;
      (*rd).c3.c3.re=-(*sd).c3.c3.re; (*rd).c3.c3.im=-(*sd).c3.c3.im;
      (*rd).c4.c1.re=-(*sd).c4.c1.re; (*rd).c4.c1.im=-(*sd).c4.c1.im;
      (*rd).c4.c2.re=-(*sd).c4.c2.re; (*rd).c4.c2.im=-(*sd).c4.c2.im;
      (*rd).c4.c3.re=-(*sd).c4.c3.re; (*rd).c4.c3.im=-(*sd).c4.c3.im;
      rd+=1;
   }
}
