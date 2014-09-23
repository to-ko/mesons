/* linalg_salg_dble.c
 *
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
 *******************************************************************************
 *
 *
 * functions that should go into
 * modules/linalg/salg_dble.c
 *
 * The gamma conventions are
 *
 * g0 = [ 0     0    -1     0 ]
 *      [ 0     0     0    -1 ]
 *      [-1     0     0     0 ]
 *      [ 0    -1     0     0 ]
 *
 * g1 = [ 0     0     0    -i ]
 *      [ 0     0    -i     0 ]
 *      [ 0     i     0     0 ]
 *      [ i     0     0     0 ]
 *
 * g2 = [ 0     0     0    -1 ]
 *      [ 0     0     1     0 ]
 *      [ 0     1     0     0 ]
 *      [-1     0     0     0 ]
 * 
 * g3 = [ 0     0    -i     0 ]
 *      [ 0     0     0     i ]
 *      [ i     0     0     0 ]
 *      [ 0    -i     0     0 ]
 *
 *
 ******************************************************************************/

#define LINALG_SALG_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "global.h"

void mulg0_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   { 
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c3.c1.re; (*s).c3.c1.re =-tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im =-(*s).c3.c1.im; (*s).c3.c1.im =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c3.c2.re; (*s).c3.c2.re =-tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im =-(*s).c3.c2.im; (*s).c3.c2.im =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c3.c3.re; (*s).c3.c3.re =-tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im =-(*s).c3.c3.im; (*s).c3.c3.im =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re =-(*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im =-(*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re =-(*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im =-(*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re =-(*s).c4.c3.re; (*s).c4.c3.re =-tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im =-(*s).c4.c3.im; (*s).c4.c3.im =-tmp;
   }
}

void mulg1_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re = (*s).c4.c1.im; (*s).c4.c1.im = tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im =-(*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re = (*s).c4.c2.im; (*s).c4.c2.im = tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im =-(*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re = (*s).c4.c3.im; (*s).c4.c3.im = tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im =-(*s).c4.c3.re; (*s).c4.c3.re =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re = (*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im =-(*s).c3.c1.re; (*s).c3.c1.re =-tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re = (*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im =-(*s).c3.c2.re; (*s).c3.c2.re =-tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re = (*s).c3.c3.im; (*s).c3.c3.im = tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im =-(*s).c3.c3.re; (*s).c3.c3.re =-tmp;
   }
}

void mulg2_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im =-(*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im =-(*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c4.c3.re; (*s).c4.c3.re =-tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im =-(*s).c4.c3.im; (*s).c4.c3.im =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re = (*s).c3.c1.re; (*s).c3.c1.re = tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im = (*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re = (*s).c3.c2.re; (*s).c3.c2.re = tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im = (*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re = (*s).c3.c3.re; (*s).c3.c3.re = tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im = (*s).c3.c3.im; (*s).c3.c3.im = tmp;
   }
}

void mulg3_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re = (*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im =-(*s).c3.c1.re; (*s).c3.c1.re =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re = (*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im =-(*s).c3.c2.re; (*s).c3.c2.re =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re = (*s).c3.c3.im; (*s).c3.c3.im = tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im =-(*s).c3.c3.re; (*s).c3.c3.re =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re =-(*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im = (*s).c4.c1.re; (*s).c4.c1.re = tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re =-(*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im = (*s).c4.c2.re; (*s).c4.c2.re = tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re =-(*s).c4.c3.im; (*s).c4.c3.im =-tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im = (*s).c4.c3.re; (*s).c4.c3.re = tmp;
   }
}

void mulg0g1_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re = (*s).c2.c1.im; (*s).c2.c1.im =-tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im =-(*s).c2.c1.re; (*s).c2.c1.re = tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re = (*s).c2.c2.im; (*s).c2.c2.im =-tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im =-(*s).c2.c2.re; (*s).c2.c2.re = tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re = (*s).c2.c3.im; (*s).c2.c3.im =-tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im =-(*s).c2.c3.re; (*s).c2.c3.re = tmp;

      tmp = (*s).c3.c1.re; (*s).c3.c1.re =-(*s).c4.c1.im; (*s).c4.c1.im = tmp;
      tmp = (*s).c3.c1.im; (*s).c3.c1.im = (*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c3.c2.re; (*s).c3.c2.re =-(*s).c4.c2.im; (*s).c4.c2.im = tmp;
      tmp = (*s).c3.c2.im; (*s).c3.c2.im = (*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c3.c3.re; (*s).c3.c3.re =-(*s).c4.c3.im; (*s).c4.c3.im = tmp;
      tmp = (*s).c3.c3.im; (*s).c3.c3.im = (*s).c4.c3.re; (*s).c4.c3.re =-tmp;
   }
}

void mulg0g2_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c2.c1.re; (*s).c2.c1.re = tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im =-(*s).c2.c1.im; (*s).c2.c1.im = tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c2.c2.re; (*s).c2.c2.re = tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im =-(*s).c2.c2.im; (*s).c2.c2.im = tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c2.c3.re; (*s).c2.c3.re = tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im =-(*s).c2.c3.im; (*s).c2.c3.im = tmp;

      tmp = (*s).c3.c1.re; (*s).c3.c1.re = (*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c3.c1.im; (*s).c3.c1.im = (*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c3.c2.re; (*s).c3.c2.re = (*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c3.c2.im; (*s).c3.c2.im = (*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c3.c3.re; (*s).c3.c3.re = (*s).c4.c3.re; (*s).c4.c3.re =-tmp;
      tmp = (*s).c3.c3.im; (*s).c3.c3.im = (*s).c4.c3.im; (*s).c4.c3.im =-tmp;
   }
}

void mulg0g3_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re = (*s).c1.c1.im; (*s).c1.c1.im =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re = (*s).c1.c2.im; (*s).c1.c2.im =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re = (*s).c1.c3.im; (*s).c1.c3.im =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re =-(*s).c2.c1.im; (*s).c2.c1.im = tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re =-(*s).c2.c2.im; (*s).c2.c2.im = tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re =-(*s).c2.c3.im; (*s).c2.c3.im = tmp;
      
      tmp = (*s).c3.c1.re; (*s).c3.c1.re =-(*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c3.c2.re; (*s).c3.c2.re =-(*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c3.c3.re; (*s).c3.c3.re =-(*s).c3.c3.im; (*s).c3.c3.im = tmp;

      tmp = (*s).c4.c1.re; (*s).c4.c1.re = (*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c4.c2.re; (*s).c4.c2.re = (*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c4.c3.re; (*s).c4.c3.re = (*s).c4.c3.im; (*s).c4.c3.im =-tmp;
   }
}

void mulg0g5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re = (*s).c3.c1.re; (*s).c3.c1.re =-tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im = (*s).c3.c1.im; (*s).c3.c1.im =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re = (*s).c3.c2.re; (*s).c3.c2.re =-tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im = (*s).c3.c2.im; (*s).c3.c2.im =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re = (*s).c3.c3.re; (*s).c3.c3.re =-tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im = (*s).c3.c3.im; (*s).c3.c3.im =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re = (*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im = (*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re = (*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im = (*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re = (*s).c4.c3.re; (*s).c4.c3.re =-tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im = (*s).c4.c3.im; (*s).c4.c3.im =-tmp;
   }
}

void mulg1g2_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c1.c1.im; (*s).c1.c1.im = tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c1.c2.im; (*s).c1.c2.im = tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c1.c3.im; (*s).c1.c3.im = tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re = (*s).c2.c1.im; (*s).c2.c1.im =-tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re = (*s).c2.c2.im; (*s).c2.c2.im =-tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re = (*s).c2.c3.im; (*s).c2.c3.im =-tmp;

      tmp = (*s).c3.c1.re; (*s).c3.c1.re =-(*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c3.c2.re; (*s).c3.c2.re =-(*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c3.c3.re; (*s).c3.c3.re =-(*s).c3.c3.im; (*s).c3.c3.im = tmp;

      tmp = (*s).c4.c1.re; (*s).c4.c1.re = (*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c4.c2.re; (*s).c4.c2.re = (*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c4.c3.re; (*s).c4.c3.re = (*s).c4.c3.im; (*s).c4.c3.im =-tmp;
   }
}

void mulg1g3_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c2.c1.re; (*s).c2.c1.re = tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im =-(*s).c2.c1.im; (*s).c2.c1.im = tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c2.c2.re; (*s).c2.c2.re = tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im =-(*s).c2.c2.im; (*s).c2.c2.im = tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c2.c3.re; (*s).c2.c3.re = tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im =-(*s).c2.c3.im; (*s).c2.c3.im = tmp;

      tmp = (*s).c3.c1.re; (*s).c3.c1.re =-(*s).c4.c1.re; (*s).c4.c1.re = tmp;
      tmp = (*s).c3.c1.im; (*s).c3.c1.im =-(*s).c4.c1.im; (*s).c4.c1.im = tmp;
      tmp = (*s).c3.c2.re; (*s).c3.c2.re =-(*s).c4.c2.re; (*s).c4.c2.re = tmp;
      tmp = (*s).c3.c2.im; (*s).c3.c2.im =-(*s).c4.c2.im; (*s).c4.c2.im = tmp;
      tmp = (*s).c3.c3.re; (*s).c3.c3.re =-(*s).c4.c3.re; (*s).c4.c3.re = tmp;
      tmp = (*s).c3.c3.im; (*s).c3.c3.im =-(*s).c4.c3.im; (*s).c4.c3.im = tmp;
   }
}

void mulg1g5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c4.c1.im; (*s).c4.c1.im = tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im = (*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c4.c2.im; (*s).c4.c2.im = tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im = (*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c4.c3.im; (*s).c4.c3.im = tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im = (*s).c4.c3.re; (*s).c4.c3.re =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re =-(*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im = (*s).c3.c1.re; (*s).c3.c1.re =-tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re =-(*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im = (*s).c3.c2.re; (*s).c3.c2.re =-tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re =-(*s).c3.c3.im; (*s).c3.c3.im = tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im = (*s).c3.c3.re; (*s).c3.c3.re =-tmp;
   }
}

void mulg2g3_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c2.c1.im; (*s).c2.c1.im = tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im = (*s).c2.c1.re; (*s).c2.c1.re =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c2.c2.im; (*s).c2.c2.im = tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im = (*s).c2.c2.re; (*s).c2.c2.re =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c2.c3.im; (*s).c2.c3.im = tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im = (*s).c2.c3.re; (*s).c2.c3.re =-tmp;

      tmp = (*s).c3.c1.re; (*s).c3.c1.re =-(*s).c4.c1.im; (*s).c4.c1.im = tmp;
      tmp = (*s).c3.c1.im; (*s).c3.c1.im = (*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c3.c2.re; (*s).c3.c2.re =-(*s).c4.c2.im; (*s).c4.c2.im = tmp;
      tmp = (*s).c3.c2.im; (*s).c3.c2.im = (*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c3.c3.re; (*s).c3.c3.re =-(*s).c4.c3.im; (*s).c4.c3.im = tmp;
      tmp = (*s).c3.c3.im; (*s).c3.c3.im = (*s).c4.c3.re; (*s).c4.c3.re =-tmp;
   }
}

void mulg2g5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re = (*s).c4.c1.re; (*s).c4.c1.re =-tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im = (*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re = (*s).c4.c2.re; (*s).c4.c2.re =-tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im = (*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re = (*s).c4.c3.re; (*s).c4.c3.re =-tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im = (*s).c4.c3.im; (*s).c4.c3.im =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re =-(*s).c3.c1.re; (*s).c3.c1.re = tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im =-(*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re =-(*s).c3.c2.re; (*s).c3.c2.re = tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im =-(*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re =-(*s).c3.c3.re; (*s).c3.c3.re = tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im =-(*s).c3.c3.im; (*s).c3.c3.im = tmp;
   }
}

void mulg3g5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;
   double tmp;

   sm=s+vol;

   for (;s<sm;s++)
   {
      tmp = (*s).c1.c1.re; (*s).c1.c1.re =-(*s).c3.c1.im; (*s).c3.c1.im = tmp;
      tmp = (*s).c1.c1.im; (*s).c1.c1.im = (*s).c3.c1.re; (*s).c3.c1.re =-tmp;
      tmp = (*s).c1.c2.re; (*s).c1.c2.re =-(*s).c3.c2.im; (*s).c3.c2.im = tmp;
      tmp = (*s).c1.c2.im; (*s).c1.c2.im = (*s).c3.c2.re; (*s).c3.c2.re =-tmp;
      tmp = (*s).c1.c3.re; (*s).c1.c3.re =-(*s).c3.c3.im; (*s).c3.c3.im = tmp;
      tmp = (*s).c1.c3.im; (*s).c1.c3.im = (*s).c3.c3.re; (*s).c3.c3.re =-tmp;

      tmp = (*s).c2.c1.re; (*s).c2.c1.re = (*s).c4.c1.im; (*s).c4.c1.im =-tmp;
      tmp = (*s).c2.c1.im; (*s).c2.c1.im =-(*s).c4.c1.re; (*s).c4.c1.re = tmp;
      tmp = (*s).c2.c2.re; (*s).c2.c2.re = (*s).c4.c2.im; (*s).c4.c2.im =-tmp;
      tmp = (*s).c2.c2.im; (*s).c2.c2.im =-(*s).c4.c2.re; (*s).c4.c2.re = tmp;
      tmp = (*s).c2.c3.re; (*s).c2.c3.re = (*s).c4.c3.im; (*s).c4.c3.im =-tmp;
      tmp = (*s).c2.c3.im; (*s).c2.c3.im =-(*s).c4.c3.re; (*s).c4.c3.re = tmp;
   }
}
