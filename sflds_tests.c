/*******************************************************************************
 *
 * Copyright (C) 2013, 2014 Tomasz Korzec
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
 * tests for random spinor fields
 *
 * 
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "archive.h"
#include "sflds.h"
#include "linalg.h"
#include "dirac.h"
#include "sap.h"
#include "dfl.h"
#include "forces.h"
#include "version.h"
#include "global.h"
#include "mesons.h"


#define VOL 1e5

int main(int argc,char *argv[])
{
   spinor_dble *s;
   complex_dble m1,m2,m3;
   int i, my_rank;
   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   start_ranlux(0,137);

   s = malloc(VOL*sizeof(spinor_dble));
   if (s==NULL)
   {
      fprintf(stderr,"Out of memory.\n");
      exit(EXIT_FAILURE);
   }

   random_Z2_sd(VOL, s);
   /* estimate moments of Z2 distribution and compare with exact results */

   m1.re=0; m1.im=0;
   m2.re=0; m2.im=0;
   m3.re=0; m3.im=0;
   
   for(i=0; i<VOL; i++)
   {
      m1.re += s[i].c1.c1.re + s[i].c1.c2.re + s[i].c1.c2.re
              +s[i].c2.c1.re + s[i].c2.c2.re + s[i].c2.c2.re
              +s[i].c3.c1.re + s[i].c3.c2.re + s[i].c3.c2.re
              +s[i].c4.c1.re + s[i].c4.c2.re + s[i].c4.c2.re;
      m1.im += s[i].c1.c1.im + s[i].c1.c2.im + s[i].c1.c2.im
              +s[i].c2.c1.im + s[i].c2.c2.im + s[i].c2.c2.im
              +s[i].c3.c1.im + s[i].c3.c2.im + s[i].c3.c2.im
              +s[i].c4.c1.im + s[i].c4.c2.im + s[i].c4.c2.im;

      m2.re += s[i].c1.c1.re*s[i].c1.c1.re+s[i].c1.c1.im*s[i].c1.c1.im
              +s[i].c1.c2.re*s[i].c1.c2.re+s[i].c1.c2.im*s[i].c1.c2.im
              +s[i].c1.c3.re*s[i].c1.c3.re+s[i].c1.c3.im*s[i].c1.c3.im
              +s[i].c2.c1.re*s[i].c2.c1.re+s[i].c2.c1.im*s[i].c2.c1.im
              +s[i].c2.c2.re*s[i].c2.c2.re+s[i].c2.c2.im*s[i].c2.c2.im
              +s[i].c2.c3.re*s[i].c2.c3.re+s[i].c2.c3.im*s[i].c2.c3.im
              +s[i].c3.c1.re*s[i].c3.c1.re+s[i].c3.c1.im*s[i].c3.c1.im
              +s[i].c3.c2.re*s[i].c3.c2.re+s[i].c3.c2.im*s[i].c3.c2.im
              +s[i].c3.c3.re*s[i].c3.c3.re+s[i].c3.c3.im*s[i].c3.c3.im
              +s[i].c4.c1.re*s[i].c4.c1.re+s[i].c4.c1.im*s[i].c4.c1.im
              +s[i].c4.c2.re*s[i].c4.c2.re+s[i].c4.c2.im*s[i].c4.c2.im
              +s[i].c4.c3.re*s[i].c4.c3.re+s[i].c4.c3.im*s[i].c4.c3.im;
   }
   for(i=0; i<VOL-1; i++)
   {
      m3.re += s[i].c1.c1.re*s[i+1].c1.c1.re+s[i].c1.c1.im*s[i+1].c1.c1.im
              +s[i].c1.c2.re*s[i+1].c1.c2.re+s[i].c1.c2.im*s[i+1].c1.c2.im
              +s[i].c1.c3.re*s[i+1].c1.c3.re+s[i].c1.c3.im*s[i+1].c1.c3.im
              +s[i].c2.c1.re*s[i+1].c2.c1.re+s[i].c2.c1.im*s[i+1].c2.c1.im
              +s[i].c2.c2.re*s[i+1].c2.c2.re+s[i].c2.c2.im*s[i+1].c2.c2.im
              +s[i].c2.c3.re*s[i+1].c2.c3.re+s[i].c2.c3.im*s[i+1].c2.c3.im
              +s[i].c3.c1.re*s[i+1].c3.c1.re+s[i].c3.c1.im*s[i+1].c3.c1.im
              +s[i].c3.c2.re*s[i+1].c3.c2.re+s[i].c3.c2.im*s[i+1].c3.c2.im
              +s[i].c3.c3.re*s[i+1].c3.c3.re+s[i].c3.c3.im*s[i+1].c3.c3.im
              +s[i].c4.c1.re*s[i+1].c4.c1.re+s[i].c4.c1.im*s[i+1].c4.c1.im
              +s[i].c4.c2.re*s[i+1].c4.c2.re+s[i].c4.c2.im*s[i+1].c4.c2.im
              +s[i].c4.c3.re*s[i+1].c4.c3.re+s[i].c4.c3.im*s[i+1].c4.c3.im;
   }
   
   m1.re /= 12*VOL;
   m1.im /= 12*VOL;
   m2.re /= 12*VOL;
   m2.im /= 12*VOL;
   m3.re /= 12*VOL;
   m3.im /= 12*VOL;

   printf("%i: Z2-noise:\n",my_rank);
   printf("%i:  <s>          = %+f %+fi   +- %f\n",my_rank,m1.re,m1.im,1/sqrt(12*VOL));
   printf("%i:  <conj(s)*s>  = %+f %+fi   +- 0\n",my_rank,m2.re,m2.im);
   printf("%i:  <conj(s1)*s2>= %+f %+fi   +- %f\n",my_rank,m3.re,m3.im,1/sqrt(12*(VOL-1)));

   if(my_rank==0) printf("\n\n");

   random_U1_sd(VOL, s);
   /* estimate moments of U1 distribution and compare with exact results */
   m1.re=0; m1.im=0;
   m2.re=0; m2.im=0;
   m3.re=0; m3.im=0;
   for(i=0; i<VOL; i++)
   {
      m1.re += s[i].c1.c1.re + s[i].c1.c2.re + s[i].c1.c2.re
              +s[i].c2.c1.re + s[i].c2.c2.re + s[i].c2.c2.re
              +s[i].c3.c1.re + s[i].c3.c2.re + s[i].c3.c2.re
              +s[i].c4.c1.re + s[i].c4.c2.re + s[i].c4.c2.re;
      m1.im += s[i].c1.c1.im + s[i].c1.c2.im + s[i].c1.c2.im
              +s[i].c2.c1.im + s[i].c2.c2.im + s[i].c2.c2.im
              +s[i].c3.c1.im + s[i].c3.c2.im + s[i].c3.c2.im
              +s[i].c4.c1.im + s[i].c4.c2.im + s[i].c4.c2.im;

      m2.re += s[i].c1.c1.re*s[i].c1.c1.re+s[i].c1.c1.im*s[i].c1.c1.im
              +s[i].c1.c2.re*s[i].c1.c2.re+s[i].c1.c2.im*s[i].c1.c2.im
              +s[i].c1.c3.re*s[i].c1.c3.re+s[i].c1.c3.im*s[i].c1.c3.im
              +s[i].c2.c1.re*s[i].c2.c1.re+s[i].c2.c1.im*s[i].c2.c1.im
              +s[i].c2.c2.re*s[i].c2.c2.re+s[i].c2.c2.im*s[i].c2.c2.im
              +s[i].c2.c3.re*s[i].c2.c3.re+s[i].c2.c3.im*s[i].c2.c3.im
              +s[i].c3.c1.re*s[i].c3.c1.re+s[i].c3.c1.im*s[i].c3.c1.im
              +s[i].c3.c2.re*s[i].c3.c2.re+s[i].c3.c2.im*s[i].c3.c2.im
              +s[i].c3.c3.re*s[i].c3.c3.re+s[i].c3.c3.im*s[i].c3.c3.im
              +s[i].c4.c1.re*s[i].c4.c1.re+s[i].c4.c1.im*s[i].c4.c1.im
              +s[i].c4.c2.re*s[i].c4.c2.re+s[i].c4.c2.im*s[i].c4.c2.im
              +s[i].c4.c3.re*s[i].c4.c3.re+s[i].c4.c3.im*s[i].c4.c3.im;
   }
   for(i=0; i<VOL-1; i++)
   {
      m3.re += s[i].c1.c1.re*s[i+1].c1.c1.re+s[i].c1.c1.im*s[i+1].c1.c1.im
              +s[i].c1.c2.re*s[i+1].c1.c2.re+s[i].c1.c2.im*s[i+1].c1.c2.im
              +s[i].c1.c3.re*s[i+1].c1.c3.re+s[i].c1.c3.im*s[i+1].c1.c3.im
              +s[i].c2.c1.re*s[i+1].c2.c1.re+s[i].c2.c1.im*s[i+1].c2.c1.im
              +s[i].c2.c2.re*s[i+1].c2.c2.re+s[i].c2.c2.im*s[i+1].c2.c2.im
              +s[i].c2.c3.re*s[i+1].c2.c3.re+s[i].c2.c3.im*s[i+1].c2.c3.im
              +s[i].c3.c1.re*s[i+1].c3.c1.re+s[i].c3.c1.im*s[i+1].c3.c1.im
              +s[i].c3.c2.re*s[i+1].c3.c2.re+s[i].c3.c2.im*s[i+1].c3.c2.im
              +s[i].c3.c3.re*s[i+1].c3.c3.re+s[i].c3.c3.im*s[i+1].c3.c3.im
              +s[i].c4.c1.re*s[i+1].c4.c1.re+s[i].c4.c1.im*s[i+1].c4.c1.im
              +s[i].c4.c2.re*s[i+1].c4.c2.re+s[i].c4.c2.im*s[i+1].c4.c2.im
              +s[i].c4.c3.re*s[i+1].c4.c3.re+s[i].c4.c3.im*s[i+1].c4.c3.im;

      m3.im += s[i].c1.c1.re*s[i+1].c1.c1.im-s[i].c1.c1.im*s[i+1].c1.c1.re
              +s[i].c1.c2.re*s[i+1].c1.c2.im-s[i].c1.c2.im*s[i+1].c1.c2.re
              +s[i].c1.c3.re*s[i+1].c1.c3.im-s[i].c1.c3.im*s[i+1].c1.c3.re
              +s[i].c2.c1.re*s[i+1].c2.c1.im-s[i].c2.c1.im*s[i+1].c2.c1.re
              +s[i].c2.c2.re*s[i+1].c2.c2.im-s[i].c2.c2.im*s[i+1].c2.c2.re
              +s[i].c2.c3.re*s[i+1].c2.c3.im-s[i].c2.c3.im*s[i+1].c2.c3.re
              +s[i].c3.c1.re*s[i+1].c3.c1.im-s[i].c3.c1.im*s[i+1].c3.c1.re
              +s[i].c3.c2.re*s[i+1].c3.c2.im-s[i].c3.c2.im*s[i+1].c3.c2.re
              +s[i].c3.c3.re*s[i+1].c3.c3.im-s[i].c3.c3.im*s[i+1].c3.c3.re
              +s[i].c4.c1.re*s[i+1].c4.c1.im-s[i].c4.c1.im*s[i+1].c4.c1.re
              +s[i].c4.c2.re*s[i+1].c4.c2.im-s[i].c4.c2.im*s[i+1].c4.c2.re
              +s[i].c4.c3.re*s[i+1].c4.c3.im-s[i].c4.c3.im*s[i+1].c4.c3.re;
   }

   
   m1.re /= 12*VOL;
   m1.im /= 12*VOL;
   m2.re /= 12*VOL;
   m2.im /= 12*VOL;
   m3.re /= 12*VOL;
   m3.im /= 12*VOL;

   printf("%i: U1-noise:\n",my_rank);
   printf("%i:  <s>          = %+f %+fi   +- %f%+fi\n",my_rank,m1.re,m1.im,1/sqrt(12*VOL),1/sqrt(12*VOL));
   printf("%i:  <conj(s)*s>  = %+f %+fi   +- 0\n",my_rank,m2.re,m2.im);
   printf("%i:  <conj(s1)*s2>= %+f %+fi   +- %f%+fi\n",my_rank,m3.re,m3.im,1/sqrt(12*(VOL-1)),1/sqrt(12*(VOL-1)));
   
   MPI_Finalize();
   exit(0);
}
