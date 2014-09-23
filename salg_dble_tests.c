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
 *
 * tests for gamma*spinor
 *
 * (needs C99 compound literals)
 * 
 *******************************************************************************/

#include "linalg.h"
#include "su3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mesons.h>

#define VOL 5
#define tol 1e-15

#define MAX(a,b) ((a)=(((a)>(b))?(a):(b)))

#define CMUL(a,b) ((complex_dble) {(a).re*(b).re-(a).im*(b).im, (a).re*(b).im+(a).im*(b).re})

#define CADD(a,b) ((complex_dble) {(a).re+(b).re, (a).im+(b).im})

double maxreldiff(spinor_dble *s1, spinor_dble *s2)
{
   double rd, mrd = 0;
   int v;
   for (v=0; v<VOL; v++)
   {
      rd=fabs((s1[v].c1.c1.re-s2[v].c1.c1.re)/s1[v].c1.c1.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c1.c1.im-s2[v].c1.c1.im)/s1[v].c1.c1.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c1.c2.re-s2[v].c1.c2.re)/s1[v].c1.c2.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c1.c2.im-s2[v].c1.c2.im)/s1[v].c1.c2.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c1.c3.re-s2[v].c1.c3.re)/s1[v].c1.c3.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c1.c3.im-s2[v].c1.c3.im)/s1[v].c1.c3.im);
      MAX(mrd,rd);

      rd=fabs((s1[v].c2.c1.re-s2[v].c2.c1.re)/s1[v].c2.c1.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c2.c1.im-s2[v].c2.c1.im)/s1[v].c2.c1.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c2.c2.re-s2[v].c2.c2.re)/s1[v].c2.c2.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c2.c2.im-s2[v].c2.c2.im)/s1[v].c2.c2.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c2.c3.re-s2[v].c2.c3.re)/s1[v].c2.c3.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c2.c3.im-s2[v].c2.c3.im)/s1[v].c2.c3.im);
      MAX(mrd,rd);

      rd=fabs((s1[v].c3.c1.re-s2[v].c3.c1.re)/s1[v].c3.c1.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c3.c1.im-s2[v].c3.c1.im)/s1[v].c3.c1.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c3.c2.re-s2[v].c3.c2.re)/s1[v].c3.c2.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c3.c2.im-s2[v].c3.c2.im)/s1[v].c3.c2.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c3.c3.re-s2[v].c3.c3.re)/s1[v].c3.c3.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c3.c3.im-s2[v].c3.c3.im)/s1[v].c3.c3.im);
      MAX(mrd,rd);

      rd=fabs((s1[v].c4.c1.re-s2[v].c4.c1.re)/s1[v].c4.c1.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c4.c1.im-s2[v].c4.c1.im)/s1[v].c4.c1.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c4.c2.re-s2[v].c4.c2.re)/s1[v].c4.c2.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c4.c2.im-s2[v].c4.c2.im)/s1[v].c4.c2.im);
      MAX(mrd,rd);
      rd=fabs((s1[v].c4.c3.re-s2[v].c4.c3.re)/s1[v].c4.c3.re);
      MAX(mrd,rd);
      rd=fabs((s1[v].c4.c3.im-s2[v].c4.c3.im)/s1[v].c4.c3.im);
      MAX(mrd,rd);
   }
   return(mrd);
}

void fullmul(complex_dble g[4][4], spinor_dble *sp1, spinor_dble *sp2)
{
   int v;
   for (v=0; v<VOL; v++)
   {
      sp2[v].c1.c1 = CADD(CADD(CADD(CMUL(g[0][0],sp1[v].c1.c1),
                                    CMUL(g[0][1],sp1[v].c2.c1)),
                               CMUL(g[0][2],sp1[v].c3.c1)),
                          CMUL(g[0][3],sp1[v].c4.c1));
      sp2[v].c1.c2 = CADD(CADD(CADD(CMUL(g[0][0],sp1[v].c1.c2),
                                    CMUL(g[0][1],sp1[v].c2.c2)),
                               CMUL(g[0][2],sp1[v].c3.c2)),
                          CMUL(g[0][3],sp1[v].c4.c2));
      sp2[v].c1.c3 = CADD(CADD(CADD(CMUL(g[0][0],sp1[v].c1.c3),
                                    CMUL(g[0][1],sp1[v].c2.c3)),
                               CMUL(g[0][2],sp1[v].c3.c3)),
                          CMUL(g[0][3],sp1[v].c4.c3));

      sp2[v].c2.c1 = CADD(CADD(CADD(CMUL(g[1][0],sp1[v].c1.c1),
                                    CMUL(g[1][1],sp1[v].c2.c1)),
                               CMUL(g[1][2],sp1[v].c3.c1)),
                          CMUL(g[1][3],sp1[v].c4.c1));
      sp2[v].c2.c2 = CADD(CADD(CADD(CMUL(g[1][0],sp1[v].c1.c2),
                                    CMUL(g[1][1],sp1[v].c2.c2)),
                               CMUL(g[1][2],sp1[v].c3.c2)),
                          CMUL(g[1][3],sp1[v].c4.c2));
      sp2[v].c2.c3 = CADD(CADD(CADD(CMUL(g[1][0],sp1[v].c1.c3),
                                    CMUL(g[1][1],sp1[v].c2.c3)),
                               CMUL(g[1][2],sp1[v].c3.c3)),
                          CMUL(g[1][3],sp1[v].c4.c3));

      sp2[v].c3.c1 = CADD(CADD(CADD(CMUL(g[2][0],sp1[v].c1.c1),
                                    CMUL(g[2][1],sp1[v].c2.c1)),
                               CMUL(g[2][2],sp1[v].c3.c1)),
                          CMUL(g[2][3],sp1[v].c4.c1));
      sp2[v].c3.c2 = CADD(CADD(CADD(CMUL(g[2][0],sp1[v].c1.c2),
                                    CMUL(g[2][1],sp1[v].c2.c2)),
                               CMUL(g[2][2],sp1[v].c3.c2)),
                          CMUL(g[2][3],sp1[v].c4.c2));
      sp2[v].c3.c3 = CADD(CADD(CADD(CMUL(g[2][0],sp1[v].c1.c3),
                                    CMUL(g[2][1],sp1[v].c2.c3)),
                               CMUL(g[2][2],sp1[v].c3.c3)),
                          CMUL(g[2][3],sp1[v].c4.c3));

      sp2[v].c4.c1 = CADD(CADD(CADD(CMUL(g[3][0],sp1[v].c1.c1),
                                    CMUL(g[3][1],sp1[v].c2.c1)),
                               CMUL(g[3][2],sp1[v].c3.c1)),
                          CMUL(g[3][3],sp1[v].c4.c1));
      sp2[v].c4.c2 = CADD(CADD(CADD(CMUL(g[3][0],sp1[v].c1.c2),
                                    CMUL(g[3][1],sp1[v].c2.c2)),
                               CMUL(g[3][2],sp1[v].c3.c2)),
                          CMUL(g[3][3],sp1[v].c4.c2));
      sp2[v].c4.c3 = CADD(CADD(CADD(CMUL(g[3][0],sp1[v].c1.c3),
                                    CMUL(g[3][1],sp1[v].c2.c3)),
                               CMUL(g[3][2],sp1[v].c3.c3)),
                          CMUL(g[3][3],sp1[v].c4.c3));
   }
}

int main(int argc,char *argv[])
{
   complex_dble g0[4][4] = {{{ 0,0},{ 0,0},{-1,0},{ 0,0}},
                            {{ 0,0},{ 0,0},{ 0,0},{-1,0}},
                            {{-1,0},{ 0,0},{ 0,0},{ 0,0}},
                            {{ 0,0},{-1,0},{ 0,0},{ 0,0}}};

   complex_dble g1[4][4] = {{{ 0,0},{ 0,0},{ 0, 0},{ 0,-1}},
                            {{ 0,0},{ 0,0},{ 0,-1},{ 0,0}},
                            {{ 0,0},{ 0,1},{ 0, 0},{ 0,0}},
                            {{ 0,1},{ 0,0},{ 0, 0},{ 0,0}}};

   complex_dble g2[4][4] = {{{ 0,0},{ 0,0},{ 0,0},{-1,0}},
                            {{ 0,0},{ 0,0},{ 1,0},{ 0,0}},
                            {{ 0,0},{ 1,0},{ 0,0},{ 0,0}},
                            {{-1,0},{ 0,0},{ 0,0},{ 0,0}}};

   complex_dble g3[4][4] = {{{ 0,0},{ 0, 0},{ 0,-1},{ 0,0}},
                            {{ 0,0},{ 0, 0},{ 0, 0},{ 0,1}},
                            {{ 0,1},{ 0, 0},{ 0, 0},{ 0,0}},
                            {{ 0,0},{ 0,-1},{ 0, 0},{ 0,0}}};

   spinor_dble sp1[VOL];
   spinor_dble sp2[VOL];
   spinor_dble sp3[VOL];

   int i;

   for (i=0;i<VOL;i++)
   {
      sp1[i].c1.c1.re=drand48(); sp1[i].c1.c2.re=drand48(); sp1[i].c1.c3.re=drand48();
      sp1[i].c2.c1.re=drand48(); sp1[i].c2.c2.re=drand48(); sp1[i].c2.c3.re=drand48();
      sp1[i].c3.c1.re=drand48(); sp1[i].c3.c2.re=drand48(); sp1[i].c3.c3.re=drand48();
      sp1[i].c4.c1.re=drand48(); sp1[i].c4.c2.re=drand48(); sp1[i].c4.c3.re=drand48();

      sp1[i].c1.c1.im=drand48(); sp1[i].c1.c2.im=drand48(); sp1[i].c1.c3.im=drand48();
      sp1[i].c2.c1.im=drand48(); sp1[i].c2.c2.im=drand48(); sp1[i].c2.c3.im=drand48();
      sp1[i].c3.c1.im=drand48(); sp1[i].c3.c2.im=drand48(); sp1[i].c3.c3.im=drand48();
      sp1[i].c4.c1.im=drand48(); sp1[i].c4.c2.im=drand48(); sp1[i].c4.c3.im=drand48();
   }

   printf("Comparison with full matrix X vector\n")
   fullmul(g0,sp1,sp2);
   mulg0_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g0 passed\n");
   else
      printf("<--------------g0 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g1,sp1,sp2);
   mulg1_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g1 passed\n");
   else
      printf("<--------------g1 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g2,sp1,sp2);
   mulg2_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g2 passed\n");
   else
      printf("<--------------g2 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g3,sp1,sp2);
   mulg3_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g3 passed\n");
   else
      printf("<--------------g3 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   /*
   fullmul(g3,sp1,sp3);
   fullmul(g2,sp3,sp2);
   fullmul(g1,sp2,sp3);
   fullmul(g0,sp3,sp2);
   mulg5_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g5 passed\n");
   else
      printf("<--------------g5 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));
   */
   
   fullmul(g1,sp1,sp3);
   fullmul(g0,sp3,sp2);
   mulg0g1_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g0g1 passed\n");
   else
      printf("<------------g0g3 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g2,sp1,sp3);
   fullmul(g0,sp3,sp2);
   mulg0g2_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g0g2 passed\n");
   else
      printf("<------------g0g2 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g3,sp1,sp3);
   fullmul(g0,sp3,sp2);
   mulg0g3_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g0g3 passed\n");
   else
      printf("<------------g0g3 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g3,sp1,sp2);
   fullmul(g2,sp2,sp3);
   fullmul(g1,sp3,sp2);
   mulg0g5_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g0g5 passed\n");
   else
      printf("<------------g0g5 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g2,sp1,sp3);
   fullmul(g1,sp3,sp2);
   mulg1g2_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g1g2 passed\n");
   else
      printf("<------------g1g2 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g3,sp1,sp3);
   fullmul(g1,sp3,sp2);
   mulg1g3_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g1g3 passed\n");
   else
      printf("<------------g1g3 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g3,sp1,sp2);
   fullmul(g0,sp2,sp3);
   fullmul(g2,sp3,sp2);
   mulg1g5_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g1g5 passed\n");
   else
      printf("<------------g1g5 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));
   
   fullmul(g3,sp1,sp3);
   fullmul(g2,sp3,sp2);
   mulg2g3_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g2g3 passed\n");
   else
      printf("<------------g2g3 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g3,sp1,sp2);
   fullmul(g1,sp2,sp3);
   fullmul(g0,sp3,sp2);
   mulg2g5_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g2g5 passed\n");
   else
      printf("<------------g2g5 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   fullmul(g2,sp1,sp2);
   fullmul(g0,sp2,sp3);
   fullmul(g1,sp3,sp2);
   mulg3g5_dble(VOL,sp1);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g3g5 passed\n");
   else
      printf("<------------g3g5 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   /* test all commutation relations */
   printf("\ncommutation relations:\n")
   /* first "copy" sp1 to sp2 */
   fullmul(g0,sp1,sp3);
   fullmul(g0,sp3,sp2);
   if(maxreldiff(sp1,sp2)>=tol)
      printf("<------------g0g0 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   mulg0_dble(VOL,sp2);
   mulg0_dble(VOL,sp2);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g0g0 passed\n");
   else
      printf("<------------g0g0 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   mulg1_dble(VOL,sp2);
   mulg1_dble(VOL,sp2);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g1g1 passed\n");
   else
      printf("<------------g1g1 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   mulg2_dble(VOL,sp2);
   mulg2_dble(VOL,sp2);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g2g2 passed\n");
   else
      printf("<------------g2g2 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   mulg3_dble(VOL,sp2);
   mulg3_dble(VOL,sp2);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g3g3 passed\n");
   else
      printf("<------------g3g3 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));

   /*
   mulg5_dble(VOL,sp2);
   mulg5_dble(VOL,sp2);
   if(maxreldiff(sp1,sp2)<tol)
      printf("g5g5 passed\n");
   else
      printf("<------------g5g5 failed! max rel diff=%f\n",maxreldiff(sp1,sp2));
   */
   
   return 0;
}
