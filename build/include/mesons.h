/*******************************************************************************
*
* File mesons.h
*
* Copyright (C) 2013 Tomasz Korzec
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef MESONS_H
#define MESONS_H

#define Z2_NOISE 0
#define GAUSS_NOISE 1
#define U1_NOISE 2

#define GAMMA0_TYPE 0
#define GAMMA1_TYPE 1
#define GAMMA2_TYPE 2
#define GAMMA3_TYPE 3
#define GAMMA5_TYPE 5
#define ONE_TYPE 6
#define GAMMA0GAMMA1_TYPE 7
#define GAMMA0GAMMA2_TYPE 8
#define GAMMA0GAMMA3_TYPE 9
#define GAMMA0GAMMA5_TYPE 10
#define GAMMA1GAMMA2_TYPE 11
#define GAMMA1GAMMA3_TYPE 12
#define GAMMA1GAMMA5_TYPE 13
#define GAMMA2GAMMA3_TYPE 14
#define GAMMA2GAMMA5_TYPE 15
#define GAMMA3GAMMA5_TYPE 16
#define MAX_TYPE 17

#define mesons_RELEASE "mesons v1.2rc"


/* LINALG_SALG_DBLE_C */
extern void mulg0_dble(int vol,spinor_dble *s);
extern void mulg1_dble(int vol,spinor_dble *s);
extern void mulg2_dble(int vol,spinor_dble *s);
extern void mulg3_dble(int vol,spinor_dble *s);
extern void mulg0g1_dble(int vol,spinor_dble *s);
extern void mulg0g2_dble(int vol,spinor_dble *s);
extern void mulg0g3_dble(int vol,spinor_dble *s);
extern void mulg0g5_dble(int vol,spinor_dble *s);
extern void mulg1g2_dble(int vol,spinor_dble *s);
extern void mulg1g3_dble(int vol,spinor_dble *s);
extern void mulg1g5_dble(int vol,spinor_dble *s);
extern void mulg2g3_dble(int vol,spinor_dble *s);
extern void mulg2g5_dble(int vol,spinor_dble *s);
extern void mulg3g5_dble(int vol,spinor_dble *s);

/* SFLDS_SFLDS_C */
extern void random_Z2_sd(int vol,spinor_dble *sd);
extern void random_U1_sd(int vol,spinor_dble *sd);
extern void assign_msd2sd(int vol,spinor_dble *sd,spinor_dble *rd);

/* UTILS_MUTILS_C */
extern long read_line_opt(char *tag,char *defaultline,char *format, ...);

#endif
