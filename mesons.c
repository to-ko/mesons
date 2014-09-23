/*******************************************************************************
*
* File mesons.c
*
* Copyright (C) 2013, 2014 Tomasz Korzec
*
* Based on openQCD, ms1 and ms4
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
* Computation of quark propagators
*
* Syntax: mesons -i <input file> [-noexp] [-a]
*
* For usage instructions see the file README.mesons
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static char line[NAME_SIZE+1];

#define MAX(n,m) \
   if ((n)<(m)) \
      (n)=(m)

static struct
{
   int ncorr;
   int nnoise;
   int tvals;
   int noisetype;
   double *kappa1;
   double *kappa2;
   int *type1;
   int *type2;
   int *x0;
   int *isreal;
} file_head;

static struct
{
   complex_dble *corr;
   complex_dble *corr_tmp;
   int nc;
} data;

static struct
{
   int nux0;     /* number of unique x0 values */
   int *ux0;     /* unique x0 values */
   int *nprop;   /* number of propagators at each x0 */
   int nmax;     /* max(nprop) */
   int **prop;   /* propagator index of each x0 and propagator */
   int **type;   /* type index of each x0 and propagator  */
} proplist;

static int my_rank,noexp,append,norng,endian;
static int first,last,step;
static int level,seed,nprop,ncorr,nnoise,noisetype,tvals;
static int *isps,*props1,*props2,*type1,*type2,*x0s;
static int ipgrd[2],*rlxs_state=NULL,*rlxd_state=NULL;
static double *kappas,*mus;


static char log_dir[NAME_SIZE],loc_dir[NAME_SIZE];
static char cnfg_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE],end_file[NAME_SIZE];
static char dat_file[NAME_SIZE],dat_save[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char rng_file[NAME_SIZE],rng_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],nbase[NAME_SIZE],outbase[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL,*fend=NULL,*fdat=NULL;


static void alloc_data(void)
{
   data.corr=malloc(file_head.ncorr*file_head.nnoise*file_head.tvals*
                                                          sizeof(complex_dble));
   data.corr_tmp=malloc(file_head.ncorr*file_head.nnoise*file_head.tvals*
                                                      sizeof(complex_dble));
   error((data.corr==NULL)||(data.corr_tmp==NULL),1,"alloc_data [mesons.c]",
         "Unable to allocate data arrays");
}


static void write_file_head(void)
{
   stdint_t istd[1];
   int iw=0;
   int i;
   double dbl[1];

   istd[0]=(stdint_t)(file_head.ncorr);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   iw=fwrite(istd,sizeof(stdint_t),1,fdat);

   istd[0]=(stdint_t)(file_head.nnoise);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   iw+=fwrite(istd,sizeof(stdint_t),1,fdat);

   istd[0]=(stdint_t)(file_head.tvals);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   iw+=fwrite(istd,sizeof(stdint_t),1,fdat);

   istd[0]=(stdint_t)(file_head.noisetype);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   iw+=fwrite(istd,sizeof(stdint_t),1,fdat);

   error_root(iw!=4,1,"write_file_head [mesons.c]",
              "Incorrect write count");
   for (i=0;i<file_head.ncorr;i++)
   {
      dbl[0] = file_head.kappa1[i];
      if (endian==BIG_ENDIAN)
      bswap_double(1,dbl);
      iw=fwrite(dbl,sizeof(double),1,fdat);

      dbl[0] = file_head.kappa2[i];
      if (endian==BIG_ENDIAN)
      bswap_double(1,dbl);
      iw+=fwrite(dbl,sizeof(double),1,fdat);

      istd[0]=(stdint_t)(file_head.type1[i]);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      iw+=fwrite(istd,sizeof(stdint_t),1,fdat);

      istd[0]=(stdint_t)(file_head.type2[i]);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      iw+=fwrite(istd,sizeof(stdint_t),1,fdat);

      istd[0]=(stdint_t)(file_head.x0[i]);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      iw+=fwrite(istd,sizeof(stdint_t),1,fdat);

      istd[0]=(stdint_t)(file_head.isreal[i]);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      iw+=fwrite(istd,sizeof(stdint_t),1,fdat);

      error_root(iw!=6,1,"write_file_head [mesons.c]",
              "Incorrect write count");
   }
}

static void check_file_head(void)
{
   int i,ir,ie;
   stdint_t istd[1];
   double dbl[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   ie=(istd[0]!=(stdint_t)(file_head.ncorr));

   ir+=fread(istd,sizeof(stdint_t),1,fdat);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   ie+=(istd[0]!=(stdint_t)(file_head.nnoise));

   ir+=fread(istd,sizeof(stdint_t),1,fdat);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   ie+=(istd[0]!=(stdint_t)(file_head.tvals));

   ir+=fread(istd,sizeof(stdint_t),1,fdat);
   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);
   ie+=(istd[0]!=(stdint_t)(file_head.noisetype));

   error_root(ir!=4,1,"check_file_head [mesons.c]",
              "Incorrect read count");
   error_root(ie!=0,1,"check_file_head [mesons.c]",
              "Unexpected value of ncorr, nnoise, tvals or noisetype");
   for (i=0;i<file_head.ncorr;i++)
   {
      ir=fread(dbl,sizeof(double),1,fdat);
      if (endian==BIG_ENDIAN)
      bswap_double(1,dbl);
      ie=(dbl[0]!=file_head.kappa1[i]);

      ir+=fread(dbl,sizeof(double),1,fdat);
      if (endian==BIG_ENDIAN)
      bswap_double(1,dbl);
      ie+=(dbl[0]!=file_head.kappa2[i]);

      ir+=fread(istd,sizeof(stdint_t),1,fdat);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      ie+=(istd[0]!=(stdint_t)(file_head.type1[i]));

      ir+=fread(istd,sizeof(stdint_t),1,fdat);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      ie+=(istd[0]!=(stdint_t)(file_head.type2[i]));

      ir+=fread(istd,sizeof(stdint_t),1,fdat);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      ie+=(istd[0]!=(stdint_t)(file_head.x0[i]));

      ir+=fread(istd,sizeof(stdint_t),1,fdat);
      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      ie+=(istd[0]!=(stdint_t)(file_head.isreal[i]));

      error_root(ir!=6,1,"check_file_head [mesons.c]",
              "Incorrect read count");
      error_root(ie!=0,1,"check_file_head [mesons.c]",
              "Unexpected value of kappa, type, x0 or isreal");
   }
}

static void write_data(void)
{
   int iw;
   int nw;
   int chunk;
   int icorr,i;

   if (my_rank==0)
   {
      fdat=fopen(dat_file,"ab");
      error_root(fdat==NULL,1,"write_data [mesons.c]",
                 "Unable to open dat file");

      nw = 1;
      if(endian==BIG_ENDIAN)
      {
         bswap_double(file_head.nnoise*file_head.tvals*file_head.ncorr*2,
                      data.corr);
         bswap_int(1,&(data.nc));
      }
      iw=fwrite(&(data.nc),sizeof(int),1,fdat);
      for (icorr=0;icorr<file_head.ncorr;icorr++)
      {
         chunk=file_head.nnoise*file_head.tvals*(2-file_head.isreal[icorr]);
         nw+=chunk;
         if (file_head.isreal[icorr])
         {
            for (i=0;i<chunk;i++)
               iw+=fwrite(&(data.corr[icorr*file_head.tvals*file_head.nnoise+i]),
                       sizeof(double),1,fdat);
         }else
         {
            iw+=fwrite(&(data.corr[icorr*file_head.tvals*file_head.nnoise]),
                       sizeof(double),chunk,fdat);
         }
      }
      if(endian==BIG_ENDIAN)
      {
         bswap_double(file_head.nnoise*file_head.tvals*file_head.ncorr*2,
                      data.corr);
         bswap_int(1,&(data.nc));
      }
      error_root(iw!=nw,1,"write_data [mesons.c]",
                 "Incorrect write count");
      fclose(fdat);
   }
}

static int read_data(void)
{
   int ir;
   int nr;
   int chunk;
   int icorr,i;
   double zero;

   zero=0;
   if(endian==BIG_ENDIAN)
      bswap_double(1,&zero);
   nr=1;
   ir=fread(&(data.nc),sizeof(int),1,fdat);

   for (icorr=0;icorr<file_head.ncorr;icorr++)
   {
      chunk=file_head.nnoise*file_head.tvals*(2-file_head.isreal[icorr]);
      nr+=chunk;
      if (file_head.isreal[icorr])
      {
         for (i=0;i<chunk;i++)
         {
            ir+=fread(&(data.corr[icorr*file_head.tvals*file_head.nnoise+i]),
                    sizeof(double),1,fdat);
            data.corr[icorr*file_head.tvals*file_head.nnoise+i].im=zero;
         }
      }else
      {
         ir+=fread(&(data.corr[icorr*file_head.tvals*file_head.nnoise]),
                    sizeof(double),chunk,fdat);
      }
   }

   if (ir==0)
      return 0;

   error_root(ir!=nr,1,"read_data [mesons.c]",
                 "Read error or incomplete data record");
   if(endian==BIG_ENDIAN)
   {
      bswap_double(nr,data.corr);
      bswap_int(1,&(data.nc));
   }
   return 1;
}

static void read_dirs(void)
{
   if (my_rank==0)
   {
      find_section("Run name");
      read_line("name","%s",nbase);
      read_line_opt("output",nbase,"%s",outbase);

      find_section("Directories");
      read_line("log_dir","%s",log_dir);

      if (noexp)
      {
         read_line("loc_dir","%s",loc_dir);
         cnfg_dir[0]='\0';
      }
      else
      {
         read_line("cnfg_dir","%s",cnfg_dir);
         loc_dir[0]='\0';
      }

      read_line("dat_dir","%s",dat_dir);

      find_section("Configurations");
      read_line("first","%d",&first);
      read_line("last","%d",&last);
      read_line("step","%d",&step);

      find_section("Random number generator");
      read_line("level","%d",&level);
      read_line("seed","%d",&seed);

      error_root((last<first)||(step<1)||(((last-first)%step)!=0),1,
                 "read_dirs [mesons.c]","Improper configuration range");
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(outbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(dat_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   if (noexp)
      error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,last,NPROC-1)>=NAME_SIZE,
                 1,"setup_files [mesons.c]","loc_dir name is too long");
   else
      error_root(name_size("%s/%sn%d",cnfg_dir,nbase,last)>=NAME_SIZE,
                 1,"setup_files [mesons.c]","cnfg_dir name is too long");

   check_dir_root(dat_dir);
   error_root(name_size("%s/%s.mesons.dat~",dat_dir,outbase)>=NAME_SIZE,
              1,"setup_files [mesons.c]","dat_dir name is too long");

   check_dir_root(log_dir);
   error_root(name_size("%s/%s.mesons.log~",log_dir,outbase)>=NAME_SIZE,
              1,"setup_files [mesons.c]","log_dir name is too long");

   sprintf(log_file,"%s/%s.mesons.log",log_dir,outbase);
   sprintf(end_file,"%s/%s.mesons.end",log_dir,outbase);
   sprintf(par_file,"%s/%s.mesons.par",dat_dir,outbase);
   sprintf(dat_file,"%s/%s.mesons.dat",dat_dir,outbase);
   sprintf(rng_file,"%s/%s.mesons.rng",dat_dir,outbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(dat_save,"%s~",dat_file);
   sprintf(rng_save,"%s~",rng_file);
}


static void read_lat_parms(void)
{
   double csw,cF;
   char tmpstring[NAME_SIZE];
   char tmpstring2[NAME_SIZE];
   int iprop,icorr;

   if (my_rank==0)
   {
      find_section("Measurements");
      read_line("nprop","%d",&nprop);
      read_line("ncorr","%d",&ncorr);
      read_line("nnoise","%d",&nnoise);
      read_line("noise_type","%s",tmpstring);
      read_line("csw","%lf",&csw);
      read_line("cF","%lf",&cF);

      error_root(nprop<1,1,"read_lat_parms [mesons.c]",
                 "Specified nprop must be larger than zero");
      error_root(ncorr<1,1,"read_lat_parms [mesons.c]",
                 "Specified ncorr must be larger than zero");
      error_root(nnoise<1,1,"read_lat_parms [mesons.c]",
                 "Specified nnoise must be larger than zero");

      noisetype=-1;
      if(strcmp(tmpstring,"Z2")==0)
         noisetype=Z2_NOISE;
      if(strcmp(tmpstring,"GAUSS")==0)
         noisetype=GAUSS_NOISE;
      if(strcmp(tmpstring,"U1")==0)
         noisetype=U1_NOISE;
      error_root(noisetype==-1,1,"read_lat_parms [mesons.c]",
                 "Unknown noise type");
   }
   MPI_Bcast(&nprop,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncorr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nnoise,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noisetype,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   kappas=malloc(nprop*sizeof(double));
   mus=malloc(nprop*sizeof(double));
   isps=malloc(nprop*sizeof(int));
   props1=malloc(ncorr*sizeof(int));
   props2=malloc(ncorr*sizeof(int));
   type1=malloc(ncorr*sizeof(int));
   type2=malloc(ncorr*sizeof(int));
   x0s=malloc(ncorr*sizeof(int));
   file_head.kappa1=malloc(ncorr*sizeof(double));
   file_head.kappa2=malloc(ncorr*sizeof(double));
   file_head.type1=type1;
   file_head.type2=type2;
   file_head.x0=x0s;
   file_head.isreal=malloc(ncorr*sizeof(int));

   error((kappas==NULL)||(mus==NULL)||(isps==NULL)||(props1==NULL)||
         (props2==NULL)||(type1==NULL)||(type2==NULL)||(x0s==NULL)||
         (file_head.kappa1==NULL)||(file_head.kappa2==NULL)||
         (file_head.isreal==NULL),
         1,"read_lat_parms [mesons.c]","Out of memory");

   if (my_rank==0)
   {
      for(iprop=0; iprop<nprop; iprop++)
      {
         sprintf(tmpstring,"Propagator %i",iprop);
         find_section(tmpstring);
         read_line("kappa","%lf",&kappas[iprop]);
         read_line("isp","%d",&isps[iprop]);
         /*TODO: read optional mu value*/
         mus[iprop]=0;
      }
      for(icorr=0; icorr<ncorr; icorr++)
      {
         sprintf(tmpstring,"Correlator %i",icorr);
         find_section(tmpstring);

         read_line("iprop","%d %d",&props1[icorr],&props2[icorr]);
         error_root((props1[icorr]<0)||(props1[icorr]>=nprop),1,"read_lat_parms [mesons.c]",
                 "Propagator index out of range");
         error_root((props2[icorr]<0)||(props2[icorr]>=nprop),1,"read_lat_parms [mesons.c]",
                 "Propagator index out of range");

         read_line("type","%s %s",tmpstring,tmpstring2);
         type1[icorr]=-1;
         type2[icorr]=-1;
         
         if(strncmp(tmpstring,"1",1)==0)
            type1[icorr]=ONE_TYPE;
         else if(strncmp(tmpstring,"G0G1",4)==0)
            type1[icorr]=GAMMA0GAMMA1_TYPE;
         else if(strncmp(tmpstring,"G0G2",4)==0)
            type1[icorr]=GAMMA0GAMMA2_TYPE;
         else if(strncmp(tmpstring,"G0G3",4)==0)
            type1[icorr]=GAMMA0GAMMA3_TYPE;
         else if(strncmp(tmpstring,"G0G5",4)==0)
            type1[icorr]=GAMMA0GAMMA5_TYPE;
         else if(strncmp(tmpstring,"G1G2",4)==0)
            type1[icorr]=GAMMA1GAMMA2_TYPE;
         else if(strncmp(tmpstring,"G1G3",4)==0)
            type1[icorr]=GAMMA1GAMMA3_TYPE;
         else if(strncmp(tmpstring,"G1G5",4)==0)
            type1[icorr]=GAMMA1GAMMA5_TYPE;
         else if(strncmp(tmpstring,"G2G3",4)==0)
            type1[icorr]=GAMMA2GAMMA3_TYPE;
         else if(strncmp(tmpstring,"G2G5",4)==0)
            type1[icorr]=GAMMA2GAMMA5_TYPE;
         else if(strncmp(tmpstring,"G3G5",4)==0)
            type1[icorr]=GAMMA3GAMMA5_TYPE;
         else if(strncmp(tmpstring,"G0",2)==0)
            type1[icorr]=GAMMA0_TYPE;
         else if(strncmp(tmpstring,"G1",2)==0)
            type1[icorr]=GAMMA1_TYPE;
         else if(strncmp(tmpstring,"G2",2)==0)
            type1[icorr]=GAMMA2_TYPE;
         else if(strncmp(tmpstring,"G3",2)==0)
            type1[icorr]=GAMMA3_TYPE;
         else if(strncmp(tmpstring,"G5",2)==0)
            type1[icorr]=GAMMA5_TYPE;
         
         if(strncmp(tmpstring2,"1",1)==0)
            type2[icorr]=ONE_TYPE;
         else if(strncmp(tmpstring2,"G0G1",4)==0)
            type2[icorr]=GAMMA0GAMMA1_TYPE;
         else if(strncmp(tmpstring2,"G0G2",4)==0)
            type2[icorr]=GAMMA0GAMMA2_TYPE;
         else if(strncmp(tmpstring2,"G0G3",4)==0)
            type2[icorr]=GAMMA0GAMMA3_TYPE;
         else if(strncmp(tmpstring2,"G0G5",4)==0)
            type2[icorr]=GAMMA0GAMMA5_TYPE;
         else if(strncmp(tmpstring2,"G1G2",4)==0)
            type2[icorr]=GAMMA1GAMMA2_TYPE;
         else if(strncmp(tmpstring2,"G1G3",4)==0)
            type2[icorr]=GAMMA1GAMMA3_TYPE;
         else if(strncmp(tmpstring2,"G1G5",4)==0)
            type2[icorr]=GAMMA1GAMMA5_TYPE;
         else if(strncmp(tmpstring2,"G2G3",4)==0)
            type2[icorr]=GAMMA2GAMMA3_TYPE;
         else if(strncmp(tmpstring2,"G2G5",4)==0)
            type2[icorr]=GAMMA2GAMMA5_TYPE;
         else if(strncmp(tmpstring2,"G3G5",4)==0)
            type2[icorr]=GAMMA3GAMMA5_TYPE;
         else if(strncmp(tmpstring2,"G0",2)==0)
            type2[icorr]=GAMMA0_TYPE;
         else if(strncmp(tmpstring2,"G1",2)==0)
            type2[icorr]=GAMMA1_TYPE;
         else if(strncmp(tmpstring2,"G2",2)==0)
            type2[icorr]=GAMMA2_TYPE;
         else if(strncmp(tmpstring2,"G3",2)==0)
            type2[icorr]=GAMMA3_TYPE;
         else if(strncmp(tmpstring2,"G5",2)==0)
            type2[icorr]=GAMMA5_TYPE;

         error_root((type1[icorr]==-1)||(type2[icorr]==-1),1,"read_lat_parms [mesons.c]",
                 "Unknown or unsupported Dirac structure");

         read_line("x0","%d",&x0s[icorr]);
         error_root((x0s[icorr]<=0)||(x0s[icorr]>=(NPROC0*L0-1)),1,"read_lat_parms [mesons.c]",
                 "Specified time x0 is out of range");
      }
   }

   MPI_Bcast(kappas,nprop,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(mus,nprop,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(isps,nprop,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(props1,ncorr,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(props2,ncorr,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(type1,ncorr,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(type2,ncorr,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(x0s,ncorr,MPI_INT,0,MPI_COMM_WORLD);

   set_lat_parms(0.0,1.0,kappas[0],0.0,0.0,csw,1.0,cF);
   set_sw_parms(sea_quark_mass(0));

   file_head.ncorr = ncorr;
   file_head.nnoise = nnoise;
   file_head.tvals = NPROC0*L0;
   tvals = NPROC0*L0;
   file_head.noisetype = noisetype;
   for(icorr=0; icorr<ncorr; icorr++)
   {
      file_head.kappa1[icorr]=kappas[props1[icorr]];
      file_head.kappa2[icorr]=kappas[props2[icorr]];
      if ((type1[icorr]==GAMMA5_TYPE)&&(type2[icorr]==GAMMA5_TYPE)&&
          (props1[icorr]==props2[icorr]))
         file_head.isreal[icorr]=1;
      else
         file_head.isreal[icorr]=0;
   }
   if (append)
      check_lat_parms(fdat);
   else
      write_lat_parms(fdat);
}


static void read_sap_parms(void)
{
   int bs[4];

   if (my_rank==0)
   {
      find_section("SAP");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   set_sap_parms(bs,1,4,5);
}


static void read_dfl_parms(void)
{
   int bs[4],Ns;
   int ninv,nmr,ncy,nkv,nmx;
   double kappa,mudfl,res;

   if (my_rank==0)
   {
      find_section("Deflation subspace");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
      read_line("Ns","%d",&Ns);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_parms(bs,Ns);

   if (my_rank==0)
   {
      find_section("Deflation subspace generation");
      read_line("kappa","%lf",&kappa);
      read_line("mu","%lf",&mudfl);
      read_line("ninv","%d",&ninv);
      read_line("nmr","%d",&nmr);
      read_line("ncy","%d",&ncy);
   }

   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mudfl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&ninv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_gen_parms(kappa,mudfl,ninv,nmr,ncy);

   if (my_rank==0)
   {
      find_section("Deflation projection");
      read_line("nkv","%d",&nkv);
      read_line("nmx","%d",&nmx);
      read_line("res","%lf",&res);
   }

   MPI_Bcast(&nkv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   set_dfl_pro_parms(nkv,nmx,res);
}


static void read_solvers(void)
{
   solver_parms_t sp;
   int i,j;
   int isap=0,idfl=0;

   for (i=0;i<nprop;i++)
   {
      j=isps[i];
      sp=solver_parms(j);
      if (sp.solver==SOLVERS)
      {
         read_solver_parms(j);
         sp=solver_parms(j);
         if (sp.solver==SAP_GCR)
            isap=1;
         if (sp.solver==DFL_SAP_GCR)
         {
            isap=1;
            idfl=1;
         }
      }
   }

   if (isap)
      read_sap_parms();

   if (idfl)
      read_dfl_parms();
}


static void read_infile(int argc,char *argv[])
{
   int ifile;

   if (my_rank==0)
   {
      flog=freopen("STARTUP_ERROR","w",stdout);
 
      ifile=find_opt(argc,argv,"-i");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [mesons.c]",
                 "Syntax: mesons -i <input file> [-noexp] [-a [-norng]]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [mesons.c]",
                 "Machine has unknown endianness");

      noexp=find_opt(argc,argv,"-noexp");
      append=find_opt(argc,argv,"-a");
      norng=find_opt(argc,argv,"-norng");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [mesons.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&norng,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();

   if (my_rank==0)
   {
      if (append)
         fdat=fopen(par_file,"rb");
      else
         fdat=fopen(par_file,"wb");

      error_root(fdat==NULL,1,"read_infile [mesons.c]",
                 "Unable to open parameter file");
   }
   read_lat_parms();
   read_solvers();

   if (my_rank==0)
   {
      fclose(fin);
      fclose(fdat);

      if (append==0)
         copy_file(par_file,par_save);
   }
}


static void check_old_log(int *fst,int *lst,int *stp)
{
   int ie,ic,isv;
   int fc,lc,dc,pc;
   int np[4],bp[4];

   fend=fopen(log_file,"r");
   error_root(fend==NULL,1,"check_old_log [mesons.c]",
              "Unable to open log file");

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;
   isv=0;

   while (fgets(line,NAME_SIZE,fend)!=NULL)
   {
      if (strstr(line,"process grid")!=NULL)
      {
         if (sscanf(line,"%dx%dx%dx%d process grid, %dx%dx%dx%d",
                    np,np+1,np+2,np+3,bp,bp+1,bp+2,bp+3)==8)
         {
            ipgrd[0]=((np[0]!=NPROC0)||(np[1]!=NPROC1)||
                      (np[2]!=NPROC2)||(np[3]!=NPROC3));
            ipgrd[1]=((bp[0]!=NPROC0_BLK)||(bp[1]!=NPROC1_BLK)||
                      (bp[2]!=NPROC2_BLK)||(bp[3]!=NPROC3_BLK));
         }
         else
            ie|=0x1;
      }
      
      if (strstr(line,"fully processed")!=NULL)
      {
         pc=lc;

         if (sscanf(line,"Configuration no %d",&lc)==1)
         {
            ic+=1;
            isv=1;
         }
         else
            ie|=0x1;

         if (ic==1)
            fc=lc;
         else if (ic==2)
            dc=lc-fc;
         else if ((ic>2)&&(lc!=(pc+dc)))
            ie|=0x2;
      }
      else if (strstr(line,"Configuration no")!=NULL)
         isv=0;
   }

   fclose(fend);

   error_root((ie&0x1)!=0x0,1,"check_old_log [mesons.c]",
              "Incorrect read count");
   error_root((ie&0x2)!=0x0,1,"check_old_log [mesons.c]",
              "Configuration numbers are not equally spaced");
   error_root(isv==0,1,"check_old_log [mesons.c]",
              "Log file extends beyond the last configuration save");

   (*fst)=fc;
   (*lst)=lc;
   (*stp)=dc;
}


static void check_old_dat(int fst,int lst,int stp)
{
   int ie,ic;
   int fc,lc,dc,pc;

   fdat=fopen(dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [mesons.c]",
              "Unable to open data file");

   check_file_head();

   fc=0;
   ic=0;
   lc=0;
   dc=0;
   pc=0;
   ie=0x0;

   while (read_data()==1)
   {
      pc=lc;
      lc=data.nc;
      ic+=1;

      if (ic==1)
         fc=lc;
      else if (ic==2)
         dc=lc-fc;
      else if ((ic>2)&&(lc!=(pc+dc)))
         ie|=0x1;
   }

   fclose(fdat);

   error_root(ic==0,1,"check_old_dat [mesons.c]",
              "No data records found");
   error_root((ie&0x1)!=0x0,1,"check_old_dat [mesons.c]",
              "Configuration numbers are not equally spaced");
   error_root((fst!=fc)||(lst!=lc)||(stp!=dc),1,"check_old_dat [mesons.c]",
              "Configuration range is not as reported in the log file");
}


static void check_files(void)
{
   int fst,lst,stp;

   ipgrd[0]=0;
   ipgrd[1]=0;
   
   if (my_rank==0)
   {
      if (append)
      {
         check_old_log(&fst,&lst,&stp);
         check_old_dat(fst,lst,stp);

         error_root((fst!=lst)&&(stp!=step),1,"check_files [mesons.c]",
                    "Continuation run:\n"
                    "Previous run had a different configuration separation");
         error_root(first!=lst+step,1,"check_files [mesons.c]",
                    "Continuation run:\n"
                    "Configuration range does not continue the previous one");
      }
      else
      {
         fin=fopen(log_file,"r");
         error_root(fin!=NULL,1,"check_files [mesons.c]",
                    "Attempt to overwrite old *.log file");
         fdat=fopen(dat_file,"r");
         error_root(fdat!=NULL,1,"check_files [mesons.c]",
                    "Attempt to overwrite old *.dat file");
         fdat=fopen(dat_file,"wb");
         error_root(fdat==NULL,1,"check_files [mesons.c]",
                    "Unable to open data file");
         write_file_head();
         fclose(fdat);
      }
   }
}


static void print_info(void)
{
   int i,isap,idfl;
   long ip;
   lat_parms_t lat;

   if (my_rank==0)
   {
      ip=ftell(flog);
      fclose(flog);

      if (ip==0L)
         remove("STARTUP_ERROR");

      if (append)
         flog=freopen(log_file,"a",stdout);
      else
         flog=freopen(log_file,"w",stdout);

      error_root(flog==NULL,1,"print_info [mesons.c]",
                 "Unable to open log file");
      printf("\n");

      if (append)
         printf("Continuation run\n\n");
      else
      {
         printf("Computation of meson correlators\n");
         printf("--------------------------------\n\n");
         printf("cnfg   base name: %s\n",nbase);
         printf("output base name: %s\n\n",outbase);
      }

      printf("openQCD version: %s, meson version: %s\n",openQCD_RELEASE,
                                                      mesons_RELEASE);
      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
      if (noexp)
         printf("Configurations are read in imported file format\n\n");
      else
         printf("Configurations are read in exported file format\n\n");

      if ((ipgrd[0]!=0)&&(ipgrd[1]!=0))
         printf("Process grid and process block size changed:\n");
      else if (ipgrd[0]!=0)
         printf("Process grid changed:\n");
      else if (ipgrd[1]!=0)
         printf("Process block size changed:\n");

      if ((append==0)||(ipgrd[0]!=0)||(ipgrd[1]!=0))
      {
         printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
         printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
         printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
         printf("%dx%dx%dx%d process block size\n",
                NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);

         if (append)
            printf("\n");
         else
            printf("SF boundary conditions on the quark fields\n\n");
      }
      
      
      if (append)
      {
         printf("Random number generator:\n");

         if (norng)
            printf("level = %d, seed = %d, effective seed = %d\n\n",
                   level,seed,seed^(first-step));
         else
         {
            printf("State of ranlxs and ranlxd reset to the\n");
            printf("last exported state\n\n");
         }
      }
      else
      {
         printf("Random number generator:\n");
         printf("level = %d, seed = %d\n\n",level,seed);

         lat=lat_parms();

         printf("Measurements:\n");
         printf("nprop     = %i\n",nprop);
         printf("ncorr     = %i\n",ncorr);
         printf("nnoise    = %i\n",nnoise);
         if (noisetype==Z2_NOISE)
            printf("noisetype = Z2\n");
         if (noisetype==GAUSS_NOISE)
            printf("noisetype = GAUSS\n");
         if (noisetype==U1_NOISE)
            printf("noisetype = U1\n");
         printf("csw       = %.6f\n",lat.csw);
         printf("cF        = %.6f\n\n",lat.cF);

         for (i=0; i<nprop; i++)
         {
            printf("Propagator %i:\n",i);
            printf("kappa  = %.6f\n",kappas[i]);
            printf("isp    = %i\n",isps[i]);
            printf("mu     = %.6f\n\n",mus[i]);
         }
         for (i=0; i<ncorr; i++)
         {
            printf("Correlator %i:\n",i);
            printf("iprop  = %i %i\n",props1[i],props2[i]);
            printf("type   = %i %i\n",type1[i],type2[i]); /*TODO: strings*/
            printf("x0     = %i\n\n",x0s[i]);
         }
      }
      print_solver_parms(&isap,&idfl);

      if (isap)
         print_sap_parms(0);

      if (idfl)
         print_dfl_parms(0);

      printf("Configurations no %d -> %d in steps of %d\n\n",
             first,last,step);
      fflush(flog);
   }
}


static void dfl_wsize(int *nws,int *nwv,int *nwvd)
{
   dfl_parms_t dp;
   dfl_pro_parms_t dpp;

   dp=dfl_parms();
   dpp=dfl_pro_parms();

   MAX(*nws,dp.Ns+2);
   MAX(*nwv,2*dpp.nkv+2);
   MAX(*nwvd,4);
}


static void make_proplist(void)
{
   int i,j,k,iprop,icorr;
   char *kappatype;

   proplist.nux0=0;
   proplist.ux0=malloc(NPROC0*L0*sizeof(int));
   error(proplist.ux0==NULL,1,"make_proplist [mesons.c]","Out of memory");

   /* unique x0 values */
   for (icorr=0;icorr<ncorr;icorr++)
   {
      for (j=0;j<proplist.nux0;j++)
      {
         if (proplist.ux0[j]==x0s[icorr])
            break;
      }
      if (j==proplist.nux0)
      {
         proplist.ux0[j]=x0s[icorr];
         proplist.nux0++;
      }
   }

   proplist.nprop=malloc(proplist.nux0*sizeof(int));
   proplist.prop=malloc(proplist.nux0*sizeof(int*));
   proplist.type=malloc(proplist.nux0*sizeof(int*));
   kappatype=malloc(MAX_TYPE*nprop*sizeof(char));
   error((proplist.nprop==NULL)||(proplist.prop==NULL)||(proplist.type==NULL)
               ||(kappatype==NULL),
               1,"make_proplist [mesons.c]","Out of memory");
   proplist.nmax=0;

   for (i=0;i<proplist.nux0;i++)
   {
      for (j=0;j<MAX_TYPE*nprop;j++)
         kappatype[j]=0;

      proplist.nprop[i]=0;
      for (icorr=0;icorr<ncorr;icorr++)
      {
         if (x0s[icorr]==proplist.ux0[i])
         {
            if (!kappatype[type1[icorr]+MAX_TYPE*props2[icorr]])
            {
               kappatype[type1[icorr]+MAX_TYPE*props2[icorr]]=1;
               proplist.nprop[i]++;
            }
            if (!kappatype[GAMMA5_TYPE+MAX_TYPE*props1[icorr]])
            {
               kappatype[GAMMA5_TYPE+MAX_TYPE*props1[icorr]]=1;
               proplist.nprop[i]++;
            }
         }
      }

      if (proplist.nprop[i]>proplist.nmax)
         proplist.nmax=proplist.nprop[i];

      proplist.prop[i]=malloc(proplist.nprop[i]*sizeof(int));
      proplist.type[i]=malloc(proplist.nprop[i]*sizeof(int));
      error((proplist.prop[i]==NULL)||(proplist.type[i]==NULL),
                 1,"make_proplist [mesons.c]","Out of memory");
      j=0;
      for (k=0;k<MAX_TYPE;k++)
      {
         for (iprop=0;iprop<nprop;iprop++)
         {
            if (kappatype[k+MAX_TYPE*iprop])
            {
               proplist.prop[i][j]=iprop;
               proplist.type[i][j]=k;
               j++;
            }
         }
      }
   }
   free(kappatype);
}

static void wsize(int *nws,int *nwsd,int *nwv,int *nwvd)
{
   int nsd;
   solver_parms_t sp;

   (*nws)=0;
   (*nwsd)=0;
   (*nwv)=0;
   (*nwvd)=0;

   sp=solver_parms(0);
   nsd=proplist.nmax+2;

   if (sp.solver==CGNE)
   {
      MAX(*nws,5);
      MAX(*nwsd,nsd+3);
   }
   else if (sp.solver==SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+1);
      MAX(*nwsd,nsd+2);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+2);
      MAX(*nwsd,nsd+4);
      dfl_wsize(nws,nwv,nwvd);
   }
   else
      error_root(1,1,"wsize [mesons.c]",
                 "Unknown or unsupported solver");
}



static void random_source(spinor_dble *eta, int x0)
{
   int y0,iy,ix;

   set_sd2zero(VOLUME,eta);
   y0=x0-cpr[0]*L0;

   if ((y0>=0)&&(y0<L0))
   {
      if (noisetype==Z2_NOISE)
      {
         for (iy=0;iy<(L1*L2*L3);iy++)
         {
            ix=ipt[iy+y0*L1*L2*L3];
            random_Z2_sd(1,eta+ix);
         }
      }
      else if (noisetype==GAUSS_NOISE)
      {
         for (iy=0;iy<(L1*L2*L3);iy++)
         {
            ix=ipt[iy+y0*L1*L2*L3];
            random_sd(1,eta+ix,1.0);
         }
      }
      else if (noisetype==U1_NOISE)
      {
         for (iy=0;iy<(L1*L2*L3);iy++)
         {
            ix=ipt[iy+y0*L1*L2*L3];
            random_U1_sd(1,eta+ix);
         }
      }
   }
}


static void solve_dirac(int prop, spinor_dble *eta, spinor_dble *psi,
                        int *status)
{
   solver_parms_t sp;
   sap_parms_t sap;

   sp=solver_parms(isps[prop]);
   set_sw_parms(0.5/kappas[prop]-4.0);

   if (sp.solver==CGNE)
   {
      mulg5_dble(VOLUME,eta);

      tmcg(sp.nmx,sp.res,mus[prop],eta,eta,status);
      if (my_rank==0)
         printf("%i\n",status[0]);
      error_root(status[0]<0,1,"solve_dirac [mesons.c]",
                 "CGNE solver failed (status = %d)",status[0]);

      Dw_dble(-mus[prop],eta,psi);
      mulg5_dble(VOLUME,psi);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      sap_gcr(sp.nkv,sp.nmx,sp.res,mus[prop],eta,psi,status);
      if (my_rank==0)
         printf("%i\n",status[0]);
      error_root(status[0]<0,1,"solve_dirac [mesons.c]",
                 "SAP_GCR solver failed (status = %d)",status[0]);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mus[prop],eta,psi,status);
      if (my_rank==0)
         printf("%i %i\n",status[0],status[1]);
      error_root((status[0]<0)||(status[1]<0),1,
                 "solve_dirac [mesons.c]","DFL_SAP_GCR solver failed "
                 "(status = %d,%d)",status[0],status[1]);
   }
   else
      error_root(1,1,"solve_dirac [mesons.c]",
                 "Unknown or unsupported solver");
}


/* xi = \gamma_5 Gamma^\dagger eta */
void make_source(spinor_dble *eta, int type, spinor_dble *xi)
{
   switch (type)
   {
      case GAMMA0_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg0g5_dble(VOLUME,xi);
         break;
      case GAMMA1_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg1g5_dble(VOLUME,xi);
         break;
      case GAMMA2_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg2g5_dble(VOLUME,xi);
         break;
      case GAMMA3_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg3g5_dble(VOLUME,xi);
         break;
      case GAMMA5_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         break;
      case GAMMA0GAMMA1_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg2g3_dble(VOLUME,xi);
         break;
      case GAMMA0GAMMA2_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg1g3_dble(VOLUME,xi);
         break;
      case GAMMA0GAMMA3_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg1g2_dble(VOLUME,xi);
         break;
      case GAMMA0GAMMA5_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg0_dble(VOLUME,xi);
         break;
      case GAMMA1GAMMA2_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg0g3_dble(VOLUME,xi);
         break;
      case GAMMA1GAMMA3_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg0g2_dble(VOLUME,xi);
         break;
      case GAMMA1GAMMA5_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg1_dble(VOLUME,xi);
         break;
      case GAMMA2GAMMA3_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg0g1_dble(VOLUME,xi);
         break;
      case GAMMA2GAMMA5_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg2_dble(VOLUME,xi);
         break;
      case GAMMA3GAMMA5_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg3_dble(VOLUME,xi);
         break;
      case ONE_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg5_dble(VOLUME,xi);
         break;
      default:
         error_root(1,1,"make_source [mesons.c]",
                 "Unknown or unsupported type");
   }
}

void make_xi(spinor_dble *eta,int type,spinor_dble *xi)
{
   /* xi = -\bar Gamma^\dagger \gamma_5 eta */
   switch (type)
   {
      case GAMMA0_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg0g5_dble(VOLUME,xi);
         break;
      case GAMMA1_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg1g5_dble(VOLUME,xi);
         break;
      case GAMMA2_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg2g5_dble(VOLUME,xi);
         break;
      case GAMMA3_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg3g5_dble(VOLUME,xi);
         break;
      case GAMMA5_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         break;
      case GAMMA0GAMMA1_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg2g3_dble(VOLUME,xi);
         break;
      case GAMMA0GAMMA2_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg1g3_dble(VOLUME,xi);
         break;
      case GAMMA0GAMMA3_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg1g2_dble(VOLUME,xi);
         break;
      case GAMMA0GAMMA5_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg0_dble(VOLUME,xi);
         break;
      case GAMMA1GAMMA2_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg0g3_dble(VOLUME,xi);
         break;
      case GAMMA1GAMMA3_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg0g2_dble(VOLUME,xi);
         break;
      case GAMMA1GAMMA5_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg1_dble(VOLUME,xi);
         break;
      case GAMMA2GAMMA3_TYPE:
         assign_sd2sd(VOLUME,eta,xi);
         mulg0g1_dble(VOLUME,xi);
         break;
      case GAMMA2GAMMA5_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg2_dble(VOLUME,xi);
         break;
      case GAMMA3GAMMA5_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg3_dble(VOLUME,xi);
         break;
      case ONE_TYPE:
         assign_msd2sd(VOLUME,eta,xi);
         mulg5_dble(VOLUME,xi);
         break;
      default:
         error_root(1,1,"make_xi [mesons.c]",
                 "Unknown or unsupported type");
   }
}


static void correlators(void)
{
   int ix0,inoise,iprop,icorr,ip1,ip2,l,stat[4],y0,iy;
   spinor_dble *eta,*xi,**zeta,**wsd;
   complex_dble tmp;

   wsd=reserve_wsd(proplist.nmax+2);
   eta=wsd[0];
   xi=wsd[1];
   zeta=malloc(proplist.nmax*sizeof(spinor_dble*));
   error(zeta==NULL,1,"correlators [mesons.c]","Out of memory");

   for (l=0;l<proplist.nmax;l++)
      zeta[l]=wsd[l+2];

   for (l=0;l<nnoise*ncorr*tvals;l++)
   {
      data.corr_tmp[l].re=0.0;
      data.corr_tmp[l].im=0.0;
   }

   if (my_rank==0)
      printf("Inversions:\n");
   for (ix0=0;ix0<proplist.nux0;ix0++)
   {
      if (my_rank==0)
         printf("   x0=%i\n",proplist.ux0[ix0]);
      for (inoise=0;inoise<nnoise;inoise++)
      {
         if (my_rank==0)
            printf("      noise vector %i\n",inoise);
         random_source(eta,proplist.ux0[ix0]);
         for (iprop=0;iprop<proplist.nprop[ix0];iprop++)
         {
            if (my_rank==0)
               printf("         type=%i, prop=%i, status:",
                   proplist.type[ix0][iprop], proplist.prop[ix0][iprop]);
            make_source(eta,proplist.type[ix0][iprop],xi);
            solve_dirac(proplist.prop[ix0][iprop],xi,zeta[iprop],stat);
         }
         /* combine propagators to correlators */
         for (icorr=0;icorr<ncorr;icorr++)
         {
            if (x0s[icorr]==proplist.ux0[ix0])
            {
               /* find the two propagators that are needed for this icorr */
               ip1=0;
               ip2=0;
               for (iprop=0;iprop<proplist.nprop[ix0];iprop++)
               {
                  if ((type1[icorr]==proplist.type[ix0][iprop])&&
                      (props2[icorr]==proplist.prop[ix0][iprop]))
                     ip1=iprop;
                  if ((GAMMA5_TYPE==proplist.type[ix0][iprop])&&
                      (props1[icorr]==proplist.prop[ix0][iprop]))
                     ip2=iprop;
               }
               make_xi(zeta[ip1],type2[icorr],xi);
               for (y0=0;y0<L0;y0++)
               {
                  for (l=0;l<L1*L2*L3;l++)
                  {
                     iy = ipt[l+y0*L1*L2*L3];
                     tmp = spinor_prod_dble(1,0,xi+iy,zeta[ip2]+iy);
                     data.corr_tmp[inoise+nnoise*(cpr[0]*L0+y0
                        +file_head.tvals*icorr)].re += tmp.re;
                     data.corr_tmp[inoise+nnoise*(cpr[0]*L0+y0
                        +file_head.tvals*icorr)].im += tmp.im;
                  }
               }
            }
         }
      }
   }

   MPI_Allreduce(data.corr_tmp,data.corr,nnoise*ncorr*file_head.tvals*2
      ,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
   free(zeta);
   release_wsd();
}

static void set_data(int nc)
{
   data.nc=nc;
   correlators();

   if (my_rank==0)
   {
      printf("G(t) =  %.4e%+.4ei",data.corr[0].re,data.corr[0].im);
      printf(",%.4e%+.4ei,...",data.corr[1].re,data.corr[1].im);
      printf(",%.4e%+.4ei",data.corr[file_head.tvals-1].re,
                           data.corr[file_head.tvals-1].im);
      printf("\n");
      fflush(flog);
   }
}

static void init_rng(void)
{
   int ic;

   if (append)
   {
      if (norng)
         start_ranlux(level,seed^(first-step));
      else
      {
         ic=import_ranlux(rng_file);
         error_root(ic!=(first-step),1,"init_rng [mesons.c]",
                    "Configuration number mismatch (*.rng file)");
      }
   }
   else
      start_ranlux(level,seed);
}


static void save_ranlux(void)
{
   int nlxs,nlxd;

   if (rlxs_state==NULL)
   {
      nlxs=rlxs_size();
      nlxd=rlxd_size();

      rlxs_state=malloc((nlxs+nlxd)*sizeof(int));
      rlxd_state=rlxs_state+nlxs;

      error(rlxs_state==NULL,1,"save_ranlux [mesons.c]",
            "Unable to allocate state arrays");
   }

   rlxs_get(rlxs_state);
   rlxd_get(rlxd_state);
}


static void restore_ranlux(void)
{
   rlxs_reset(rlxs_state);
   rlxd_reset(rlxd_state);
}


static void check_endflag(int *iend)
{
   if (my_rank==0)
   {
      fend=fopen(end_file,"r");

      if (fend!=NULL)
      {
         fclose(fend);
         remove(end_file);
         (*iend)=1;
         printf("End flag set, run stopped\n\n");
      }
      else
         (*iend)=0;
   }

   MPI_Bcast(iend,1,MPI_INT,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int nc,iend,status[4];
   int nws,nwsd,nwv,nwvd;
   double wt1,wt2,wtavg;
   dfl_parms_t dfl;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   alloc_data();
   check_files();
   print_info();
   dfl=dfl_parms();

   geometry();
   init_rng();

   make_proplist();
   wsize(&nws,&nwsd,&nwv,&nwvd);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);

   iend=0;
   wtavg=0.0;

   for (nc=first;(iend==0)&&(nc<=last);nc+=step)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      if (my_rank==0)
         printf("Configuration no %d\n",nc);

      if (noexp)
      {
         save_ranlux();
         sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,nc,my_rank);
         read_cnfg(cnfg_file);
         restore_ranlux();
      }
      else
      {
         sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,nc);
         import_cnfg(cnfg_file);
      }

      if (dfl.Ns)
      {
         dfl_modes(status);
         error_root(status[0]<0,1,"main [mesons.c]",
                    "Deflation subspace generation failed (status = %d)",
                    status[0]);

         if (my_rank==0)
            printf("Deflation subspace generation: status = %d\n",status[0]);
      }

      set_data(nc);
      write_data();

      export_ranlux(nc,rng_file);
      error_chk();
      
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wtavg+=(wt2-wt1);
      

      if (my_rank==0)
      {
         printf("Configuration no %d fully processed in %.2e sec ",
                nc,wt2-wt1);
         printf("(average = %.2e sec)\n\n",
                wtavg/(double)((nc-first)/step+1));
      }
      check_endflag(&iend);

      if (my_rank==0)
      {
         fflush(flog);
         copy_file(log_file,log_save);
         copy_file(dat_file,dat_save);
         copy_file(rng_file,rng_save);
      }
   }


   if (my_rank==0)
   {
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
