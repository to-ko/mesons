/*******************************************************************************
 * utils_mutils.c
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
 * modules/utils/mutils.c
 *
 *
 ******************************************************************************/

#define UTILS_MUTILS_C

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "global.h"

static char line[NAME_SIZE+1];

/* this is an exact copy of /modules/utils/mutils.c's version */
static char *get_line(void)
{
   char *s,*c;

   s=fgets(line,NAME_SIZE+1,stdin);

   if (s!=NULL)
   {
      error_root(strlen(line)==NAME_SIZE,1,"get_line [mutils.c]",
                 "Input line is longer than NAME_SIZE-1");

      c=strchr(line,'#');
      if (c!=NULL)
         c[0]='\0';
   }

   return s;
}

/* this is an exact copy of /modules/utils/mutils.c's version */
static void check_tag(char *tag)
{
   if (tag[0]=='\0')
      return;

   error_root((strspn(tag," 0123456789.")!=0L)||
              (strcspn(tag," \n")!=strlen(tag)),1,
              "check_tag [mutils.c]","Improper tag %s",tag);
}


static long find_tag_opt(char *tag)
{
   int ie;
   long tofs,lofs,ofs;
   char *s,*pl,*pr;

   ie=0;
   tofs=-1L;
   lofs=ftell(stdin);
   rewind(stdin);
   ofs=ftell(stdin);
   s=get_line();

   while (s!=NULL)
   {
      pl=strchr(line,'[');
      pr=strchr(line,']');

      if ((pl==(line+strspn(line," \t")))&&(pr>pl))
      {
         if (ofs<lofs)
         {
            ie=0;
            tofs=-1L;
         }
         else
            break;
      }
      else
      {
         pl=line+strspn(line," \t");
         pr=pl+strcspn(pl," \t\n");
         pr[0]='\0';

         if (strcmp(pl,tag)==0)
         {
            if (tofs!=-1L)
               ie=1;
            tofs=ofs;
         }
      }

      ofs=ftell(stdin);
      s=get_line();
   }

   if (tofs==-1L)
      return tofs;
   else
   {
      error_root(ie!=0,1,"find_tag_opt [mutils.c]",
              "Tag %s occurs more than once in the current section",tag);

      ie=fseek(stdin,tofs,SEEK_SET);
      error_root(ie!=0,1,"find_tag_opt [mutils.c]",
              "Unable to go to line with tag %s",tag);
   }
   return tofs;
}


/*     long read_line_opt(char *tag,char *defaultline,char *format,...)
 *
 *     On process 0, this program reads a line of text and data from stdin
 *     in a controlled manner. The tag cannot
 *     be the empty string "", the program searches
 *     for the tag in the current section. If the tag is not
 *     found, default values are assigned to the variables. The program returns
 *     the offset of the line from the beginning
 *     of the file and positions the file pointer to the next line. On
 *     processes other than 0, the program does nothing and returns -1L.
 *
 */
long read_line_opt(char *tag,char *defaultline,char *format, ...)
{
   int my_rank,is,ic;
   long tofs;
   char *pl,*p;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      check_tag(tag);

      error_root(tag[0]=='\0',1,"read_line_opt [mutils.c]",
              "optional parameters must have a non empty tag");

      tofs=find_tag_opt(tag);
      if(tofs != -1L)
      {
         get_line();
         pl=line+strspn(line," \t");
         pl+=strcspn(pl," \t\n");
      }else
      {
         /* optional parameter gets default value */
         error_root(strlen(defaultline)+strlen(tag)+1>NAME_SIZE,1,
             "read_line_opt [mutils.c]",
             "default line too long");
         sprintf(line,"%s %s",tag,defaultline);
         pl=line+strspn(line," \t");
         pl+=strcspn(pl," \t\n");
      }

      va_start(args,format);

      for (p=format;;)
      {
         p+=strspn(p," ");
         ic=0;
         is=2;

         if ((p[0]=='\0')||(p[0]=='\n'))
            break;
         else if (p==strstr(p,"%s"))
            ic=sscanf(pl,"%s",va_arg(args,char*));
         else if (p==strstr(p,"%d"))
            ic=sscanf(pl,"%d",va_arg(args,int*));
         else if (p==strstr(p,"%f"))
            ic=sscanf(pl,"%f",va_arg(args,float*));
         else if (p==strstr(p,"%lf"))
         {
            is=3;
            ic=sscanf(pl,"%lf",va_arg(args,double*));
         }
         else
            error_root(1,1,"read_line_opt [mutils.c]",
                       "Incorrect format string %s on line with tag %s",
                       format,tag);

         error_root(ic!=1,1,"read_line_opt [mutils.c]",
                    "Missing data item(s) on line with tag %s",tag);

         p+=is;
         pl+=strspn(pl," \t");
         pl+=strcspn(pl," \t\n");
      }

      va_end(args);

      return tofs;
   }
   else
      return -1L;
}
