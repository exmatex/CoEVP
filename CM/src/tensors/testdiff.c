/*
Program to compare two files
ignore white space
ignore "small" differences in numbers, like round-off errors and
typical differences on different machines like "numerical zeros"
*/

#include "stdio.h"
#include "math.h"
#include "string.h"
#include "ctype.h"

#define MAX_DIFFERENCE 1E-10
#define ZERO 1E-10
#define REL_PRECISION 1E-5

int  readgenfile(FILE *tmpfile1, FILE *tmpfile2);
int  compare(double *refcomp, double *rescomp);

int main(int argc, char *argv[])
{

  FILE *res, *ref;
  int  j;
  
  if (argc==1)
  {
     printf("Usage: testdiff reference_file file_to_check\n");
  }

  ref = fopen(argv[1], "r");
  if (ref==NULL)
  {
     printf("testdiff: Can't open reference file %s\n", argv[1]);
     return(1);
  }

  res = fopen(argv[2], "r");
  if (res==NULL)
  {
     printf("testdiff: Can't open result file %s\n", argv[2]);
     return(1);
  }

  j=readgenfile(ref, res);

  if (j==1)
     return(1);
  else
     return(0);

} /* main */      

int readgenfile(FILE *tmpfile1, FILE *tmpfile2) /* This function reads general files */
{
   int    line, rtnval, status, err_in_line;
   char   wstr1[160], wstr2[160];   /* lines read */
   char   *p1, *p2, *r1, *r2, *begin1, *begin2;
   double d1, d2;

   status=0;
   line = 0;

   do
   {
      line++;
      err_in_line = 0;

      /* reading the ref and res file */
      r1 = fgets(wstr1, 160, tmpfile1);
      r2 = fgets(wstr2, 160, tmpfile2);

      if ( r1==NULL ) wstr1[0] = '\0';
      if ( r2==NULL ) wstr2[0] = '\0';

      if ( r1==NULL && r2 == NULL ) break;

      p1 = wstr1;
      p2 = wstr2;

      while ( *p1 != '\0' || *p2 != '\0' )
      {
         /* go to next non space character */
         while ( isspace(*p1) && *p1 != '\0' ) p1++;
         while ( isspace(*p2) && *p2 != '\0' ) p2++;

         if ( *p1=='\0' && *p2=='\0' ) break;

         begin1 = p1;
         begin2 = p2;

         /* check whether next sub-strings match */

         do
         {
            if ( *p1=='\0' && *p2=='\0' ) break;

            if ( *p1=='\0' || *p2=='\0' )
            {
               err_in_line++;
               break;
            }

            if ( *p1 != *p2 )
            {
               /* there is some difference */

               /* check whether it are almost equal doubles */
               rtnval = sscanf( begin1, "%lf", &d1 );
               if ( rtnval==1 )
               {
                  rtnval = sscanf( begin2, "%lf", &d2 );
                  if ( rtnval==0 || !compare( &d1, &d2 ) )
                  {
                     err_in_line++;
                     break;
                  }

                  /* go to next separator character */
                  while ( !isspace(*p1) && *p1 != '\0' ) p1++;
                  while ( !isspace(*p2) && *p2 != '\0' ) p2++;
                  break;
               }
               else
               {
                  err_in_line++;
                  break;
               }
            }
            p1++;
            p2++;
         }
         while ( !( isspace(*p1) && isspace(*p2) ) );

         if ( err_in_line != 0 )
         {
            status = 1;
            printf( "%d\n<<<%s>>>%s\n", line, wstr1, wstr2 );
            break;
         }
      }

   } while ( r1!=NULL || r2!=NULL );

   return status;
} /* readgenfile */

int compare(double *refcomp, double *rescomp) /* This function filters and compares two numbers */
{
   double d, f, g, h;

   if (fabs(*refcomp)<ZERO) *refcomp=0.0;
   if (fabs(*rescomp)<ZERO) *rescomp=0.0;
   if (*refcomp==0.0 && *rescomp==0.0) return (1);

   d=fabs(*refcomp-*rescomp);
   if (d<MAX_DIFFERENCE)
   {
      *rescomp = *refcomp;
      return(1);
   }

   if (*refcomp!=0.0 && *rescomp!=0.0)
   {  f=*refcomp;
      g=*rescomp;
	  h=REL_PRECISION*pow(10,ceil(log10(fabs(g))));
      if (fabs(f-g)<h)
      {
         *rescomp = *refcomp;
         return(1);
      }
      else return (0);
   }
   else if (*refcomp==0.0 || *rescomp==0.0)
   { /* one of both is zero, but not both are zero */
      return (0);
   }

   return (1);
} /* compare */
