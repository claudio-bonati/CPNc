#ifndef ENDIANNESS_C
#define ENDIANNESS_C

#include<stdio.h>
#include<string.h>

#include"../include/endianness.h"

int endian(void)
   {
   int i = 1;
   char *p = (char *)&i;

   if (p[0] == 1)
        return 0; // LITTLE ENDIAN
   else
        return 1; // BIG ENDIAN
   }


void SwapBytesInt(void *pv)
    {
    char *p = pv;
    size_t lo, hi;
    char tmp;

    for(lo=0, hi=sizeof(int)-1; hi>lo; lo++, hi--)
       {
       tmp=p[lo];
       p[lo] = p[hi];
       p[hi] = tmp;
       }
    }

void SwapBytesFloat(void *pv)
    {
    char *p = pv;
    size_t lo, hi;
    char tmp;

    for(lo=0, hi=sizeof(float)-1; hi>lo; lo++, hi--)
       {
       tmp=p[lo];
       p[lo] = p[hi];
       p[hi] = tmp;
       }
    }

void SwapBytesDouble(void *pv)
    {
    char *p = pv;
    size_t lo, hi;
    char tmp;

    for(lo=0, hi=sizeof(double)-1; hi>lo; lo++, hi--)
       {
       tmp=p[lo];
       p[lo] = p[hi];
       p[hi] = tmp;
       }
    }


int print_on_binary_file_bigen_double(FILE *fp, double r)
  {
  size_t err;

  if(endian()==0) // little endian machine
    {
    SwapBytesDouble(&r);
    err=fwrite(&r, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problem in binary writing on file a double (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }
    }
  else
    {
    err=fwrite(&r, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problem in binary writing on file a double (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }
    }

  return 0;
  }


int read_from_binary_file_bigen_double(FILE *fp, double *r)
  {
  size_t err;

  if(endian()==0) // little endian machine
    {
    err=fread(r, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problems reading double from file (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }

    SwapBytesDouble(r);
    }
  else
    {
    err=fread(r, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problems reading double from file (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }
    }

  return 0;
  }


int print_on_binary_file_bigen_doublecomplex(FILE *fp, double complex r)
  {
  size_t err;
  double re, im;

  if(endian()==0) // little endian machine
    {
    re=creal(r);
    im=cimag(r);

    SwapBytesDouble(&re);
    SwapBytesDouble(&im);

    err=fwrite(&re, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problem in binary writing on file a double complex (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }

    err=fwrite(&im, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problem in binary writing on file a double complex (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }
    }
  else
    {
    re=creal(r);
    im=cimag(r);

    err=fwrite(&re, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problem in binary writing on file a double complex (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }

    err=fwrite(&im, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problem in binary writing on file a double complex (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }
    }

  return 0;
  }


int read_from_binary_file_bigen_doublecomplex(FILE *fp, double complex *r)
  {
  size_t err;
  double aux[2], re, im;

  if(endian()==0) // little endian machine
    {
    err=fread(&re, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problems reading double complex from file (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }
    err=fread(&im, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problems reading double complex from file (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }

    SwapBytesDouble(&re);
    SwapBytesDouble(&im);
    aux[0]=re;
    aux[1]=im;

    memcpy((void *)r, (void*)aux, sizeof(aux));
    //equivalent to r=re+im*I;
    }
  else
    {
    err=fread(&re, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problems reading double complex from file (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }
    err=fread(&im, sizeof(double), 1, fp);
    if(err!=1)
      {
      fprintf(stderr, "Problems reading double complex from file (%s, %d)\n", __FILE__, __LINE__);
      return 1;
      }

    aux[0]=re;
    aux[1]=im;

    memcpy((void *)r, (void*)aux, sizeof(aux));
    //equivalent to r=re+im*I;
    }

  return 0;
  }



#endif
