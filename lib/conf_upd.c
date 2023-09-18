#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdlib.h>

#include"../include/conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

// compute the staple for the phi field, such that
// \sum_{mu>0 and <0} lambda_{x,mu} conj(\phi_x).\phi_{x+\mu} = conj(\phi_x).staple
void calcstaples_for_phi(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         Vec *staple)
  {
  int i;
  Vec v1, v2;

  zero_Vec(staple);

  for(i=0; i<STDIM; i++)
     {
     // forward
     #ifdef CSTAR_BC
       if(bcsitep(geo, r, i)==1)
         {
         equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]) );
         }
       else
         {
         equal_cc_Vec(&v1, &(GC->phi[nnp(geo, r, i)]) );
         }
     #else
       equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]) );
     #endif
     times_equal_complex_Vec(&v1, chargepow(GC->lambda[r][i]));
     plus_equal_Vec(staple, &v1);

     // backward
     #ifdef CSTAR_BC
       if(bcsitem(geo, r, i)==1)
         {
         equal_Vec(&v2, &(GC->phi[nnm(geo, r, i)]) );
         times_equal_complex_Vec(&v2, chargepow(conj(GC->lambda[nnm(geo, r, i)][i])));
         }
       else
         {
         equal_cc_Vec(&v2, &(GC->phi[nnm(geo, r, i)]) );
         times_equal_complex_Vec(&v2, chargepow(GC->lambda[nnm(geo, r, i)][i]));
         }
     #else
       equal_Vec(&v2, &(GC->phi[nnm(geo, r, i)]) );
       times_equal_complex_Vec(&v2, chargepow(conj(GC->lambda[nnm(geo, r, i)][i])));
     #endif

     plus_equal_Vec(staple, &v2);
     }
  }



// perform an update of the phi field with metropolis
// retrn 0 if the trial state is rejected and the number of accepted trial
int metropolis_for_phi(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r)
  {
  int acc=0;
  double old_energy, new_energy;
  Vec staple, new_vector;

  calcstaples_for_phi(GC, geo, r, &staple);

  old_energy=-2.0 * (double) NFLAVOUR * param->d_J * creal(scal_prod_Vec(&(GC->phi[r]), &staple));

  rand_rot_single_Vec(&new_vector, &(GC->phi[r]), param->d_epsilon_metro_site);

  new_energy=-2.0 * (double) NFLAVOUR * param->d_J * creal(scal_prod_Vec(&new_vector, &staple));

  #ifdef DEBUG
  double old_energy_aux, new_energy_aux;
  Vec old_vector;
  old_energy_aux= -2.0 * (double) NFLAVOUR * param->d_J * higgs_interaction(GC, geo, param) * (double) STDIM * (double) param->d_volume;
  equal_Vec(&old_vector, &(GC->phi[r]));
  equal_Vec(&(GC->phi[r]), &new_vector);
  new_energy_aux= -2.0 * (double) NFLAVOUR * param->d_J * higgs_interaction(GC, geo, param) * (double) STDIM * (double) param->d_volume;
  equal_Vec(&(GC->phi[r]), &old_vector);
  //printf("%g %g\n", old_energy-new_energy, old_energy-new_energy -(old_energy_aux-new_energy_aux));
  if(fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux))>1.0e-10 )
    {
    fprintf(stderr, "Problem in energy in metropolis for phi (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  if(old_energy>new_energy)
    {
    equal_Vec(&(GC->phi[r]), &new_vector);
    acc+=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         equal_Vec(&(GC->phi[r]), &new_vector);
         acc+=1;
         }

  #ifdef DEBUG
  if(fabs(norm_Vec(&(GC->phi[r]))-1)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in metropolis (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  return acc;
  }


// perform an update of the phi field with overrelaxation
void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            long r)
  {
  double norm;
  double complex aux1;
  Vec staple, newlink;

  calcstaples_for_phi(GC, geo, r, &staple);
  norm=norm_Vec(&staple);

  #ifdef DEBUG
  double prod_before=creal(scal_prod_Vec(&(GC->phi[r]), &staple));
  #endif

  if(norm>MIN_VALUE)
    {
    equal_Vec(&newlink, &staple);
    times_equal_real_Vec(&newlink, 1./norm);
    aux1=scal_prod_Vec(&newlink, &(GC->phi[r]) );
    times_equal_complex_Vec(&newlink, 2.0*aux1);

    minus_equal_Vec(&newlink, &(GC->phi[r]));

    equal_Vec(&(GC->phi[r]), &newlink);
    }
  else
    {
    rand_vec_Vec(&(GC->phi[r]) );
    }

  #ifdef DEBUG
  double prod_after=creal(scal_prod_Vec(&(GC->phi[r]), &staple));

  if(fabs(prod_before-prod_after)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in overrelaxation (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  if(fabs(norm_Vec(&(GC->phi[r]))-1)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in overrelaxation (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif
  }


// staples for the plaquette component of the action
// sum (plaq) = lambda_{x,mu}*staple + independent of lambda_{x,mu}
double complex plaqstaples_for_link(Conf *GC,
                                    Geometry const * const geo,
                                    long r,
                                    int i)
  {
  int j, k;
  double complex ris, tmp;
  long r1;

  ris=0.0;

  for(j=i+1; j<i+STDIM; j++)
     {
     k=j%STDIM;

//              ^ i
//         (6)  |  (3)
//      +---<---+--->---+
//      |       |       |
//   (5)V       ^       V (2)
//      |       |       |
//      +--->---+---<---+---> k
//     r1  (4)  r  (1)

       #ifdef CSTAR_BC
         tmp=conj(GC->lambda[r][k]); // (1)

         if(bcsitep(geo, r, k)==1){tmp*=conj(GC->lambda[nnp(geo, r, k)][i]); } //(2)
         else {tmp*=GC->lambda[nnp(geo, r, k)][i]; }

         if(bcsitep(geo, r, i)==1){tmp*=GC->lambda[nnp(geo,r,i)][k]; } //(3)
         else { tmp*=conj(GC->lambda[nnp(geo,r,i)][k]); }

         ris+=tmp;

         r1=nnm(geo, r, k);

         if(bcsitem(geo,r,k)==1){tmp=GC->lambda[r1][k]; } //(4)
         else{tmp=conj(GC->lambda[r1][k]);}

         if(bcsitem(geo,r,k)==1) {tmp*=conj(GC->lambda[r1][i]); } //(5)
         else {tmp*=GC->lambda[r1][i]; }

         if(bcsitep(geo, r1, i)*bcsitem(geo,r,k)==1){tmp*=conj(GC->lambda[nnp(geo, r1, i)][k]); } //(6)
         else {tmp*=GC->lambda[nnp(geo, r1, i)][k]; }

         ris+=tmp;
       #else
         tmp=conj(GC->lambda[r][k]); //(1)
         tmp*=conj(GC->lambda[nnp(geo, r, k)][i]); //(2)
         tmp*=GC->lambda[nnp(geo,r,i)][k]; //(3)
         ris+=tmp;

         r1=nnm(geo, r, k);
         tmp= GC->lambda[r1][k]; //(4)
         tmp*=conj(GC->lambda[r1][i]); //(5)
         tmp*=conj(GC->lambda[nnp(geo, r1, i)][k]); //(6)
         ris+=tmp;
       #endif
       }

    return ris;
    }


// perform an update with metropolis of the link variables
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i)
  {
  double old_energy, new_energy;
  double complex old_lambda, new_lambda;
  double complex sc, pstaple;
  int acc=0;

  Vec v1;

  #ifdef CSTAR_BC
    if(bcsitep(geo, r, i)==1)
      {
      equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
    else
      {
      equal_cc_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
  #else
    equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
  #endif

  sc=scal_prod_Vec(&(GC->phi[r]), &v1);

  old_lambda=GC->lambda[r][i];

  if(fabs(param->d_K)>MIN_VALUE)
    {
    pstaple=plaqstaples_for_link(GC, geo, r, i);
    }
  else
    {
    pstaple=0.0;
    }

  old_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*chargepow(old_lambda) );
  old_energy-=2.0*param->d_K*creal(old_lambda*pstaple);
  old_energy-= param->d_masssq * creal(old_lambda);
  new_lambda = old_lambda*cexp(I*param->d_epsilon_metro_link*(2.0*casuale()-1));

  new_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*chargepow(new_lambda) );
  new_energy-=2.0*param->d_K*creal(new_lambda*pstaple);
  new_energy-= param->d_masssq * creal(new_lambda);

  #ifdef DEBUG
  double old_energy_aux, new_energy_aux;
  old_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
  old_energy_aux -= 2.0 * (param->d_K)*plaquette(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
  old_energy_aux -= param->d_masssq * creal(old_lambda);

  GC->lambda[r][i] = new_lambda;
  new_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
  new_energy_aux -= 2.0 * (param->d_K)*plaquette(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
  new_energy_aux -= param->d_masssq * creal(new_lambda);
  GC->lambda[r][i] = old_lambda;

  //printf("%g %g\n", old_energy-new_energy, old_energy-new_energy -(old_energy_aux-new_energy_aux));
  if(fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux))>1.0e-10 )
    {
    fprintf(stderr, "Problem in energy in metropolis for link (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  if(old_energy>new_energy)
    {
    GC->lambda[r][i] = new_lambda;
    acc=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         GC->lambda[r][i] = new_lambda;
         acc=1;
         }

  return acc;
  }


// perform a discrete update with metropolis of the link variables
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_link_big(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i)
  {
  #if CHARGE == 1 // if CHARGE==1 the move is trivial
    (void) GC;
    (void) geo;
    (void) param;
    (void) r;
    (void) i;

    int acc=1;
  #else
  double old_energy, new_energy;
  double complex old_lambda, new_lambda;
  double complex sc, pstaple;
  int acc=0, tmp;

  Vec v1;

  #ifdef CSTAR_BC
    if(bcsitep(geo, r, i)==1)
      {
      equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
    else
      {
      equal_cc_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
  #else
    equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
  #endif

  sc=scal_prod_Vec(&(GC->phi[r]), &v1);

  old_lambda=GC->lambda[r][i];

  if(fabs(param->d_K)>MIN_VALUE)
    {
    pstaple=plaqstaples_for_link(GC, geo, r, i);
    }
  else
    {
    pstaple=0.0;
    }

  old_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*chargepow(old_lambda) );
  old_energy-=2.0*param->d_K*creal(old_lambda*pstaple);
  old_energy-= param->d_masssq * creal(old_lambda);

  tmp=1+(int)((CHARGE-1)*casuale());
  new_lambda = old_lambda*cexp(I*PI2/CHARGE*(double)tmp);

  new_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*chargepow(new_lambda) );
  new_energy-=2.0*param->d_K*creal(new_lambda*pstaple);
  new_energy-= param->d_masssq * creal(new_lambda);

  #ifdef DEBUG
  double old_energy_aux, new_energy_aux;
  old_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
  old_energy_aux -= 2.0 * (param->d_K)*plaquette(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
  old_energy_aux -= param->d_masssq * creal(old_lambda);

  GC->lambda[r][i] = new_lambda;
  new_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
  new_energy_aux -= 2.0 * (param->d_K)*plaquette(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
  new_energy_aux -= param->d_masssq * creal(new_lambda);
  GC->lambda[r][i] = old_lambda;

  //printf("%g %g\n", old_energy-new_energy, old_energy-new_energy -(old_energy_aux-new_energy_aux));
  if(fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux))>1.0e-10 )
    {
    fprintf(stderr, "Problem in energy in metropolis for link (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  if(old_energy>new_energy)
    {
    GC->lambda[r][i] = new_lambda;
    acc=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         GC->lambda[r][i] = new_lambda;
         acc=1;
         }
  #endif
  return acc;
  }


// perform a complete update
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc_site,
            double *acc_link,
            double *acc_link_big)

   {
   long r, asum_site, asum_link, asum_link_big;
   int j, dir;
   double complex norm;

   // metropolis on links
   asum_link=0;
   asum_link_big=0;
   #ifndef LINKS_FIXED_TO_ONE
     #ifndef TEMPORAL_GAUGE
     for(r=0; r<param->d_volume; r++)
        {
        for(dir=0; dir<STDIM; dir++)
           {
           asum_link+=metropolis_for_link(GC, geo, param, r, dir);
           asum_link_big+=metropolis_for_link_big(GC, geo, param, r, dir);
           }
        }
     #else
      for(r=0; r<param->d_volume; r++)
        {
        for(dir=1; dir<STDIM; dir++)
           {
           asum_link+=metropolis_for_link(GC, geo, param, r, dir);
           asum_link_big+=metropolis_for_link_big(GC, geo, param, r, dir);
           }
        }
     #endif
   #endif
   *acc_link=((double)asum_link)*param->d_inv_vol;
   *acc_link_big=((double)asum_link_big)*param->d_inv_vol;

   #ifndef TEMPORAL_GAUGE
   *acc_link/=(double)STDIM;
   *acc_link_big/=(double)STDIM;
   #else
   *acc_link/=(double)(STDIM-1);
   *acc_link_big/=(double)(STDIM-1);
   #endif

   // metropolis on phi
   asum_site=0;
   for(r=0; r<param->d_volume; r++)
      {
      asum_site+=metropolis_for_phi(GC, geo, param, r);
      }
   *acc_site=((double)asum_site)*param->d_inv_vol;

   // overrelax on phi
   for(j=0; j<param->d_overrelax; j++)
      {
      for(r=0; r<(param->d_volume); r++)
         {
         overrelaxation_for_phi(GC, geo, r);
         }
      }

   // final unitarization
   for(r=0; r<(param->d_volume); r++)
      {
      unitarize_Vec(&(GC->phi[r]));
      for(dir=0; dir<STDIM; dir++)
         {
         norm=cabs(GC->lambda[r][dir]);
         GC->lambda[r][dir]/=norm;
         }
      }

   GC->update_index++;
   }

#endif
