/*
 * water molecule with sto-3g basis
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "cint.h"



void run_all(int *atm, int natm, int *bas, int nbas, double *env);

int cint2e_ip1_sph(double *buf, int *shls,
                   int *atm, int natm, int *bas, int nbas, double *env,
                   CINTOpt *opt);
void cint2e_ip1_sph_optimizer(CINTOpt **opt, int *atm, int natm,
                              int *bas, int nbas, double *env);

int main()
{
        int natm = 3;
        int nbas = natm*20;
        /* ATM_SLOTS = 6; BAS_SLOTS = 8; */
        int *atm = malloc(sizeof(int) * natm * ATM_SLOTS);
        int *bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
        double *env = malloc(sizeof(double) * 10000);

        int i, j, ia, n, off;
        off = PTR_ENV_START; 

        /* Geometry */
        atm(CHARGE_OF,0)=8; atm(PTR_COORD,0)=off; env[off+0]= 0.000; env[off+1]=-0.143; env[off+2]= 0.000; off+=3;
        atm(CHARGE_OF,1)=1; atm(PTR_COORD,1)=off; env[off+0]= 1.638; env[off+1]= 1.136; env[off+2]= 0.000; off+=3;
        atm(CHARGE_OF,2)=1; atm(PTR_COORD,2)=off; env[off+0]=-1.638; env[off+1]= 1.136; env[off+2]= 0.000; off+=3;

        /* STO-3G basis for O,H (basissetexchange) 
         * 
         * H     0
         * S    3   1.00
         *       0.3425250914D+01       0.1543289673D+00
         *       0.6239137298D+00       0.5353281423D+00
         *       0.1688554040D+00       0.4446345422D+00
         * ****
         * O     0
         * S    3   1.00
         *       0.1307093214D+03       0.1543289673D+00
         *       0.2380886605D+02       0.5353281423D+00
         *       0.6443608313D+01       0.4446345422D+00
         * SP   3   1.00
         *       0.5033151319D+01      -0.9996722919D-01       0.1559162750D+00
         *       0.1169596125D+01       0.3995128261D+00       0.6076837186D+00
         *       0.3803889600D+00       0.7001154689D+00       0.3919573931D+00
         * ****
        */
        /* OXYGEN S 3 */
        env[off+0 ] = 130.70932; env[off+3 ] = 0.1543289*CINTgto_norm(0,env[off+0 ]);
        env[off+1 ] = 23.808866; env[off+4 ] = 0.5353281*CINTgto_norm(0,env[off+1 ]);
        env[off+2 ] = 6.4436083; env[off+5 ] = 0.4446345*CINTgto_norm(0,env[off+2 ]);
        /* OXYGEN S 3 */
        env[off+6 ] = 5.0331513; env[off+9 ] =-0.0999672*CINTgto_norm(0,env[off+6 ]);
        env[off+7 ] = 1.1695961; env[off+10] = 0.3995128*CINTgto_norm(0,env[off+7 ]);
        env[off+8 ] = 0.3803889; env[off+11] = 0.7001154*CINTgto_norm(0,env[off+8 ]);
        /* OXYGEN P 3 */
        env[off+11] = 5.0331513; env[off+14] = 0.1559162*CINTgto_norm(1,env[off+11]);
        env[off+12] = 1.1695961; env[off+15] = 0.6076837*CINTgto_norm(1,env[off+12]);
        env[off+13] = 0.3803889; env[off+16] = 0.3919573*CINTgto_norm(1,env[off+13]);
        /* H S 3 */
        env[off+17] = 3.4252509; env[off+20] = 0.1543289*CINTgto_norm(0,env[off+17]);
        env[off+18] = 0.6239137; env[off+21] = 0.5353281*CINTgto_norm(0,env[off+18]);
        env[off+19] = 0.1688554; env[off+22] = 0.4446345*CINTgto_norm(0,env[off+19]);
        /* Assign basis to atoms, angular momentum, etc. */
        /* first loop runs over waters (single water), second one over hydrogens */
        for (i = 0, ia = 0, n = 0; i < 1; i++) {
                /* OXYGEN S 3 */
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+0;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+3;
                n++;
                /* OXYGEN S 3 */
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+6;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+9;
                n++;
                /* OXYGEN P 3 */
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+11;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+14;
                n++;
                ia++;
                for (j = 0; j < 2; j++) {
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+17;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+20;
                        n++;
                        ia++;
                }
        }
        nbas = n;
        printf("STO-3G basis\n");
        run_all(atm, natm, bas, nbas, env);
        free(atm);
        free(bas);
        free(env);
}

void run_all(int *atm, int natm, int *bas, int nbas, double *env)
{
        int i, j, k, l, ij, kl;
        int di, dj, dk, dl;
        int kl_max;
        int shls[4];
        double *buf;

        int *ishls = malloc(sizeof(int)*nbas*nbas);
        int *jshls = malloc(sizeof(int)*nbas*nbas);

        /* Build shell pair list */
        for (i = 0, ij = 0; i < nbas; i++) {
                for (j = 0; j <= i; j++, ij++) {
                        ishls[ij] = i;
                        jshls[ij] = j;
                }
        }

        int ncgto = CINTtot_cgto_spheric(bas, nbas);
        printf("\tshells = %d\n", nbas);
        printf("\ttotal cGTO = %d\n", ncgto);
        printf("\ttotal pGTO = %d\n",
               CINTtot_pgto_spheric(bas, nbas));

        /* Build AO offsets for global indexing */
        int *ao_loc = malloc(sizeof(int)*(nbas+1));
        ao_loc[0] = 0;
        for (i = 0; i < nbas; i++) {
                ao_loc[i+1] = ao_loc[i] + CINTcgto_spheric(i, bas);
        }

        printf("\nPrinting ERIs (μ ν | λ σ)\n\n");

        CINTOpt *opt_for_cint2e = NULL;
        cint2e_sph_optimizer(&opt_for_cint2e, atm, natm, bas, nbas, env);

        for (ij = 0; ij < nbas*(nbas+1)/2; ij++) {

                i = ishls[ij];
                j = jshls[ij];

                di = CINTcgto_spheric(i, bas);
                dj = CINTcgto_spheric(j, bas);

                kl_max = (i+1)*(i+2)/2;

                for (kl = 0; kl < kl_max; kl++) {

                        k = ishls[kl];
                        l = jshls[kl];

                        dk = CINTcgto_spheric(k, bas);
                        dl = CINTcgto_spheric(l, bas);

                        shls[0] = i;
                        shls[1] = j;
                        shls[2] = k;
                        shls[3] = l;

                        buf = malloc(sizeof(double) * di*dj*dk*dl);

                        cint2e_sph(buf, shls,
                                   atm, natm,
                                   bas, nbas,
                                   env,
                                   opt_for_cint2e);

                        int ii, jj, kk, ll;
                        int idx = 0;

                        for (ii = 0; ii < di; ii++) {
                                for (jj = 0; jj < dj; jj++) {
                                        for (kk = 0; kk < dk; kk++) {
                                                for (ll = 0; ll < dl; ll++, idx++) {

                                                        if (fabs(buf[idx]) > 1e-12) {

                                                                int shrek = ao_loc[i] + ii;
                                                                int donkey = ao_loc[j] + jj;
                                                                int fiona = ao_loc[k] + kk;
                                                                int ass = ao_loc[l] + ll;

                                                                printf("(%2d %2d | %2d %2d) = % .12e\n",
                                                                       shrek,donkey,fiona,ass,
                                                                       buf[idx]);
                                                        }
                                                }
                                        }
                                }
                        }

                        free(buf);
                }
        }

        CINTdel_optimizer(&opt_for_cint2e);

        free(ao_loc);
        free(ishls);
        free(jshls);
}

