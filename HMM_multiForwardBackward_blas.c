/*
 * HMM_multiForwardBackward.c
 * =========================================================================
 *
 * Copyright (C) 2014 Martin Lindén, E-mail: bmelinden@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or any later
 * version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * Additional permission under GNU GPL version 3 section 7
 *
 * If you modify this Program, or any covered work, by linking or combining it
 * with Matlab or any Matlab toolbox, the licensors of this Program grant you
 * additional permission to convey the resulting work.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/* [lnZ,wA,pst]=VB_multiForwardBackward(Q,H,iEnd)
 *
 * M.L. 2013-05-24
 *
 */

/* This is a modified version where BLAS routines have been implemented
 * and optimized for use with vbSPT. Only for tests.
 *
 * Marcus Näslund, 2015-06-01 */

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/* Input Arguments */
#define	Q_IN    prhs[0]
#define	H_IN    prhs[1]
#define	IEND_IN prhs[2]
#define SIZE_IN prhs[3]

/* output arguments */
#define	LNZ_OUT     plhs[0]
#define	WA_OUT      plhs[1]
#define	PST_OUT     plhs[2]

/* Prototype for the DGEMM routine */
extern void dgemm(char*, char*, int*, const int*, const int*, double*, double*, int*, const double*, const int*, double*, double*, int*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* check number of input/output arguments */
    static const bool doBackward=true;
    static const bool doOccProb=true;     
    
    /* size of input variables */
    const int T = mxGetM(H_IN); /* number of rows    */
    const int N = mxGetN(H_IN); /* number of columns */
       
    /* retrieve input data */
    const double * const H=mxGetPr(H_IN);
    const double * const Q=mxGetPr(Q_IN);
    const double * const iEnd=mxGetPr(IEND_IN);
    const double * const sizes=mxGetPr(SIZE_IN);
    const int Nends = mxGetN(IEND_IN); /* number of rows */
    const int sizesmax = mxGetN(SIZE_IN);
    
    /* Create an mxArray for the output data */
    LNZ_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    double * const lnZ=mxGetPr(LNZ_OUT);
    lnZ[0]=0.0;

    WA_OUT  = mxCreateDoubleMatrix(N, N, mxREAL);
    double * const wA=mxGetPr(WA_OUT);

    PST_OUT=mxCreateDoubleMatrix(T, N, mxREAL);
    double * const pst=mxGetPr(PST_OUT);
    
    /* create temporary matrices */
    mxArray * const mx_alpha=mxCreateDoubleMatrix(T, N, mxREAL);
    double * const alpha=mxGetPr(mx_alpha);

    mxArray * const mx_beta=mxCreateDoubleMatrix(T, N, mxREAL);
    double * const beta=mxGetPr(mx_beta);

    mxArray * const mx_P=mxCreateDoubleMatrix(N, N, mxREAL);
    double * const P=mxGetPr(mx_P);
    
    /* temporary variables */
    int tStart,tEnd,n,t,i,j,k,z,index;
    double *Za,*Zb,ZW;

    index = 0;

    for (n=0; n<sizesmax; n++) {
    	/* For length n+1... */

    	/* Number of trajectories */
    	int m = (sizes[n+1] - sizes[n])/n;

    	/* Init alpha-zero */
    	/*Za=sum(H(tStart,:));*/
        Za=(double*)calloc(m, sizeof(double));
        for(i=0;i<n;i++) {
            for(j=0;j<N;j++) {
                Za[i] += H[index + (n+1)*j + T*j]; 
            }
        }        

        /*alpha(tStart,:)=H(tStart,:)/Za;*/
        for(i=0;i<n;i++) {
            for(j=0;j<N;j++) {
                alpha[index + (n+1)*j + T*j]=H[index + (n+1)*j + T*j]/Za[i];
            }
            lnZ[0]=lnZ[0]+log(Za[i]);
        }

    	for (t=0; t<n+1; t++) {
    		/* Time step t=0,1,...,n */

    		/* Chug everything from index to index+m*n */

    		/* Forward sweep alpha */

            /*mwSignedIndex m,n,p;*/

            /* call BLAS function */
            /* alpha(t,:)=(alpha(t-1,:)*Q);*/
            /* Za=sum(alpha(t,:));*/
            char *chn = "N";
            double one = 1.0, zero = 0.0;
            dgemm(chn, chn, &m, &N, &N, &one, &alpha[index], &m, Q, &N, &zero, &alpha[index+m*n], &m);

            /* alpha(t,:) = alpha(t,:).*qst(t,:) */
            for(j=0;j<N;j++) {
                for(k=0;k<N;k++) {
                    alpha[index+t+j*T]=alpha[index+t+j*T]*Q[k+j*N]*H[index+t+j*T];
                }
            }

            /* alpha(t,:)=alpha(t,:)/Za(t); */
            for(i=0;i<n;i++){
                Za[i] = 0;

                for (j=0;j<N;j++) {
                    Za[i] += alpha[index + (n+1)*j + j*T];
                }

                for(j=0;j<N;j++) {
                    alpha[index + (n+1)*j + j*T]=alpha[index + (n+1)*j + j*T]/Za[i];
                }
            }
            /* lnZ=lnZ+log(Za);*/
            lnZ[0]=lnZ[0]+log(Za[i]);

    		/* Increase index */
    		index += n*m;
    	}

        free(Za);

    	/* Always true */
    	if (doBackward) {

    		/* Decrease index */
    		index -= m*n;

    		/* Init beta-end */
    		/*beta(tEnd,:)=ones(1,N)/N;*/
            for(j=index;j<index + m*n;j++) {
            	beta[j] = 1.0/n;
            }

            Zb = (double*)calloc(m, sizeof(double));
    		
    		for (t=n-1; t>=0; t--) {
    			/* Time step t=n-1,n-2,...,0 */

    			/* Decrease index */
    			index -= m*n;

    			/* Chug everything from index to index + m*n

    			/* Backward sweep beta */

	            /* beta(t,:)=(beta(t+1,:).*H(t+1,:)); */
                for(j=0;j<N;j++) {
                    for(k=0;k<N;k++) {
                        beta[t+j*T]=beta[t+j*T]+beta[t+1+k*T]*H[t+1+k*T];
                    }
                }

                /* call BLAS function */
                /* beta(t,:)=(beta(t,:)*QT);*/
                char *chn = "N";
                double one = 1.0, zero = 0.0;
                dgemm(chn, chn, &m, &N, &N, &one, &beta[index], &m, Q, &N, &zero, &beta[index], &m);

                /* Zb=sum(beta(t,:)); */
                for(i=0;i<n;i++) {
                    Zb[i] = 0;

                    for(j=0;j<N;j++) {
                        Zb[i] += beta[index + (n+1)*j + j*T];
                    }

                    for(j=0;j<N;j++) {
                        beta[index + (n+1)*j + j*T] = beta[index + (n+1)*j + j*T]/Zb[i]; 
                    }
                }
    		}

            free(Zb);

    		for (i=0; i<m; t++) {
    			/* Trajectory i=0,1,...,m-1 */

    			/* Impossible to BLAS, must do for each trajectory */

    			/* Transition counts */
    			for(t=1;t<n+1;t++) {
                    /* Time step t=1,2,...,n-1 */

	                ZW=0.0;

	                for(j=0;j<N;j++) {
	                    for(k=0;k<N;k++) {
	                        P[j+k*N]=alpha[index + (t-1) + j*T] * Q[j+k*N] * H[index + t + k*T] * beta[index + t + k*T];
	                        ZW=ZW+P[j+k*N];
	                    }
	                }

	                for(j=0;j<N;j++) {
	                    for(k=0;k<N;k++) {
	                        wA[j+k*N]=wA[j+k*N]+P[j+k*N]/ZW;
	                    }
	                }
	            }

    			/* Increase index */
    			index += n*m;
    		}
    	}
    }

    /* do occupation probability if asked for */
    if(doOccProb) {
        /* pst = rowNormalize(alpha.*beta);*/
        for(t=0;t<T*N;t++) {
            pst[t]=alpha[t]*beta[t];
        }

        for(t=0;t<T;t++) {
            ZW=0.0;
            for(j=0;j<N;j++) {
                ZW=ZW+pst[t+j*T];
            }

            for(j=0;j<N;j++) {
                pst[t+j*T]=pst[t+j*T]/ZW;
            }
        }
    }
    
    /* destroy temporary vriables */
    mxDestroyArray(mx_alpha);
    mxDestroyArray(mx_beta);
    mxDestroyArray(mx_P);
}
