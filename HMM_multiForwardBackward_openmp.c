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

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/* Input Arguments */
#define	Q_IN    prhs[0]
#define	H_IN    prhs[1]
#define	IEND_IN prhs[2]

/* output arguments */
#define	LNZ_OUT     plhs[0]
#define	WA_OUT      plhs[1]
#define	PST_OUT     plhs[2]

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    /* input parameters */
    double *H,*Q;
    double *iEnd;
    /* output parameters */
    double *lnZ,*wA,*pst;
    /* temporary variables */
    int Nends,Mends;
    int N,T,tStart,tEnd,n,t,j,k;
    double Za,Zb,ZW;
    mxArray *mx_alpha,*mx_beta,*mx_P;
    double  *alpha,*beta,*P,*watemp;
    unsigned int doBackward,doOccProb;
    
    /* check number of input/output arguments */
    if (nlhs > 3 )
        mexErrMsgTxt("At most three output arguments can be handled.");
    if (nlhs == 1){
        doBackward=0;doOccProb=0;}
    if (nlhs == 2){
        doBackward=1;doOccProb=0;}
    if (nlhs == 3){
        doBackward=1;doOccProb=1;}        
    
    if (nrhs != 3)
        mexErrMsgTxt("Three input arguments required.");
    
    /* size of input variables */
    T = mxGetM(H_IN); /* number of rows    */
    N = mxGetN(H_IN); /* number of columns */
    /* printf("input H is [ %d %d ] \n",T,N); */
    
    /* require iEnd to be a 1 x N matrix */
    Mends = mxGetM(IEND_IN); /* number of rows */
    Nends = mxGetN(IEND_IN); /* number of rows */
    /* printf("input iEnds is [ %d %d ] \n",Mends,Nends); */
    if( (Mends==1) && (Nends >=1) ){ /* then do nothing */
    }else if( (Nends == 1) && (Mends >=1) ){ /* a N x 1 is also OK */
        Nends=Mends;
        Mends=1;
    }else{ /*something weird going on */
        mexErrMsgTxt("iEnds must be a 1-by-N or N-by-1 matrix.");
    }
    
    /* retrieve input data */
    H=mxGetPr(H_IN);
    Q=mxGetPr(Q_IN);
    iEnd=mxGetPr(IEND_IN);
    /* check H and Q have consistent sizes */
    if ( mxGetM(Q_IN) != N || mxGetN(Q_IN) != N)
        mexErrMsgTxt("Q is not an N by N matrix.");
    
    if( iEnd[Nends-1] > T)
        mexErrMsgTxt("iEnds outside range of data.");
    
    /* Create an mxArray for the output data */
    LNZ_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    lnZ=mxGetPr(LNZ_OUT); /* create pointer to real data in output */
    lnZ[0]=0.0;
    if(doBackward>0){
        WA_OUT  = mxCreateDoubleMatrix(N, N, mxREAL);
        wA=mxGetPr(WA_OUT);
    }
    if (doOccProb>0){
        PST_OUT=mxCreateDoubleMatrix(T, N, mxREAL);
        pst=mxGetPr(PST_OUT);
    }
    
    /* create temporary matrices */
    alpha = (double*)calloc(T*N, sizeof(double));
    beta = (double*)calloc(T*N, sizeof(double));
    

    /*mx_alpha=mxCreateDoubleMatrix(T, N, mxREAL);
    alpha=mxGetPr(mx_alpha);
    mx_beta=mxCreateDoubleMatrix(T, N, mxREAL);
    beta=mxGetPr(mx_beta);
    mx_P=mxCreateDoubleMatrix(N, N, mxREAL);
    P=mxGetPr(mx_P);*/
    
    const int threads = 4;
    omp_set_num_threads(threads);
    int id;
    double lnztemp[4] = {0.0, 0.0, 0.0, 0.0};
    P = (double*)calloc(N*N*threads, sizeof(double));
    watemp = (double*)calloc(N*N*threads, sizeof(double));
    
    /* actual forward-backward sweeps */
    #pragma omp parallel for schedule(dynamic) private(Za,Zb,j,k,t,ZW,tStart,tEnd,n,id) shared(H,Q,iEnd)
    for(n=0;n<Nends;n++){
    	id = omp_get_thread_num();
     	/*mexPrintf("n=%d, id=%d. ", n, id);*/

    	if (n==0)
    		tStart=0;
    	else
    		tStart=(int)(iEnd[n-1]);

        tEnd=(int)(iEnd[n]);
        
        /* forward sweep */
        /*Za=sum(H(tStart,:));*/
        /*alpha(tStart,:)=H(tStart,:)/Za;*/
        Za=0.0;
        for(j=0;j<N;j++){
            Za=Za+H[tStart+T*j];}
        for(j=0;j<N;j++){
            alpha[tStart+T*j]=H[tStart+T*j]/Za;}
        lnztemp[id]=lnztemp[id]+log(Za);
        for(t=tStart+1;t<tEnd;t++){
            /* alpha(t,:)=(alpha(t-1,:)*Q).*qst(t,:);*/
            /* Za=sum(alpha(t,:));*/
            Za=0;
            for(j=0;j<N;j++){
                for(k=0;k<N;k++){
                    alpha[t+j*T]=alpha[t+j*T]+alpha[t-1+k*T]*Q[k+j*N]*H[t+j*T];
                }
                Za=Za+alpha[t+j*T];
            }
            /* alpha(t,:)=alpha(t,:)/Za(t); */
            for(j=0;j<N;j++){
                alpha[t+j*T]=alpha[t+j*T]/Za;}
            /* lnZ=lnZ+log(Za);*/
            lnztemp[id]=lnztemp[id]+log(Za);
        }

        /*  backward sweep */
        if(doBackward>0){
            /*beta(tEnd,:)=ones(1,N)/N;*/
            for(j=0;j<N;j++){
                beta[tEnd-1+j*T]=1.0/N;
            }
            for (t=tEnd-2;t>=tStart;t--){
	      /* beta(t,:)=(beta(t+1,:).*H(t+1,:))*QT; */
	      /* Zb=sum(beta(t,:)); */
	      /* beta(t,:)=beta(t,:)/Zb;*/
                Zb=0;
                for(j=0;j<N;j++){
                    for(k=0;k<N;k++){
                        beta[t+j*T]=beta[t+j*T]+beta[t+1+k*T]*H[t+1+k*T]*Q[j+k*N];
                    }
                    Zb=Zb+beta[t+j*T];
                }
                for(j=0;j<N;j++){
                    beta[t+j*T]=beta[t+j*T]/Zb;
                }
            }

            /* transition counts */
	        for(t=tStart+1;t<tEnd;t++){
	            ZW=0.0;
	            for(j=0;j<N;j++){
	                for(k=0;k<N;k++){
	                    P[id*N*N+j+k*N]=alpha[(t-1)+j*T]*Q[j+k*N]*H[t+k*T]*beta[t+k*T];
	                    ZW=ZW+P[id*N*N+j+k*N];
	                }
	            }
	            for(j=0;j<N;j++){
	                for(k=0;k<N;k++){
	                    watemp[id*N*N+j+k*N]=watemp[id*N*N+j+k*N]+P[id*N*N+j+k*N]/ZW;
	                }
	            }
	        }
        }
    }

    lnZ[0] = lnztemp[0] + lnztemp[1] + lnztemp[2] + lnztemp[3];

    for (id=0;id<threads;id++){
        for(j=0;j<N;j++){
            for(k=0;k<N;k++){
				wA[j+k*N] += watemp[id*N*N+j+k*N];
			}
		}
	}

    /* do occupation probability if asked for */
    if(doOccProb>0){
        /* pst = rowNormalize(alpha.*beta);*/
        for(t=0;t<T*N;t++){
            pst[t]=alpha[t]*beta[t];
        }
        for(t=0;t<T;t++){
            ZW=0.0;
            for(j=0;j<N;j++){
                ZW=ZW+pst[t+j*T];
            }
            for(j=0;j<N;j++){
                pst[t+j*T]=pst[t+j*T]/ZW;
            }
        }
    }
    
    
    /* destroy temporary vriables */
    free(alpha);
    free(beta);
    free(P);
    free(watemp);
    /*mxDestroyArray(mx_alpha);
    mxDestroyArray(mx_beta);
    mxDestroyArray(mx_P);*/
}
