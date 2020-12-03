#include <stdio.h>
#include <stdlib.h>
/**#include <math.h>**/
#include <Rcpp.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif



using namespace Rcpp;


double* R2C_mat(NumericMatrix a) {
	int n = a.nrow();
	int m = a.ncol();
	double* a_c = (double*)malloc(sizeof(double)*n*m);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
        		a_c[i*m + j] = a(i,j);
		}
  	}
	return a_c;
}


double* R2C_vec(NumericVector a) {
	int n = a.size();
	double* a_c = (double*)malloc(sizeof(double)*n);
	for(int i = 0; i < n; ++i) {
        	a_c[i] = a[i];
  	}
	return a_c;
}



void rcpp_information_matrix_c(double* __restrict__ B, double* __restrict__ W, double* __restrict__ expb, double* __restrict__ T, double* __restrict__ R, double* __restrict__  D,double* __restrict__  X, int n, int r,double* __restrict__ ret) {
	double S0 = 0;
	double S1[r];
	double S2[r][r];
	memset(S1,0, sizeof(double)*r);
	memset(S2,0, sizeof(double)*r*r);

	for (int i=0;i<n;i++) {
		double expbw = W[i]*expb[i];
		if (R[i] == 1) {
			S0 += expbw;
			for (int j = 0; j<r;j++) {
				S1[j] += expbw*X[i*r + j];
				for (int k = j; k<r;k++) {
					S2[j][k] += expbw*X[i*r + j]*X[i*r + k];
				}
			}
		} else {
			if (D[i]) {
				for (int j = 0; j<r;j++) {
					for (int k = j; k<r;k++) {
						ret[j*r + k] -= ((S2[j][k]/S0) - (S1[j]/S0) * (S1[k]/S0));
					}
				}		
			}

			S0 -= expbw;
			for (int j = 0; j<r;j++) {
				S1[j] -= expbw*X[i*r + j];
				for (int k = j; k<r;k++) {
					S2[j][k] -= expbw*X[i*r + j]*X[i*r + k];
				}
			}
		}
	}
	
}

// [[Rcpp::export]]
NumericMatrix rcpp_information_matrix(NumericVector B, NumericVector W,  NumericVector expb,NumericVector T,NumericVector R , NumericVector D, NumericMatrix X) {
	int n = X.nrow();
	int r = X.ncol();
	double* B_c = R2C_vec(B);
	double* W_c = R2C_vec(W);
	double* expb_c = R2C_vec(expb);
	double* T_c = R2C_vec(T);
	double* R_c = R2C_vec(R);
	double* D_c = R2C_vec(D);
	double* X_c = R2C_mat(X);

	
	double* ret = (double*)malloc(sizeof(double)*r * r);
	memset(ret,0,sizeof(double)*r * r);
	rcpp_information_matrix_c(B_c,W_c,expb_c,T_c,R_c,D_c,X_c,n,r,ret);
	NumericMatrix out(r,r);

	for(int i = 0; i <r; ++i) {
		for(int j = 0; j  < r; ++j) {
	    		out(i,j) = ret[i*r + j];
		}
	}
	free(B_c);
	free(W_c);
	free(expb_c);
	free(T_c);
	free(D_c);
	free(R_c);
	free(X_c);
	free(ret);
	return out;
}





void rcpp_score_residual_c(double* __restrict__ B, double* __restrict__ W, double* __restrict__ expb, double* __restrict__ T,double* __restrict__ R,double* __restrict__ RL, double* __restrict__  D,double* __restrict__  X, int n, int c, int r,double* __restrict__ ret) {
	double S0 = 0;
	double part_a = 0;
	double part_a_at_trunc[n];
	double* part_b_at_trunc = (double*)malloc(sizeof(double)*r*n);
	double S1[r];
	double part_b[r];
	memset(S1,0, sizeof(double)*r);
	memset(part_b,0, sizeof(double)*r);
	memset(part_a_at_trunc,0, sizeof(double)*n);
	memset(part_b_at_trunc,0, sizeof(double)*r*n);
	int found_c = 0;
	for (int i=0;i<n;i++) {
		if (R[i] == 1) {
			double expbw = W[i]*expb[i];
			S0 += expbw;
			part_a_at_trunc[i] = part_a;
			for (int j = 0; j<r;j++) {
				S1[j] += expbw*X[i*r + j];
				part_b_at_trunc[i*r + j] =  part_b[j];
			}

		}
		else {
			if (D[i]) {
				part_a += (1/S0);
				for (int k = 0; k<r;k++) {
					part_b[k] +=  (S1[k]/pow(S0,2));
				}			
			}
			else {
				int trunc_loc = int(RL[i]);
				for (int k = 0; k<r;k++) {
					ret[found_c*r + k] += expb[i]*(X[i*r +k]*(part_a-part_a_at_trunc[trunc_loc]) - (part_b[k]-part_b_at_trunc[trunc_loc*r + k]));
				}
				found_c += 1;		
			}
			double expbw = W[i]*expb[i];
			S0 -= expbw;
			for (int j = 0; j<r;j++) {
				S1[j] -= expbw*X[i*r + j];
			}
		}
	}
	free(part_b_at_trunc);
	
}

// [[Rcpp::export]]
NumericMatrix rcpp_score_residual(NumericVector B, NumericVector W,  NumericVector expb,NumericVector T,NumericVector R,NumericVector RL , NumericVector D, NumericMatrix X, int c) {
	
	int n = X.nrow();
	int r = X.ncol();

	double* B_c = R2C_vec(B);
	double* W_c = R2C_vec(W);
	double* expb_c = R2C_vec(expb);
	double* T_c = R2C_vec(T);
	double* R_c = R2C_vec(R);
	double* RL_c = R2C_vec(RL);
	double* D_c = R2C_vec(D);
	double* X_c = R2C_mat(X);
	
	
	double* ret = (double*)malloc(sizeof(double)*c * r);
	memset(ret,0,sizeof(double)*c * r);
	
	rcpp_score_residual_c(B_c,W_c,expb_c,T_c,R_c,RL_c,D_c,X_c,n,c,r,ret);
	NumericMatrix out(c,r);

	for(int i = 0; i <c; i++) {
		for(int j = 0; j  < r; j++) {
	    		out(i,j) = ret[i*r + j];
		}
	}
	free(B_c);
	free(W_c);
	free(expb_c);
	free(T_c);
	free(R_c);
	free(RL_c);
	free(D_c);
	free(X_c);
	free(ret);
	return out;
}

