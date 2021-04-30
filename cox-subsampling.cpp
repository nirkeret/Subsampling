#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
/**#include <math.h>**/
#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#if defined(_OPENMP)
#include <omp.h>
#endif



using namespace Rcpp;
using namespace std;

double* R2C_mat(NumericMatrix a) {
	int n = a.nrow();
	int m = a.ncol();
	double* a_c = (double*)malloc(sizeof(double)*n*m);

	for(long i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
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


struct less_than_key
{
    double *T;
    double *D;
    int n;
    less_than_key(double* T, double* D, int n) {
	this->T = T;
	this->D = D;
	this->n = n;
    }
    inline bool operator() (const int i, const int j)
    {
        if (T[i]<T[j]) return(1);
	if (T[i] == T[j]) {
		if ((i < n && D[i]) && (j >=n || D[j] == 0)) {
			return 1;
		}
		if (i <n && j >= n) {
			return 1;
		}		
	}
	return(0);
    }
};


double* L2_norm_norm(double* score, int n, int m) {
	double* ret = (double*)malloc(sizeof(double)*n);
	memset(ret,0,sizeof(double)*n);
	for (long i = 0; i < n; i++) {
		for (int k = 0; k < m; k++) {
			ret[i] += score[i*m + k]*score[i*m + k];
		}
		ret[i] = sqrt(ret[i]);
	}
	double sum = 0;
	for (int i = 0; i < n; i++) {
		/**if (i % 1000000 == 0) {
			printf("%f\n",sum);
		}**/
		sum += ret[i];
	}
	for (int i = 0; i < n; i++) {
		ret[i] /= sum;
	}	
	return(ret);
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_L_OPT(NumericMatrix res) {
	int n =res.nrow();
	int m = res.ncol();
	double* res_c = R2C_mat(res);
	double *ret_val = L2_norm_norm(res_c,n,m);
	NumericVector out(n);
	for (int i = 0; i < n; i++) {
		out(i) = ret_val[i];
	}
	free(res_c);
	free(ret_val);
	return(out);
}



// [[Rcpp::export]]
Rcpp::NumericVector rcpp_A_OPT(NumericMatrix res,NumericMatrix inv_v) {
	int n =res.nrow();
	int m = res.ncol();
	double* res_c = R2C_mat(res);
	double* inv_v_c = R2C_mat(inv_v);
	double* res_inv = (double*)malloc(sizeof(double)*n*m);
	memset(res_inv,0,sizeof(double)*n*m);
	for (long i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < m; k++) {
				res_inv[i*m + j] += res_c[i*m + k]*inv_v_c[j*m + k];
			}
		}
	}

	double *ret_val = L2_norm_norm(res_inv,n,m);
	NumericVector out(n);
	for (int i = 0; i < n; i++) {
		out(i) = ret_val[i];
	}
	free(res_c);
	free(ret_val);
	free(inv_v_c);
	free(res_inv);
	return(out);
}


// [[Rcpp::export]]
Rcpp::List rcpp_score_wrapper(NumericVector beta, NumericVector weights,NumericVector times,NumericVector truncs, NumericVector status, NumericMatrix covariates,int hess,int residuals) {

	double* W = NULL;
	double* B = R2C_vec(beta);
	double* times_c = R2C_vec(times);
	double* truncs_c = NULL;
	double* D = R2C_vec(status);
	double* X = R2C_mat(covariates);
	int n = covariates.nrow();
	int r = covariates.ncol();
	int c = 0;
	int use_truncs = 0;
	for (int i = 0; i < n;i++) {
		c += D[i] == 0;
	}

	if (weights[0] > 0) {
		W = R2C_vec(weights);
	}
	else {
		W = (double*)malloc(sizeof(double)*n);
		for (int i = 0; i < n; i++) {
			W[i] = 1;
		}
	}
	if (truncs[0] >= 0) {
		truncs_c = R2C_vec(truncs);
		use_truncs = 1;
	}

	double* expb = (double*)malloc(sizeof(double)*n);
	for (long i = 0; i < n;i++) {
		double sm = 0;
		for (int j = 0; j < r;j++) {
			sm += beta[j] * X[i*r + j];
		}
		expb[i] = exp(sm);
	}
	double* ret_residual = NULL;
	double* ret_hess = NULL;
	if (residuals) {
		ret_residual = (double*)malloc(sizeof(double)*c * r);
		memset(ret_residual,0,sizeof(double)*c * r);
	}
	if (hess) {
		ret_hess = (double*)malloc(sizeof(double)*r * r);
		memset(ret_hess,0,sizeof(double)*r * r);
	}
	double* T = (double*)malloc(sizeof(double)*2 * n);
	int* ord = (int*)malloc(sizeof(int)*2 * n);
	memcpy(T,times_c,sizeof(double)*(n));
	if (use_truncs) {
		memcpy(T+n,truncs_c,sizeof(double)*(n));
	} else {
		memset(T+n,0,sizeof(double)*(n));
	}
	for (int i = 0; i< 2*n; i++) {
		ord[i] = i;
	}
	less_than_key cmp(T,D,n);
	//sort(ord,ord + (2*n), [&](int i,int j){return T[i]<T[j];} );
	sort(ord,ord + (2*n), cmp);
	double S0 = 0;
	double part_a = 0;
	double* part_a_at_trunc = NULL;
	double* part_b_at_trunc = NULL;
	if (residuals && use_truncs) {
		part_b_at_trunc = (double*)malloc(sizeof(double)*r*n);
		memset(part_b_at_trunc,0, sizeof(double)*r*n);
		part_a_at_trunc = (double*)malloc(sizeof(double)*n);
		memset(part_a_at_trunc,0, sizeof(double)*n);
	}
	double S1[r];
	double S2[r][r];
	double part_b[r];
	memset(S1,0, sizeof(double)*r);
	memset(S2,0, sizeof(double)*r*r);
	memset(part_b,0, sizeof(double)*r);

	long found_c = 0;
	long i = 0;
	for (long o=0;o<n*2;o++) {
		i = ord[o];
		if (i >= n) {
			long ip = i - n;
			double expbw = W[ip]*expb[ip];
			S0 += expbw;
			if (residuals && use_truncs) part_a_at_trunc[ip] = part_a;
			for (int j = 0; j<r;j++) {
				S1[j] += expbw*X[ip*r + j];
				if (residuals && use_truncs) part_b_at_trunc[ip*r + j] =  part_b[j];
				if (hess) {
					for (int k = j; k<r;k++) {
						S2[j][k] += expbw*X[ip*r + j]*X[ip*r + k];
					}
				}
			}
		}
		else {
			double cur_T  = T[i];
			long end_o = o;
			while (T[i] == cur_T && i < n && end_o < n*2) {
				if (D[i]) {
					if (hess) {
						for (int j = 0; j<r;j++) {
							for (int k = j; k<r;k++) {
								ret_hess[j*r + k] -= ((S2[j][k]/S0) - (S1[j]/S0) * (S1[k]/S0));
							}
						}
					}
					if (residuals) {		
						part_a += (1/S0);
						for (int k = 0; k<r;k++) {
							part_b[k] +=  (S1[k]/pow(S0,2));
						}
					}			
				}
				else {
					if (residuals) {
						for (int k = 0; k<r;k++) {
							if (use_truncs) {
								ret_residual[found_c*r + k] += expb[i]*(X[i*r +k]*(part_a-part_a_at_trunc[i]) - (part_b[k]-part_b_at_trunc[i*r + k]));
							} else {
								ret_residual[found_c*r + k] += expb[i]*(X[i*r +k]*(part_a) - (part_b[k]));
							}
						}
						found_c += 1;		
					}
				}
				end_o += 1;
				if (end_o < n*2) i = ord[end_o];
			}
			for (;o < end_o; o++) {
				i = ord[o];
				double expbw = W[i]*expb[i];
				S0 -= expbw;
				for (int j = 0; j<r;j++) {
					S1[j] -= expbw*X[i*r + j];
					if(hess) {
						for (int k = j; k<r;k++) {
							S2[j][k] -= expbw*X[i*r + j]*X[i*r + k];
						}
					}
				}
			}
			o--;
		}
	}
	if (residuals && use_truncs) {	
		free(part_b_at_trunc);
		free(part_a_at_trunc);
	}
	

	NumericMatrix out_residuals(c,r);
	NumericMatrix out_hess(r,r);
	NumericVector out_ord(n);

	if (residuals)
	for(long i = 0; i <c; i++) {
		for(int j = 0; j  < r; j++) {
	    		out_residuals(i,j) = -ret_residual[i*r + j];
		}
	}
	if (hess) 
	for(int i = 0; i <r; ++i) {
		for(int j = 0; j  < r; ++j) {
	    		out_hess(i,j) = ret_hess[i*r + j];
		}
	}

	int count = 0;
	for(int i = 0; i <2*n; ++i) {
		if (ord[i] < n) {
			out_ord(count) = ord[i] + 1;
			count+=1;
		}
	}
	
	free(W);
	free(times_c);
	free(truncs_c);
	free(D);
	free(X);
	free(ret_hess);
	free(ret_residual);
	free(T);
	free(ord);
	free(B);
	free(expb);
	return Rcpp::List::create(Rcpp::Named("residual") = out_residuals,
                          Rcpp::Named("hess") = out_hess,Rcpp::Named("ord") = out_ord);
}

// [[Rcpp::export]]
NumericVector ids_to_index(NumericVector ids_table, NumericVector ids_search) {

	int n = ids_table.length();
	int m = ids_search.length();		
	int* ret = (int*)malloc(sizeof(int) * n);
	int ret_size = 0;

	int max_id = -1;
	for (int i = 0; i < m; i++) {
		if (ids_search[i] > max_id) max_id = ids_search[i];
	}
	max_id = max_id + 1;
	int* map = (int*)malloc(sizeof(int) * max_id);

	memset(map,0, sizeof(int)*max_id);


	for (int i = 0; i < m; i++) {
		map[int(ids_search[i])] += 1;
	}

	for (int i = 0; i < n; i++) {
		int id = ids_table[i];
		if (id < max_id) {
			for (int j = 0; j < map[id]; j++) {
				ret[ret_size] = i;
				ret_size += 1;
			}
		}
	}

	NumericVector out(ret_size);
	for (int i = 0; i < ret_size; i++) {
		out[i] = ret[i] + 1;
	}
	

	free(ret);
	free(map);
	return out;
}

// [[Rcpp::export]]
NumericVector get_weights(NumericVector ids,NumericVector W, NumericVector ids_all,double q) {

	int n = ids_all.length();	
	double* a = R2C_vec(ids);
	double* b = R2C_vec(ids_all);
	double* ret = (double*)malloc(sizeof(double) * n);
	
	int max_id = -1;
	for (int i = 0; i <  ids.length(); i++) {
		if (a[i] > max_id) max_id = a[i];
	}
	max_id = max_id + 1;
	double* map = (double*)malloc(sizeof(double) * max_id);
	memset(map,0, sizeof(double)*max_id);

	for (int i = 0; i < ids.length(); i++) {
		map[int(a[i])] = W[i];
	}

	for (int i = 0; i < ids_all.length(); i++) {
		if (int(b[i]) < max_id && map[int(b[i])] > 0) {
			ret[i] = 1/(map[int(b[i])]*q);
		} else {
			ret[i] = 1;
		}
	}

	NumericVector out(n);
	for (int i = 0; i < n; i++) {
		out[i] = ret[i];
	}
	
	free(a);
	free(b);
	free(ret);
	free(map);
	return out;
}






























