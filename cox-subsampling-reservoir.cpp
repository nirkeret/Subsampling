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
	/**double sum = 0;
	for (int i = 0; i < n; i++) {
		if (i % 1000000 == 0) {
			printf("%f\n",sum);
		}
		sum += ret[i];
	}
	for (int i = 0; i < n; i++) {
		ret[i] /= sum;
	}	**/
	return(ret);
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_L_OPT(NumericMatrix res) {
	int n = res.nrow();
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
Rcpp::NumericVector matRows2Norms(NumericMatrix res) {
	double n = res.nrow();
	double nCov = res.ncol();
	NumericVector resVec(n);	
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < nCov; j++) {
			sum += res(i,j) * res(i,j);
		}
		resVec(i) = sqrt(sum);
	}
	return resVec;
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
		if (i >= n) { //checking if truncation or event;
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

#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace std;

const int mcols = 5;
const char* columns[mcols] = {"ID","truncTime","sortTime","isTrunc","D"};
const char* alt_columns[mcols] = {"\"ID\"","\"truncTime\"","\"sortTime\"","\"isTrunc\"","\"D\""};

// [[Rcpp::export]]
Rcpp::List get_file_param(string filename, int use_truncs=1) {
    //File pointer
    std::ifstream fin;
    //open an existing file
    fin.open(filename, ios::in);
    string line;
    int m = 0;
    int n = 0;
    int d = 0;
    int rs = 0;
    if (!getline(fin, line)) {
        throw "empty or non existant file\n";
    }


    std::istringstream sin(line);
    string col;

    for (int i = 0; i < mcols; i++) {
        if (!getline(sin, col, ',')) {
            cout << "missing column " << columns[i] << "\n";
            throw "missing column\n";
        }
        if (col != columns[i] && col != alt_columns[i]) {
            cout << "missing column " << columns[i] << "\n";
            throw "missing column\n";
        }
    }

    while (getline(sin, col, ',')) m+=1;


    string tmp;
    string D;
    while (getline(fin, line)) {
        n += 1;
        // Parse line
        std::istringstream sin(line);
        getline(sin, tmp, ',');
        getline(sin, tmp, ',');
        getline(sin, tmp, ',');
        getline(sin, tmp, ',');
        getline(sin, D, ',');
        if (D == "1") d += 1;
    }

    if (use_truncs) n /= 2;

    return Rcpp::List::create(Rcpp::Named("n") = n,Rcpp::Named("d") = d,Rcpp::Named("m") = m);
}


// [[Rcpp::export]]
Rcpp::List score_from_file(string filename, NumericVector beta, int n, int m, int d, int use_truncs, string op, NumericMatrix inv_hess=NumericMatrix(1,1)) {
    double* B = R2C_vec(beta);
    //File pointer
    std::ifstream fin;
    //open an existing file
    fin.open(filename, ios::in);
    string line;
    getline(fin, line); //Skip headers
    int id = 0;
    double truncTime = 0;
    double sortTime = 0;
    string isTrunc;
    string D;
    double X[m];
    string tmp;

    int found_c = 0;


    int residuals = 0;
    int hess = 0;
    int L = 1;
    double inv_hess_c[m][m];
    if (op == "hess") {
        hess =1;
    } 
    if (op == "L") {
        L = 1;
        residuals = 1;
    }
    if (op == "A") {
        L = 0;
        residuals = 1;
		for (int j = 0; j<m;j++) {
		    for (int k = 0; k<m;k++) {
			    inv_hess_c[j][k] = inv_hess(j,k);
		    }
		}
    }




	double* W = (double*)malloc(sizeof(double)*n);
	for (int i = 0; i < n; i++) {
		W[i] = 1;
	}



    int cur_event_number = 0;


	double* part_a_at_event = NULL;
	double* part_b_at_event = NULL;
	int* trunc_to_prev_event = NULL;

	if (residuals && use_truncs) {
		part_b_at_event = (double*)malloc(sizeof(double)*m*(d+1));
		memset(part_b_at_event,0, sizeof(double)*m*(d+1));
		part_a_at_event = (double*)malloc(sizeof(double)*(d+1));
		memset(part_a_at_event,0, sizeof(double)*(d+1));
	    trunc_to_prev_event = (int*)malloc(sizeof(int)*n);
        memset(trunc_to_prev_event,0, sizeof(int)*n);
	}

	double S0 = 0;
	double part_a = 0;
	double S1[m];
	double S2[m][m];
	double part_b[m];
	memset(S1,0, sizeof(double)*m);
	memset(S2,0, sizeof(double)*m*m);
	memset(part_b,0, sizeof(double)*m);


    // When there are ties we only update these after all ties have been prosseced.
    // We keep the data in this buffers and than add them to the global counters when ties are done?
	double buff_S0 = 0;
	double buff_S1[m];
	double buff_S2[m][m];
	memset(buff_S1,0, sizeof(double)*m);
	memset(buff_S2,0, sizeof(double)*m*m);
    double prev_t = -1;
    int flush = 0; 


    //Return values
    int c = n-d;
	double* ret_probs = NULL;
	int* ret_pos = NULL;
	double* ret_hess = NULL;
	if (residuals) {
		ret_probs = (double*)malloc(sizeof(double)*c);
		memset(ret_probs,0,sizeof(double)*c);
		ret_pos = (int*)malloc(sizeof(int)*c);
		memset(ret_pos,0,sizeof(int)*c);
	}
	if (hess) {
		ret_hess = (double*)malloc(sizeof(double)*m * m);
		memset(ret_hess,0,sizeof(double)*m * m);
	}

    double tmp_res_for_col = 0;
    int line_n = 0;
    while (getline(fin, line)) {
        line_n += 1;
        // Parse line
        std::istringstream sin(line);
        getline(sin, tmp, ',');
        sscanf(tmp.c_str(),"%d",&id);
        id -=1; //Zero index

        getline(sin, tmp, ',');
        sscanf(tmp.c_str(),"%lf",&truncTime);

        getline(sin, tmp, ',');
        sscanf(tmp.c_str(),"%lf",&sortTime);
        

        getline(sin, isTrunc, ',');
        getline(sin, D, ',');

        double expb = 0;
        for (int l = 0; l < m; l++) {
            getline(sin, tmp, ',');
            sscanf(tmp.c_str(),"%lf",X + l);
            expb += X[l]*B[l];
        }
        
        expb = exp(expb);

        //If we are not tied with previous event/cen or it is a trunc we flush the buff
        if (flush && (prev_t != sortTime || isTrunc == "1")) {
			S0 += buff_S0;
			for (int j = 0; j<m;j++) {
				S1[j] += buff_S1[j];
				if(hess) {
					for (int k = j; k<m;k++) {
						S2[j][k] += buff_S2[j][k];
					}
				}
			}
            flush = 0;
            buff_S0 = 0;
        	memset(buff_S1,0, sizeof(double)*m);
            if(hess) memset(buff_S2,0, sizeof(double)*m*m);
        }

		if (isTrunc == "1") { //checking if truncation or event;
			double expbw = W[id]*expb;
			S0 += expbw;
			if (residuals && use_truncs) trunc_to_prev_event[id] = cur_event_number;
			for (int j = 0; j<m;j++) {
				S1[j] += expbw*X[j];
				if (hess) {
					for (int k = j; k<m;k++) {
						S2[j][k] += expbw*X[j]*X[k];
					}
				}
			}
		}
		else {
            prev_t  = sortTime;
			if (D == "1") {
                cur_event_number += 1;
				if (hess) {
					for (int j = 0; j<m;j++) {
						for (int k = j; k<m;k++) {
							ret_hess[j*m + k] -= ((S2[j][k]/S0) - (S1[j]/S0) * (S1[k]/S0));
						}
					}
				}
				if (residuals) {		
					part_a += (1/S0);
        			if (use_truncs) part_a_at_event[cur_event_number] = part_a;
					for (int k = 0; k<m;k++) {
						part_b[k] +=  (S1[k]/pow(S0,2));
                        if (use_truncs) part_b_at_event[cur_event_number*m + k] =  part_b[k];
					}
				}
			}


			else {
				if (residuals) {
                    double resid[m];
					for (int k = 0; k<m;k++) {
						if (use_truncs) {
                            int ev_index = trunc_to_prev_event[id];                           
                            resid[k] = expb*(X[k]*(part_a-part_a_at_event[ev_index]) - (part_b[k]-part_b_at_event[ev_index*m + k]));
						} else {
							resid[k] = expb*(X[k]*(part_a) - (part_b[k]));
						}
                        if (L) {
                            ret_probs[found_c] += resid[k]*resid[k];
                        }   
					}
                    if (!L) {
                        double residA[m];
	                    for (int j = 0; j < m; j++) {
                            residA[j] = 0;
		                    for (int k = 0; k < m; k++) {
			                    residA[j] += resid[k]*inv_hess_c[j][k];
		                    }
                            ret_probs[found_c] += residA[j]*residA[j];
	                    }
                    }               
                    //ret_probs[found_c] += tmp_res_for_col*tmp_res_for_col;
                    //printf("%d %d\n",found_c,c);
                    ret_probs[found_c] = sqrt(ret_probs[found_c]);
                    //if (found_c == 0) printf("%f %d %d %d %d\n",ret_probs[found_c],found_c,trunc_to_prev_event[id],residuals && use_truncs,id);
                    ret_pos[found_c] = line_n;
                    found_c += 1;
				}
			}
			double expbw = W[id]*expb;
			buff_S0 -= expbw;
			for (int j = 0; j<m;j++) {
				buff_S1[j] -= expbw*X[j];
				if(hess) {
					for (int k = j; k<m;k++) {
						buff_S2[j][k] -= expbw*X[j]*X[k];
					}
				}
			}
            flush = 1;
		}
	}

    free(W);
	NumericVector out_prob(c);
	NumericVector out_prob_ord(c);
	NumericMatrix out_hess(m,m);
    if (residuals) {
	    double sum = 0;
	    for (int i = 0; i < c; i++) {
		    sum += ret_probs[i];
	    }
	    for (int i = 0; i < c; i++) {
		    ret_probs[i] /= sum;
	    }	
        

	    for (int i = 0; i < c; i++) {
		    out_prob[i] = ret_probs[i];
            out_prob_ord[i] = ret_pos[i];
	    }
        free(ret_probs);
        free(ret_pos);
	    if (use_truncs) {
            free(part_a_at_event);
            free(part_b_at_event);
            free(trunc_to_prev_event);
        }

        return Rcpp::List::create(Rcpp::Named("samp_prob") = out_prob,Rcpp::Named("ord") = out_prob_ord);
    }
    else {
	    for(int i = 0; i <m; ++i) {
		    for(int j = 0; j  < m; ++j) {
	        		out_hess(i,j) = ret_hess[i*m + j];
		    }
	    }

        free(ret_hess);
	    return Rcpp::List::create(Rcpp::Named("hess") = out_hess);
    }
}


// [[Rcpp::export]]
void select_from_file(string in_filename, string out_file, NumericVector ind) {
    //File pointer
    std::ifstream fin;
    std::ofstream fout;
    //open an existing file
    fin.open(in_filename, ios::in);
    fout.open(out_file, ios::out);
    string line;
    int n = 0;
    int m = 0;
    if (!getline(fin, line)) {
        throw "empty or non existant file\n";
    }
    fout << line << "\n";


    while (getline(fin, line)) {
        n += 1;
        while (n == ind[m]) {
            fout << line << "\n";
            m+=1;
        }
    }
}



















