// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
using namespace arma;
using namespace Rcpp;
using namespace std;

double min_sc(double a, double b) {
  double c;
  if (a <= b) {
    c = a;
  } else {
    c = b;
  }
	return(c);
}

 // [[Rcpp::export]]
mat extr_coef(mat a, int K, int type) {
	
	int arow = a.n_rows;
	
	mat mu = a.submat(0, 0, K - 1, 0);
	
	if (type == 1) {
		return(mu);
	}
	
	mat z = a.submat(K, 0, arow - 1, 0);
	
	mat block_of_z(K*K, 1, fill::zeros);
	int n_blocks = (arow - K) / (K * K);
	mat z_out(K, n_blocks * K, fill::zeros);
	
	
	int r1_s, r2_s, c1_e, c2_e;
	for (int i = 0; i < n_blocks; i++) {
			r1_s = K * K * i;
			r2_s = K * K * (i + 1) - 1;
			
			c1_e = i * K;
			c2_e = (i + 1) * K - 1;
			block_of_z = z.submat(r1_s, 0, r2_s, 0);
			block_of_z.reshape(K, K);
			z_out.submat(0, c1_e, K - 1, c2_e) = trans(block_of_z);

			block_of_z.reshape(K*K, 1);
	}
	
	return(z_out);
}


// [[Rcpp::export]]
mat repmat_cpp(mat Amat, int rows, int cols) {
	
	mat Bmat = repmat(Amat, rows, cols);
	
	return(Bmat);
}

double pow_cust(double a, int pow) {
	mat c(pow, 1, fill::value(a));
	double res = as_scalar(prod(c));
	return(res);
}

// [[Rcpp::export]]
arma::mat inv_c (arma::mat x) {
  arma::mat XX = inv(x);
  return(XX);
}

// [[Rcpp::export]]
arma::mat bols_c (arma::mat x,
                  arma::mat y) {
  arma::mat bols = inv(trans(x) * x) * trans(x) * y;
  return(bols);
}

// [[Rcpp::export]]
arma::mat m_c (arma::mat V,
               arma::mat H,
               arma::mat inv_sig,
               arma::mat B,
               arma::mat X,
               arma::mat bols
) {
  arma::mat m = V * (H * B + kron(inv_sig, trans(X) * X) * bols);
  return(m);
}

// [[Rcpp::export]]
arma::mat v_c (arma::mat H,
               arma::mat inv_sig,
               arma::mat X
) {
  arma::mat v = inv(H + kron(inv_sig, trans(X) * X));
  return(v);
}

// [[Rcpp::export]]
mat lag0(mat x,
         int p
		 ) {
		
		int x_rows = x.n_rows;
		int x_cols = x.n_cols;
		
		mat A(x_rows, x_cols, fill::value(R_NaN));
		
		A.submat(p, 0, x_rows - 1, x_cols - 1) = x.submat(0, 0, x_rows - p - 1, x_cols - 1);
		return(A);
}


// [[Rcpp::export]]
mat remSC(mat x,
          mat rho
		  ) {
		
		int x_rows = x.n_rows;
		int x_cols = x.n_cols;
		
		mat out(x_rows, x_cols, fill::zeros);
		
		for (int j = 0; j < x_cols; j++) {
			for (int i = 1; i < x_rows; i++) {
				out(i, j) = x(i,j) - rho(i,0) * x(i - 1, j);
			}
		}
		
		return(out);
}

// [[Rcpp::export]]
mat do_regs(mat dat, 
			int cut, 
			int p) {
	
	int t = dat.n_rows;
	int M = dat.n_cols;
	
	int K = M + p * M * M;
	
	int ztempc1, ztempc2;
	
	mat Z(t * M, K, fill::zeros);
	mat Ztemp(M, K, fill::zeros);
	Ztemp.submat(0, 0, M - 1, M - 1) = eye(M, M);
	for (int i = cut + p; i < t; i++) {
		for (int j = 0; j < p; j++) {
			ztempc1 = M + j * M * M;
			ztempc2 = M + (j + 1) * M * M - 1;
			Ztemp.submat(0, ztempc1, M - 1, ztempc2) = kron(eye(M, M), dat.row(i - j - 1));
		}
		
		Z.submat((i - 1) * M, 0, i * M - 1, K - 1) = Ztemp;
	}
	return(Z.rows((p - 1)*M, t * M - M - 1));
}

int stability (mat beta,
               int n,
               int l,
			   int ex
			   ) {
	
	mat FF(n*l, n*l, fill::zeros);
	int i = n;
	int j = 0;
	
	int iter = n * (l - 1);
	
	for (int z = 0; z < iter; z++) {
		
		FF(i, j) = 1;
		
		i = i + 1;
		j = j + 1;
	}
	
	mat temp = reshape(beta, n*l+ex, n);
	mat temp1 = temp.submat(0, 0, n*l - 1, n - 1).t();
	FF.submat(0, 0, n - 1, n*l - 1) = temp1;
	cx_vec eigval = eig_gen(FF);
	vec z = abs(eigval);
	double b = max(z);
	
	if (b > 1) {
		return(1);
	} else {
		return(0);
	}
}

double extract_diag(mat a, int pos) {
	vec b = a.diag();
	double b1 = as_scalar(b[pos - 1]);
	return(b1);
}

// [[Rcpp::export]]
cube calc_irf_cpp(mat Bt, mat Ht, int K, 
			      mat Fload, int sh_pos, 
				  int horizon) {
	
	int Ht_rows = Ht.n_rows;
	int Bt_rows = Bt.n_rows;
	int t = Ht_rows / K;
	int L = (Bt_rows - K) / (K*K);
	int start_r, end_r;
	mat Ht_submat(K, K, fill::zeros);
	mat choled_Ht(K, K, fill::zeros);
	mat start_sh(K, 1, fill::zeros);
	start_sh(sh_pos - 1, 0) = 1;
	int response_cols = Fload.n_rows;
	mat response(horizon, response_cols, fill::zeros);
	int beta_mat_dim = (Bt.n_rows - K) / K;
	mat beta_mat(beta_mat_dim, beta_mat_dim, fill::zeros);
	mat state_cont(K*L, 1, fill::zeros);
	cube response_cube(horizon, response_cols, t);
	
	for (int i = 1; i <= t; i++) {
		state_cont.zeros();
		start_r = (i - 1) * K;
		end_r = i * K - 1;
		
		Ht_submat = Ht.submat(start_r, 0, end_r, K - 1);
		choled_Ht = trans(chol(Ht_submat));
		
		beta_mat.submat(0, 0, K - 1, beta_mat_dim - 1) = extr_coef(Bt.col(i - 1), K, 2);
		beta_mat.submat(K, 0, beta_mat_dim - 1, beta_mat_dim - K - 1).diag().ones();
		
		state_cont.submat(0, 0, K - 1, 0) = choled_Ht * start_sh;
		response.row(0) = trans(Fload * state_cont.submat(0, 0, K - 1, 0));
		for (int j = 1; j < horizon; j++) {
			state_cont = beta_mat * state_cont;
			response.row(j) = trans(Fload * state_cont.submat(0, 0, K - 1, 0));
		}
		
		response_cube.slice(i - 1) = response;
	}
	return(response_cube);
}


// [[Rcpp::export]]
List calc_irf_mean_and_quantile(List irf_list,
								int t_,
								int horizon,
								int NN,
								vec quantile_vec) {
	int list_l = irf_list.length();
	
	vec V(list_l, fill::zeros);
	cube res_cube_mean(horizon, NN, t_);
	cube res_cube_quant1(horizon, NN, t_);
	cube res_cube_quant2(horizon, NN, t_);
	vec quant(2, fill::zeros);
		
	for (int j1 = 0; j1 < horizon; j1++) {
		Rcout << "Irf calculation: " << 100 * (j1 + 1) / horizon << "%\n";
		for (int j2 = 0; j2 < NN; j2++) {
			for (int j3 = 0; j3 < t_; j3++) {
				
				for (int i = 0; i < list_l; i++) {
					V[i] = as<cube>(irf_list(i))(j1, j2, j3);
				}
				
				res_cube_mean(j1, j2, j3) = mean(V);
				quant = quantile(V, quantile_vec);
				res_cube_quant1(j1, j2, j3) = quant[0];
				res_cube_quant2(j1, j2, j3) = quant[1];
			}
		}
	}
	
	return(List::create(res_cube_mean, 
						res_cube_quant1,
						res_cube_quant2));
}

// [[Rcpp::export]]
arma::mat alphahelper(arma::mat y, arma::mat Z, arma::mat Btdraw) {

double M = y.n_rows;
double t = y.n_cols;
arma::mat yhat = arma::zeros(M,t);

for (int i = 1; i < (t+1); i++) {
  yhat.col(i-1) = y.col(i-1) - Z.rows((i-1)*M,i*M-1)*Btdraw.col(i-1);
}

return yhat;
}

// [[Rcpp::export]]
mat prepare_const(mat mu, mat ar_coef) {
	int nrows = ar_coef.n_rows;
	int ncols = ar_coef.n_cols;
	mat mu_new(nrows, ncols, fill::zeros);
	for (int i = 0; i < ncols; i++) {
		for (int j = 0; j < nrows; j++) {
			mu_new(j, i) = mu(i, 0) * (1 - ar_coef(j, i));
		}
	}
	return(mu_new);
}

// [[Rcpp::export]]
List fast_reg(mat y, mat x) {
	mat x_ext(x.n_rows, x.n_cols + 1, fill::ones);
	x_ext.submat(0, 1, x.n_rows - 1, x.n_cols) = x;
	mat beta = inv(trans(x_ext) * x_ext) * trans (x_ext) * y;
	return(List::create(beta.submat(1, 0, x.n_cols, 0), beta(0,0)));
}

// [[Rcpp::export]]
List AR_KFS(mat Y,
			mat X,
			mat Q,
			mat hlast,
			mat beta0,
			mat P00,
			int L,
			int check,
			int maxdraws,
			int EX) {
	
	int T = Y.n_rows;
	int N = Y.n_cols;
	int ns = beta0.n_cols;
	mat beta_tt(T, ns, fill::zeros);
	cube ptt(ns, ns, T, fill::zeros);
	mat beta11 = beta0;
	mat p11 = P00;
	mat p10;
	mat x;
	double R = 0;
	mat beta10;
	mat yhat;
	mat feta;
	mat K;
	mat eta;
	mat A;
	
	for (int i = 0; i < T; i++) {
		
		x = X.row(i);
		R = as_scalar(hlast.submat(i + 1, 0, i + 1, 0));
		
		beta10 = beta11;
		p10 = p11 + Q;
		
		yhat = x * beta10.t();
		
		eta = Y.submat(i, 0, i, N - 1) - yhat.t();
		feta = x * p10 * x.t() + R;
		
		K = p10 * x.t() * inv(feta);
		
		beta11 = trans(beta10.t() + K * eta.t());
		A = eye(ns, ns) - K * x;
		p11 = A * p10 * A.t() + K * R * K.t();
		
		ptt.slice(i) = p11;
		beta_tt.row(i) = beta11;
	}
	
	int problem = 0;
	int trys = 1;
	int subcheck = -1;
	mat beta2(T, ns, fill::zeros);
	mat wa(T, ns, fill::randn);
	mat error(T, N, fill::zeros);	
	mat roots(T, 1, fill::zeros);	
	int i = T - 1;
	
	
	mat iFptF;
	mat bm;
	mat pm;
	mat p00;
	mat pt;
	
//	string tmp;
	
//	ofstream myfile;
	while (subcheck < 0 && trys <= maxdraws) {	
//	tmp = to_string(trys);
//	myfile.open("test_dir\\" + tmp + ".txt");
	p00 = ptt.slice(i);
	
//	myfile << '1' << ' ';

	beta2.row(i) = beta_tt.row(i) + wa.row(i) * chol(p00);
//	myfile << '2' << ' ';
	error.row(i) = Y.row(i) - X.row(i) * beta2.row(i).t();
	
	
	for (int i = T - 2; i >= 0; i--) {
		
//		myfile << '3' << ' ';
		pt = ptt.slice(i);
//		myfile << '4' << ' ';
		iFptF = inv(pt + Q);
//		myfile << '5' << ' ';
		bm = beta_tt.row(i) + (pt * iFptF * (beta2.row(i + 1).t() - beta_tt.row(i).t())).t();
//		myfile << '6' << ' ';
		pm = pt - pt * iFptF * pt;
		
//		myfile << '7' << ' ';
		beta2.row(i) = bm + wa.row(i) * chol(pm);
//		myfile << '8' << ' ';
		error.row(i) = Y.row(i) - X.row(i) * beta2.row(i).t();
//		if (i == 3) {
//			Rcout << Y.row(i) << " " << X.row(i) << " " << beta2.row(i);
//		}
//		myfile << '9' << ' ';
		roots(i, 0) = stability(beta2.row(i).t(), N, L, EX);
//		myfile << "10" << '\n';
	}
//	myfile.close();

//	return(List::create(beta2));
//	Rcout << trys << " ";
	if (check != 1) {
		return(List::create(beta2, error, roots, problem));
	}
	
	if (as_scalar(sum(roots)) == 0) {
		return(List::create(beta2, error, trys, problem));
	} else if (trys == maxdraws) {
		problem = 1;
		return(List::create(beta2, error, trys, problem));
	}
		
	trys = trys + 1;
	}
}


// [[Rcpp::export]]
mat iwpq(int v, mat ixpx) {
	int k = ixpx.n_rows;
	mat z(v, k, fill::zeros);
	
	mat cixpx = trans(chol(ixpx));
	mat randmat(k, 1, fill::zeros);
	
	for (int i = 0; i < v; i++) {
		randmat.randn();
		z.row(i) = trans(cixpx * randmat);
	}
	mat out = inv(z.t() * z);
	return(out);
}


// [[Rcpp::export]]
mat do_svol(mat hlast, double g, double mubar, double sigmabar, mat errors) {
	
	int T = errors.n_rows;
	mat hnew(T + 1, 1, fill::zeros);
	double hlag = 0;
	int i = 0;
	double yt;
	double htrial;
	double hlead = hlast(i + 1);
	double ss = sigmabar * g / (g + sigmabar);
	double mu = ss * (mubar/sigmabar + log(hlead)/g);
	double u = 0;
	double lp1, lp0;
	double h = exp(mu + sqrt(ss) * randn<double>());
	double L1 = 1;
	double accept;
	hnew(i, 0) = h;
	
	for (int i = 1; i < T; i++) {
		hlead = hlast(i + 1, 0);
		hlag = hnew(i - 1, 0);
		yt = errors(i - 1, 0);
		
		mu = (log(hlead) + log(hlag))/2;
		ss = g / 2;
		
		htrial = exp(mu + sqrt(ss) * randn<double>());
		
		lp1 = -0.5*log(htrial) - (yt * yt)/(2*htrial);
		lp0 = -0.5*log(hlast(i)) - (yt * yt)/(2*hlast(i));
		accept = min_sc(L1, exp(lp1 - lp0));
		
		u = randu<double>();
		
		if (u <= accept) {
			h = htrial;
		} else {
			h = hlast(i, 0);
		}
		
		hnew(i, 0) = h;
	}
	
	i = T;
	
	yt = errors(i - 1, 0);
	hlag = hnew(i - 1, 0);
	
	mu = log(hlag);
	ss = g;
	
	htrial = exp(mu + sqrt(ss) * randn<double>());
	
	lp1 = -0.5 * log(htrial) - (yt * yt)/(2 * htrial);
	lp0 = -0.5 * log(hlast(i)) - (yt * yt)/(2 * hlast(i));
	accept = min_sc(1, exp(lp1 - lp0));
	
	u = randu<double>();
	
	if (u <= accept) {
		h = htrial;
	} else {
		h = hlast(i, 0);
	}
	
	hnew(i, 0) = h;
	
	return(hnew);
}

// [[Rcpp::export]]
mat getreg(mat Y, 
		   mat X,
		   mat B0,
		   mat Sigma0,
		   double sigma2) {
	
	int colx = X.n_cols;
	mat V = inv(inv(Sigma0) + (1/sigma2) * (X.t() * X));
	mat M = V * (inv(Sigma0)*B0 + (1/sigma2) * X.t() * Y);
	mat rand_mat(1, colx, fill::randn);
	mat B = M + trans(rand_mat * chol(V));
	return(B);
}

// [[Rcpp::export]]
double IG(double v0, 
		  double d0,
		  mat x) {
	
	double T = x.n_rows;
	double v1 = v0 + T;
	double d1 = d0 + as_scalar(trans(x) * x);
	mat z(v1, 1, fill::randn);
	double x1 = as_scalar(trans(z) * z);
	double v = d1 / x1;
	return(v);
}


// [[Rcpp::export]]
List preparex(mat a_mat, int L, int Const) {
	mat Y = a_mat;
	
	int c_int = 0;
	if (Const == 1) {
		c_int = 1;
	}
	
	int colsn = a_mat.n_cols * L + c_int;
	int rowsn = a_mat.n_rows;
	
	mat X(rowsn, colsn, fill::zeros);
	mat lag_mat;
	
	int acc_c_l = 0;
	int acc_c_r = 0;
	for (int i = 1; i <= L; i++) {
		acc_c_l = a_mat.n_cols * (i - 1);
		acc_c_r = a_mat.n_cols * i - 1;
		
		X.submat(0, acc_c_l, rowsn - 1, acc_c_r) = lag0(a_mat, i);
	}
	
	if (Const == 1) {
		X.col(colsn - 1).fill(1);
	}
	
	mat Y_hat = Y.submat(L, 0, rowsn - 1, a_mat.n_cols - 1);
	mat X_hat = X.submat(L, 0, rowsn - 1, colsn - 1);
	return(List::create(Y_hat, X_hat));
}

// [[Rcpp::export]]
mat prepare_for_factors_CK(mat data,
						   mat beta2e,
						   int Lx
						   ) {
	int T = data.n_rows;
	int NN = data.n_cols;
	
	mat dataF(T, NN, fill::zeros);
	
	for (int i = 0; i < NN; i++) {
		dataF.col(i) = remSC(data.submat(0, i, T - 1, i), beta2e.submat(0, i, T - 1, i));
	}
	
	dataF.submat(0, 0, Lx - 1, NN - 1) = repmat_cpp(dataF.submat(Lx, 0, Lx, NN - 1), Lx, 1);
	
	return(dataF);
}

// [[Rcpp::export]]
mat test_mat_ret(mat a, 
				 int Lx) {
	int NN = a.n_cols;
	mat b = a.submat(0, 0, Lx - 1, NN - 1);
	return(b);
}

// [[Rcpp::export]]
List Factors_KFS(
				 mat dataF,
				 mat fload,
				 mat beta2e,
				 int Lx,
				 int L,
				 mat pmat00,
				 mat vmat00,
				 List AR_fact_list,
				 int K,
				 mat hlaste,
				 mat hlast,
				 List a_ij_list
				 ) {
	int T = dataF.n_rows;
	int NN = dataF.n_cols;

	int ns = pmat00.n_cols;
	mat beta_tt(T, ns, fill::zeros);
	cube ptt(ns, ns, T);
	mat beta11 = pmat00;
	mat p11 = vmat00;
	mat p00(K, K, fill::zeros);
	mat rho;
	mat H(NN, K*2, fill::zeros);
	H.submat(0, 0, NN - 1, K - 1) = fload;
	mat R(NN, NN, fill::zeros);
	mat Q(ns, ns, fill::zeros);
	mat F(ns, ns, fill::zeros);
	mat tmp_mat;
	mat mu(1, ns, fill::zeros);
	mat x;
	mat K_gain;
	mat f_ad, q_ad, mu_ad;
	mat ifptfq;
	mat bm, pm, beta10, p10, yhat, eta, feta, A, pt;
	mat A_mat(K, K, fill::zeros);
	A_mat.diag().ones();	
	mat chol_t;
	mat Q_bf;
	
	for (int i = 0; i < T; i++) {
		Rcout << " " << 1;
		// Observation equation
		rho = trans(beta2e.row(i));
		mat H2 = -repmat(rho, 1, K) % fload;
		H.submat(0, K, NN - 1, K*2 - 1) = H2;
		
		Rcout << " " << 2;
		R.diag() = hlaste.row(i + 1);
		
		// Transition equation
		Q.zeros();
		Q.submat(0, 0, K - 1, K - 1).diag() = hlast.row(i + 1);
		
		Rcout << " " << 3;
		for (int j = 1; j < K; j++) {
			
			A_mat.submat(j, 0, j, j - 1) = as<mat>(a_ij_list[j]).row(i);
			
		}

		Q.submat(0, 0, K - 1, K - 1) = inv(A_mat) * Q.submat(0, 0, K - 1, K - 1) * trans(inv(A_mat));
		Q_bf = Q;
		Rcout << " " << 4;
		for (int j = 0; j < K; j++) {
			tmp_mat = as<mat>(AR_fact_list[j]);
			F.submat(j, 0, j, ns - 1) = tmp_mat.submat(i, 0, i, ns - 1);
			mu(0, j) = tmp_mat(i, ns);
		}
		
		Rcout << " " << "4_";
		Rcout << " " << "4_";

		F.submat(K, 0, ns - 1, ns - K - 1).diag().ones();
		
		Rcout << " " << 5;
		x = H;
		beta10 = mu + beta11 * F.t();
		p10 = F * p11 * F.t() + Q;
		
		try {
			chol_t = chol(p10);
		} catch (...) {
			Rcout << 'z';
			Rcout << ' ' << i << ' ';
			return(List::create(1, A_mat, Q_bf));
		}
		
		Rcout << " " << 6;
		yhat = trans(x * beta10.t());
		
		eta = dataF.submat(i, 0, i, NN - 1) - yhat;
		
		feta = x * p10 * x.t() + R;
		
		K_gain = p10 * x.t() * inv(feta);
		
		beta11 = trans(beta10.t() + K_gain * eta.t());
		
		Rcout << " " << 7;
		A = eye(ns, ns) - K_gain * x;
		
		p11 = A * p10 * A.t() + K_gain * R * K_gain.t();
		p11 = 0.5 * (p11.t() + p11);
		ptt.slice(i) = p11;
		beta_tt.row(i) = beta11;
	}
	
	mat beta2(T, ns, fill::zeros);
	mat wa(T, ns, fill::randn);
	
	int i = T - 1;
	Rcout << 'a';
	p00 = ptt.slice(i).submat(0, 0, K - 1, K - 1);
	beta2.row(i) = beta_tt.row(i);
	beta2.submat(i, 0, i, K - 1) = beta_tt.submat(i, 0, i, K - 1) + wa.submat(i, 0, i, K - 1) * chol(p00);
	Rcout << 'b';
	for (int i = T - 2; i >= 0; i--) {
		Q.zeros();
		Q.submat(0, 0, K - 1, K - 1).diag() = hlast.row(i + 2);
		
		for (int j = 1; j < K; j++) {
			
			A_mat.submat(j, 0, j, j - 1) = as<mat>(a_ij_list[j]).row(i + 1);
		}
		
		Q.submat(0, 0, K - 1, K - 1) = inv(A_mat) * Q.submat(0, 0, K - 1, K - 1) * trans(inv(A_mat));
		
		for (int j = 0; j < K; j++) {
			tmp_mat = as<mat>(AR_fact_list[j]);
			F.submat(j, 0, j, ns - 1) = tmp_mat.submat(i + 1, 0, i + 1, ns - 1);
			mu(0, j) = tmp_mat(i + 1, ns);
		}
		
	F.submat(K, 0, ns - 1, ns - K - 1).diag().ones();
	
	f_ad = F.submat(0, 0, K - 1, ns - 1);
	q_ad = Q.submat(0, 0, K - 1, K - 1);
	mu_ad = mu.submat(0, 0, 0, K - 1);
	
	pt = ptt.slice(i);
	ifptfq = inv(f_ad * pt * f_ad.t() + q_ad);
	
	bm = beta_tt.row(i) + trans(pt * f_ad.t() * ifptfq * trans(
		beta2.submat(i + 1, 0, i + 1, K - 1) - mu_ad - beta_tt.row(i) * f_ad.t())
		);
		
	pm = pt - pt * f_ad.t() * ifptfq * f_ad * pt;
	beta2.row(i) = bm;
	
	try{
	beta2.submat(i, 0, i, K - 1) = bm.submat(0, 0, 0, K - 1) + wa.submat(i, 0, i, K - 1) * chol(pm.submat(0, 0, K - 1, K - 1));
	} catch(...) {
		return(List::create(1, pt));
	}
	}
	Rcout << 'c';

	return(List::create(beta2));
}


// [[Rcpp::export]]
mat Innov_cor_mat_KFS(mat beta0,
					  mat p00,
					  mat hlast,
					  mat Q,
					  mat Y,
					  mat X) {
	
	int ns = beta0.n_cols;
	int T_obs = X.n_rows;
	mat beta_tt(T_obs, ns, fill::zeros);
	
	cube ptt(ns, ns, T_obs);
	mat F(ns, ns, fill::zeros);
	F.diag().ones();
	
	mat yhat(1, 1, fill::zeros);
	mat eta(1, 1, fill::zeros);
	mat feta(1, 1, fill::zeros);
	mat p10(ns, ns, fill::zeros);
	mat beta10(1, ns, fill::zeros);
	mat beta11 = beta0;
	mat p11 = p00;
	mat x(1, ns, fill::zeros);
	mat R(1, 1, fill::zeros);
	mat K;
	mat pt(ns, ns, fill::zeros);
	mat A;
	
	for (int i = 0; i < T_obs; i++) {
		x = X.row(i);
		R(0, 0) = hlast(i + 1, 0);
		
		beta10 = beta11 * F.t();
		p10 = F * p11 * F.t() + Q;
		
		yhat = trans(x * beta10.t());
		eta = Y.row(i) - yhat;
		feta = x * p10 * x.t() + R;
		
		K = p10 * x.t() * inv(feta);
		
		beta11 = trans(beta10.t() + K * eta.t());
		
		A = eye(ns, ns) - K * x;
		
		p11 = A * p10 * A.t() + K * R * K.t();
		
		ptt.slice(i) = p11;
		beta_tt.row(i) = beta11;
	}
	
	mat beta2(T_obs, ns, fill::zeros);
	mat wa(T_obs, ns, fill::randn);
	mat bm;
	mat pm;
	
	int i = T_obs - 1;
	p00 = ptt.slice(i);
	
	beta2.row(i) = beta_tt.row(i) + wa.row(i) * chol(p00);
	
	for (int i = T_obs - 2; i >= 0; i--) {
		pt = ptt.slice(i);
		
		bm = beta_tt.row(i) + trans(
		pt * F.t() * inv(F * pt * F.t() + Q) * trans(beta2.row(i + 1) - beta_tt.row(i) * F.t()) 
		);
		pm = pt - pt * F.t() * inv(F * pt * F.t() + Q) * F * pt;
		
		beta2.row(i) = bm + wa.row(i) * chol(pm);
	}
	
	return(beta2);
}

// [[Rcpp::export]]
mat ext_vareps(mat fact_AR_error, 
			   List a_ij_list,
			   int K) {
	
	mat A(K, K, fill::zeros);
	A.diag().ones();
	
	int s1 = fact_AR_error.n_rows;
	int s2 = fact_AR_error.n_cols;
	
	mat varepsilon(s1, s2, fill::zeros);
	
	for (int i = 0; i < s1; i++) {
		
		for (int j = 1; j < K; j++) {
			
			A.submat(j, 0, j, j - 1) = as<mat>(a_ij_list[j]).row(i);
		}
	
	varepsilon.row(i) = fact_AR_error.row(i) * A.t();
	}
	
	return(varepsilon);
}

// [[Rcpp::export]]
arma::colvec mvndrawC(arma::colvec mu, arma::mat sig) {

double k = mu.size();
arma::colvec aux = as<arma::colvec>(rnorm(k));
arma::mat csig = arma::chol(sig).t();
arma::colvec out = mu + csig*aux;
return(out);

}

// [[Rcpp::export]]
mat test_chol(mat a) {
	mat b = chol(a);
	return(b);
}

// [[Rcpp::export]]
mat test_list(List a) {
	return(a["a_mat"]);
}

 // [[Rcpp::export]]
List carterkohn(arma::mat y, arma::mat Z, arma::mat Ht, arma::mat Qt, 
				double m, double p, double t, arma::colvec B0, arma::mat V0) {
 
arma::colvec bp = B0; // reuses memory and avoids extra copy
int ns = B0.n_rows;
arma::mat Vp = V0;
arma::colvec btt = B0;	// initialize now s/t that their scope extends beyond loop
arma::mat Vtt = V0;
arma::mat f = arma::zeros(p, Ht.n_cols);
arma::mat invf = f;
arma::colvec cfe = arma::zeros(y.n_rows);
arma::mat bt = arma::zeros(t,m);
arma::mat Vt = arma::zeros(pow(m,2),t);
arma::colvec loglik = arma::zeros(1);
arma::mat R = arma::zeros(p,p);
arma::mat H = arma::zeros(p,m);
arma::mat A, K_gain;
List bbb = List::create();
List cov_l = List::create();

for (int i = 1; i < (t+1); i++) {
  R = Ht.rows((i - 1) * p, i * p - 1);
  H = Z.rows((i - 1) * p, i * p - 1);

  cfe = y.col(i-1) - H*bp;   // conditional forecast error
  f = H*Vp*H.t() + R;    	  // variance of the conditional forecast error
  invf = f.i();		  // invert only once
  loglik = loglik + log(det(f)) + cfe.t()*invf*cfe;
  
  bbb.push_back(cfe);
  btt = bp + Vp*H.t()*invf*cfe;
  K_gain = Vp*H.t()*invf;
  A = eye(ns, ns) - K_gain*H;
  Vtt = A*Vp*A.t() + K_gain * R * K_gain.t();

  if (i < t) {
    bp = btt;
    Vp = Vtt + Qt;
  }
  cov_l.push_back(Vtt);
  bt.row(i-1) = btt.t();
  Vt.col(i-1) = arma::vectorise(Vtt);
}

arma::mat bdraw = arma::zeros(t,m);
bdraw.row(t-1) = mvndrawC(btt,Vtt).t();

for (int i = 1; i < t; i++) {	// Backward recursions
    arma::colvec bf = bdraw.row(t-i).t();
    btt = bt.row(t-i-1).t();
    Vtt = arma::reshape(Vt.col(t-i-1),m,m);
    f = Vtt + Qt;
    invf = f.i();
    cfe = bf - btt;
    arma::colvec bmean = btt + Vtt*invf*cfe;
    arma::mat bvar = Vtt - Vtt*invf*Vtt;
    bdraw.row(t-i-1) = mvndrawC(bmean, bvar).t();
}

return List::create(Named("loglik") = loglik, 
					Named("bdraws") = bdraw.t(),
					bbb, bt, cov_l);
}

// [[Rcpp::export]]
arma::mat sigmahelper1(arma::mat Atdraw, double M) {
  
double t = Atdraw.n_cols; 
  
arma::mat capAt = arma::zeros(M*t,M);
for (int i = 1; i < (t+1); i++) {
  arma::mat capatemp(M,M);
  capatemp.eye();
  arma::colvec aatemp = Atdraw.col(i-1);
  double ic = 1;
  for (int j = 2; j < (M+1); j++) {
    capatemp(arma::span(j - 1, j - 1), arma::span(0,j - 2)) = aatemp.rows(ic - 1, ic + j - 3).t();
    ic = ic + j - 1;
  }
  capAt.rows((i - 1) * M,(i * M) - 1) = capatemp;
}

return capAt;
}

// [[Rcpp::export]]
List sigmahelper2(arma::mat capAt, arma::mat yhat, arma::colvec qs, arma::colvec ms, 
				  arma::colvec u2s, arma::mat Sigtdraw, arma::mat Zs, arma::mat Wdraw, 
				  arma::colvec sigma_prmean, arma::mat sigma_prvar) {

double M = capAt.n_cols;
double t = capAt.n_rows / M; 

arma::mat y2 = arma::zeros(M,t);
for (int i = 1; i < (t+1); i++) {
  y2.col(i-1) = pow(capAt.rows((i-1)*M,(i*M)-1) * yhat.col(i-1),2);
}

arma::mat aux = 0.001 * arma::ones(t,M);
arma::mat yss = log( aux + y2.t() );

arma::colvec cprw = arma::zeros(7,1);
arma::mat statedraw = arma::zeros(t,M);
for (int jj = 1; jj < (M+1); jj++) {
  for (int i = 1; i < (t+1); i++) {
  arma::colvec prw = arma::zeros(7,1);
    for (int k = 1; k < 8; k++) {
      prw(k-1) = qs(k-1) * (1/sqrt(2*M_PI*u2s(k-1)))*exp(-0.5*((pow(yss(i-1,jj-1) - Sigtdraw(jj-1,i-1) - ms(k-1) + 1.2704,2))/u2s(k-1)));
    }
  cprw = arma::cumsum(prw/arma::sum(prw));
  double trand = as<double>(runif(1));
  double imix = 0;
  if (trand < cprw[0]){
    imix = 1;
  } else if (trand < cprw[1]) {
    imix = 2;
  } else if (trand < cprw[2]) {
    imix = 3;
  } else if (trand < cprw[3]) {
    imix = 4;
  } else if (trand < cprw[4]) {
    imix = 5;
  } else if (trand < cprw[5]) {
    imix = 6;
  } else if (trand < cprw[6]) {
    imix = 7;
  }
  statedraw(i-1,jj-1) = imix;  
  }
}

arma::mat vart = arma::zeros(t*M,M);
arma::mat yss1 = arma::zeros(t,M);
for (int i = 1; i < (t+1); i++) {
  for (int j = 1; j < (M+1); j++) {
    double imix = statedraw(i-1,j-1);
    vart(((i-1)*M+j-1),j-1) = u2s(imix-1);
    yss1(i-1,j-1) = yss(i-1,j-1) - ms(imix-1) + 1.2704;
  }
}

arma::mat Sigtdraw_new = carterkohn(yss1.t(),Zs,vart,Wdraw,M,M,t,sigma_prmean,sigma_prvar)["bdraws"];

arma::mat sigt = arma::zeros(M*t,M);
for (int i = 1; i < (t+1); i++) {
  arma::mat sigtemp = arma::zeros(M,M);
  sigtemp.diag() = exp(0.5*Sigtdraw_new.col(i-1));
  sigt.rows((i-1)*M,(i*M)-1) = sigtemp;
}

return List::create(Named("Sigtdraw") = Sigtdraw_new, Named("sigt") = sigt);
}

// [[Rcpp::export]]
List sigmahelper3(arma::mat capAt, arma::mat sigt){

double M = sigt.n_cols;
double t = sigt.n_rows/M;
arma::mat Ht = arma::zeros(M*t,M);
arma::mat Htsd = arma::zeros(M*t,M);

for (int i = 1; i < (t+1); i++) {
  arma::mat inva = capAt.rows((i-1)*M, (i*M)-1).i();
  arma::mat stem = sigt.rows((i-1)*M, (i*M)-1);
  arma::mat Hsd = inva*stem;
  Ht.rows((i-1)*M, (i*M)-1) = Hsd * Hsd.t();
  Htsd.rows((i-1)*M, (i*M)-1) = Hsd;
}

return List::create(Named("Ht") = Ht, Named("Htsd") = Htsd);
}


// [[Rcpp::export]]
List Factors_KFS1(
				  mat dataF,
				  mat fload,
				  mat beta2e,
				  int Lx,
				  int L,
				  mat pmat00,
				  mat vmat00,
				  mat Ar_trans_mat,
				  int K,
				  mat hlaste,
				  mat Ht,
				  int s
				  ) {
	int T = dataF.n_rows;
	int NN = dataF.n_cols;

	int ns = pmat00.n_cols;
	mat beta_tt(T, ns, fill::zeros);
	cube ptt(ns, ns, T);
	mat beta11 = pmat00;
	mat p11 = vmat00;
	mat p00(K, K, fill::zeros);
	mat rho;
	mat H(NN, ns, fill::zeros);
	H.submat(0, 0, NN - 1, K - 1) = fload;
	mat R(NN, NN, fill::zeros);
	mat Q(ns, ns, fill::zeros);
	mat F(ns, ns, fill::zeros);
	mat tmp_mat;
	mat mu(1, ns, fill::zeros);
	mat x;
	mat K_gain;
	mat f_ad, q_ad, mu_ad;
	mat ifptfq;
	mat bm, pm, beta10, p10, yhat, eta, feta, A, pt;
	mat A_mat(K, K, fill::zeros);
	A_mat.diag().ones();	
	mat chol_t;
	mat Q_bf;
	mat z_out;
	int z_col, Ht_row1, Ht_row2;
	
	for (int i = 0; i < T; i++) {
		//Rcout << " " << 1;
		// Observation equation
		rho = trans(beta2e.row(i));
		mat H2 = -repmat(rho, 1, K) % fload;
		H.submat(0, K, NN - 1, K * 2 - 1) = H2;
		R.diag() = hlaste.row(i + 1);
		
		// Transition equation
		Q.zeros();
		
		//Под замену!!!!!!!!!!!!!!!!!!!!
		// Done
		Ht_row1 = i * K;
		Ht_row2 = (i + 1) * K - 1;
		
		//Под замену!!!!!!!!!!!!!!!!!!!!
		// Done
		//Rcout << " " << "2_1";
		Q.submat(0, 0, K - 1, K - 1) = Ht.submat(Ht_row1, 0, Ht_row2, K - 1);
		
		//Rcout << " " << 3;
		Q_bf = Q;
		//Rcout << " " << 4;
		//Под замену!!!!!!!!!!!!!!!!!!!!
		//Rcout << " " << "4_1";
		mu.submat(0, 0, 0, K - 1) = trans(extr_coef(Ar_trans_mat.col(i), K, 1));
		//Rcout << " " << "4_2";
		z_out = extr_coef(Ar_trans_mat.col(i), K, 2);
		//Rcout << " " << "4_3";
		z_col = z_out.n_cols;
		//Rcout << " " << "4_4";
		F.submat(0, 0, K - 1, z_col - 1) = z_out;
		
		//Rcout << " " << "4_5";
		F.submat(K, 0, ns - 1, ns - K - 1).diag().ones();
		
		//Rcout << " " << 5;
		x = H;
		beta10 = mu + beta11 * F.t();
		p10 = F * p11 * F.t() + Q;
		
		//try {
		//	chol_t = chol(p10);
		//} catch (...) {
		//	Rcout << 'z';
		//	Rcout << ' ' << i << ' ';
		//	return(List::create(1, p10, Q_bf));
		//}
		
		//Rcout << " " << 6;
		yhat = trans(x * beta10.t());
		
		eta = dataF.submat(i, 0, i, NN - 1) - yhat;
		
		feta = x * p10 * x.t() + R;
		
		K_gain = p10 * x.t() * inv(feta);
		
		beta11 = trans(beta10.t() + K_gain * eta.t());
		
		//Rcout << " " << 7;
		A = eye(ns, ns) - K_gain * x;
		
		p11 = A * p10 * A.t() + K_gain * R * K_gain.t();
		p11 = 0.5 * (p11.t() + p11);
		ptt.slice(i) = p11;
		beta_tt.row(i) = beta11;
	}
	
	mat beta2(T, ns, fill::zeros);
	mat wa(T, ns, fill::randn);
	
	int i = T - 1;
	//Rcout << 'a';
	p00 = ptt.slice(i).submat(0, 0, K - 1, K - 1);
	beta2.row(i) = beta_tt.row(i);
	beta2.submat(i, 0, i, K - 1) = beta_tt.submat(i, 0, i, K - 1) + wa.submat(i, 0, i, K - 1) * chol(p00);
	//Rcout << 'b';
	for (int i = T - 2; i >= 0; i--) {
	Q.zeros();
	Ht_row1 = i * K;
	Ht_row2 = (i + 1) * K - 1;
	Q.submat(0, 0, K - 1, K - 1) = Ht.submat(Ht_row1, 0, Ht_row2, K - 1);
	
	mu.submat(0, 0, 0, K - 1) = trans(extr_coef(Ar_trans_mat.col(i), K, 1));
	z_out = extr_coef(Ar_trans_mat.col(i), K, 2);
	F.submat(0, 0, K - 1, z_col - 1) = z_out;
		
	F.submat(K, 0, ns - 1, ns - K - 1).diag().ones();
	
	f_ad = F.submat(0, 0, K - 1, ns - 1);
	q_ad = Q.submat(0, 0, K - 1, K - 1);
	mu_ad = mu.submat(0, 0, 0, K - 1);
	
	pt = ptt.slice(i);
	ifptfq = inv(f_ad * pt * f_ad.t() + q_ad);
	
	bm = beta_tt.row(i) + trans(pt * f_ad.t() * ifptfq * trans(
		beta2.submat(i + 1, 0, i + 1, K - 1) - mu_ad - beta_tt.row(i) * f_ad.t())
		);
		
	pm = pt - pt * f_ad.t() * ifptfq * f_ad * pt;
	beta2.row(i) = bm;
	
	try{
	beta2.submat(i, 0, i, K - 1) = bm.submat(0, 0, 0, K - 1) + wa.submat(i, 0, i, K - 1) * chol(pm.submat(0, 0, K - 1, K - 1));
	} catch(...) {
		return(List::create(1, pt));
	}
	}
	//Rcout << 'c';

	return(List::create(beta2));
}