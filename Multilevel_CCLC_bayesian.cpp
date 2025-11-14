// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include </opt/homebrew/Cellar/boost/1.81.0_1/include/boost/random.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <armadillo>
#include</opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/unsupported/Eigen/CXX11/Tensor>
#define _USE_MATH_DEFINES
#include <chrono>


using namespace Eigen;
using namespace std;
using namespace arma;

int fattoriale(int a){
    int f=1;
    if(a==0){
        f=1;
    }
    if(a>0){
        for(int i=0;i<a;i++){
            f=f*(a-i);
        }
    }
    return f;
}


// sample from dirichlet distribution
arma::vec rdirichlet(arma::vec alpha_m)
{
    int distribution_size = alpha_m.n_elem;
    arma::vec distribution = arma::zeros(distribution_size);

    double sum_term = 0 ;
    // draw Gamma variables
    for (int j = 0; j < distribution_size; j++){

        std::random_device rd;
        std::mt19937 gen(rd());
        std::gamma_distribution<> d(alpha_m[j], 1.0);
        double cur = d(gen);

        distribution(j) = cur;
        sum_term += cur;
    }
    // normalize
    for (int j = 0; j < distribution_size; j++){
        distribution(j) = distribution(j)/sum_term;
    }
    return(distribution);
}


tuple <mat,int,int> vector_to_matrix(const std::vector< std::vector<int> > &student_vector)
{
  mat y(student_vector.size(),student_vector[0].size());
  
  for(int i = 0; i < student_vector.size(); i++)
  {
    for(int j = 0; j < student_vector[i].size(); j++)
    {
      y(i,j) = student_vector[i][j];
    }
  }
  
  int _n=student_vector.size();
  int _G=student_vector[0].size();
  
  return {y,_n,_G};
}


// [[Rcpp::export]]
int main2(vector<vector<int>> yy){
  mat y;
  int _n;
  int _G;
  tie(y,_n,_G)=vector_to_matrix(yy);
  
    auto start = chrono::steady_clock::now();
    std::vector< std::vector<int> > student_vector;

    int _L, _H, _R, _K, _Q, _iter, _iter0, _iter2, _iter3, _iter4, _starting, _starting2, _burn, _thinn;
    _L = 4;
    _H = 3;
    _R = 2;
    _iter = 150;
    _iter0 = 1000;
    _iter2 = 100;
    _iter3 = 100;
    _iter4 = 0;
    _starting = 5;
    _starting2 = 15;
    _burn = 50;
    _thinn = 3;
    double eps0 = 0.0000000000001;
    double eps4 = 0.000001;

    _K = y(0, 1);
    for (int j = 0; j < _n; j++) {
        if (y(j, 1) > _K) {
            _K = y(j, 1);
        }
    }
    _Q = y(0, 2);
    for (int j = 0; j < _n; j++) {
        if (y(j, 2) > _Q) {
            _Q = y(j, 2);
        }
    }

    //thinning
    int nn = 0;
    int nnn = _iter - _burn;
    while (_thinn <= nnn) {
        nn = nn + 1;
        nnn = nnn - _thinn;
    }
    cout << nn << endl;
    vector<int> it(nn);
    it[0] = _burn + _thinn - 1;
    for (int j = 1; j < nn; j++) {
        it[j] = it[j - 1] + _thinn;
    }
    int mu = 0;

    vector<int> _vettoreqk(_n);
    int xx = 0;
    for (int q = 0; q < _Q; q++) {
        for (int j = 0; j < _n; j++) {
            if (y(j, 2) == q + 1) {
                _vettoreqk[xx] = y(j, 0);
                xx = xx + 1;
            }
        }
    }

    vector<int> _indvecqk(_n);
    int xxx = 0;
    for (int j = 0; j < _n; j++) {
        if (y(j, 0) == _vettoreqk[xxx]) {
            _indvecqk[xxx] = j;
            xxx = xxx + 1;
            if (xxx <= _n - 1) {
                j = -1;
            } else if (xxx == _n) {
                j = _n - 1;
            }
        }
    }


    vector<int> _nk(_K);
    vector<int> _nq(_Q);
    mat _nkq(_K, _Q);

    for (int k = 0; k < _K; k++) {
        _nk[k] = 0;
    }
    for (int k = 0; k < _K; k++) {
        for (int j = 0; j < _n; j++) {
            if (y(j, 1) == k + 1) {
                _nk[k] += 1;
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        _nq[q] = 0;
    }
    for (int q = 0; q < _Q; q++) {
        for (int j = 0; j < _n; j++) {
            if (y(j, 2) == q + 1) {
                _nq[q] += 1;
            }
        }
    }

    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            _nkq(k, q) = 0;
        }
    }

    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int j = 0; j < _n; j++) {
                if (y(j, 1) == k + 1 && y(j, 2) == q + 1) {
                    _nkq(k, q) += 1;
                }
            }
        }
    }

    int ckq = 0;
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            if (_nkq(k, q) != 0) {
                ckq += 1;
            }
        }
    }

    arma::cube _pi_l(_L, _H, _R);
    for (int h = 0; h < _H; h++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pi_l(l, h, r) = (double) 1 / _L;
            }
        }
    }


    //local dependence
    int nlocdep = 0;
    mat locdep(nlocdep, 2);

    int _I = _G - 3;
    vector<int> _vC(_I);
    int ma;
    for (int i = 0; i < _I; i++) {
        ma = y(0, i + 3);
        for (int j = 0; j < _n; j++) {
            if (y(j, i + 3) > ma) {
                ma = y(j, i + 3);
            }
        }
        _vC[i] = ma + 1;
    }

    int C;
    int maa;
    maa = _vC[0];
    for (int i = 0; i < _I; i++) {
        if (_vC[i] > maa) {
            maa = _vC[i];
        }
    }
    C = maa;


    Tensor<double, 5> _p_finale(_I, C, _L, _H, _R);
    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _p_finale(i, c, l, h, r) = 0;
                    }
                }
            }
        }
    }

    Tensor<double, 6> _p_finale2(_starting2, _I, C, _L, _H, _R);
    for (int st2 = 0; st2 < _starting2; st2++) {
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _p_finale2(st2, i, c, l, h, r) = 0;
                        }
                    }
                }
            }
        }
    }

    mat _Prob1(_K, _H);
    mat _Prob2(_Q, _R);
    //boost::mt19937 generator2(1234);
    std::random_device rdtest;
    std::mt19937 generator(rdtest());
    mat _w(_K, _Q);
    mat _z(_Q, _R);
    vector<int> _numh(_H);
    vector<double> _p_h(_H);
    vector<int> _numr(_R);
    vector<double> _p_r(_R);

    mat _w_hat(_K, _H);
    mat _z_hat(_Q, _R);
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            _w_hat(k, h) = 0;
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _z_hat(q, r) = 0;
        }
    }


    Tensor<double, 4> _f(_n, _L, _H, _R);
    //Tensor<double,4> _f(_n,_L,_H,_R);
    Tensor<double, 4> _res(_K, _Q, _H, _R);
    mat _ff(_K, _H);
    mat _num(_n, _L);
    vector<double> _den(_n);
    vector<double> _denom(_n);
    mat _Prob3(_n, _L);
    mat _num2(_K, _H);
    vector<double> _den2(_K);
    mat _ff2(_Q, _R);
    Tensor<double, 4> _res2(_Q, _K, _H, _R);
    mat _num3(_Q, _R);
    vector<double> _den3(_Q);
    mat _wjk(_n, _H);//estensione w in cui ripeto la k-sima riga di w _nk[k] volte
    mat _zjq(_n, _R);//estensione z in cui ripeto la q-sima riga di z _nq[q] volte
    //vector<int> _copiank(_K);
    //vector<int> _copianq(_Q);
    //Tensor<int,3> _A(_H,_n,_R);
    arma::cube _B(_n, _H, _R);
    mat _den4(_H, _R);
    //Tensor<double,4> _aa(_n,_L,_H,_R);
    arma::cube _num4(_L, _H, _R);
    vector<double> _den5(_L);
    mat _den6(_L, _R);
    mat _den7(_L, _H);
    //arma::cube _den5(_L,_H,_R);
    //mat _num5(_I,_L);

    mat _Prob4l(_n,_L);
    mat _Prob4lvar(_n,_L);
    Tensor<double,4> _Pxwzyvar(_n,_L,_H,_R);

    mat _PL(_starting, _L);
    arma::cube _Tau(_starting, _n, _L);
    arma::cube _F(_starting, _n, _L);
    //vector<int> Seed(_starting);
    vector<double> logLik(_starting);
    Tensor<double, 4> _pf(_n, _L, _H, _R);
    arma::cube _fykq(_n, _H, _R);
    Tensor<double, 4> _res3(_n, _L, _H, _R);
    mat _xjkq(_n, _L);
    mat _tt(_K, _Q);
    double min1, temp1;
    double min2;

    mat _x_hat(_n, _L);
    _x_hat.zeros(_n, _L);

    mat _x_hat_int(_n, _L);
    _x_hat_int.zeros(_n, _L);

    vector<int> temp2(_K);
    mat temp3(_n, _R);
    arma::cube temp4(_n, _L, _R);

    vector<int> _indh(_K);
    vector<int> _indr(_Q);
    mat _classif2(_n, 3);

    vector<int> _indl(_n);

    mat _classif(_n, 4);

    double _maxL;

    Tensor<double, 4> _pi_l_hat(_starting2, _L, _H, _R);
    //Eigen::Tensor<double, 4, 0>::operator()<double,int> _pi_l_hat(_starting2, _L, _H, _R);
    //in instantiation of function template specialization 'Eigen::Tensor<double, 4, 0>::operator()<double, int>' requested here
    for (int st2 = 0; st2 < _starting2; st2++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _pi_l_hat(st2, l, h, r) = 0;
                }
            }
        }
    }
    arma::cube _pi_lc_hat(_L, _H, _R);
    for (int l = 0; l < _L; l++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _pi_lc_hat(l, h, r) = 0;
            }
        }
    }

    mat _p_h_hat(_starting2, _H);
    for (int st2 = 0; st2 < _starting2; st2++) {
        for (int h = 0; h < _H; h++) {
            _p_h_hat(st2, h) = 0;
        }
    }
    mat _p_hc_hat(_starting2,_H);
    vector<double> _p_hcc_hat(_H);
    //std::fill(_p_hc_hat.begin(), _p_hc_hat.end(), 0);
    mat _p_r_hat(_starting2, _R);
    for (int st2 = 0; st2 < _starting2; st2++) {
        for (int r = 0; r < _R; r++) {
            _p_r_hat(st2, r) = 0;
        }
    }
    mat _p_rc_hat(_starting2,_R);
    vector<double> _p_rcc_hat(_R);
    //std::fill(_p_rc_hat.begin(), _p_rc_hat.end(), 0);


    mat _Prob1_hat(_K, _H);
    mat _Prob2_hat(_Q, _R);
    mat _Prob3_hat(_n, _L);
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            _Prob1_hat(k, h) = 0;
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _Prob2_hat(q, r) = 0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3_hat(j, l) = 0;
        }
    }

    Tensor<double, 4> _Prob32(_n, _L, _H, _R);
    Tensor<double, 4> _Prob32_hat(_n, _L, _H, _R);
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32_hat(j, l, h, r) = 0;
                }
            }
        }
    }

    Tensor<double, 4> empiricalpwzy1(_K, _Q, _H, _R);
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    empiricalpwzy1(k, q, h, r) = 0;
                }
            }
        }
    }
    Tensor<double, 4> empiricalpwzy2(_K, _Q, _H, _R);
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    empiricalpwzy2(k, q, h, r) = 0;
                }
            }
        }
    }

    mat _p_h_all(_iter, _H);
    mat _p_r_all(_iter, _R);
    //vector<vector<vector<double>>> _p_h_all;
    mat _logL_all(_starting, _iter);

    mat _Pwy2(_K, _H);
    mat _Pzy2(_Q, _R);
    double explogL2;
    double Entropy2;
    double logLapprox2;

    vector<double> logLtot(_starting2);
    vector<double> explogLtot(_starting2);

    arma_rng::set_seed_random();

    vector<double> _pl(_L);
    for (int l = 0; l < _L; l++) {
        _pl[l] = (double) 1 / _L;
    }

    for (int h = 0; h < _H; h++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pi_l(l, h, r) = (double) 1 / _L;
            }
        }
    }

    vector<double> _de(_n);
    vector<double> _de2(_L);
    mat _tau(_n, _L);

    mat BVR(_I,_I);
    vector<double> BVRG(_I);
    vector<double> BVRG1(_I);
    vector<double> BVRP(_I);
    vector<double> BVRP1(_I);
    vector<double> bvag(_I);
    vector<double> bvag1(_I);
    vector<double> bvap(_I);
    vector<double> bvap1(_I);
    mat bva(_I, _I);
    mat m0(C, C);
    mat m1(_K, C);
    mat m11(_Q, C);
    mat m2(C, C);
    mat m22(C, C);
    mat n(C, C);
    mat n1(_K, C);
    mat n11(_Q, C);
    mat n2(C, C);
    mat n22(C, C);
    arma::cube BVAG(_starting2, _iter, _I);
    arma::cube BVAG1(_starting2, _iter, _I);
    arma::cube BVAP(_starting2, _iter, _I);
    arma::cube BVAP1(_starting2, _iter, _I);
    Tensor<double, 4> BVA(_starting2, _iter, _I, _I);
    mat yy(_n, _I + 3);
    int vp = 0;

    for (int i1 = 0; i1 < _I; i1++) {
        for (int i2 = 0; i2 < _I; i2++) {
            bva(i1, i2) = 0;
        }
    }
    for (int i1 = 0; i1 < _I; i1++) {
        for (int i2 = i1 + 1; i2 < _I; i2++) {

            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    n(c, cc) = 0;
                }
            }
            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    for (int j = 0; j < _n; j++) {
                        if (y(j, i1 + 3) == c && y(j, i2 + 3) == cc) {
                            n(c, cc) += 1;
                        }
                    }
                }
            }


            vector<int> nr(C);
            vector<int> nc(C);
            for (int c = 0; c < C; c++) {
                nr[c] = 0;
            }
            for (int cc = 0; cc < C; cc++) {
                nc[cc] = 0;
            }
            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    nr[c] += n(c, cc);
                }
            }
            for (int cc = 0; cc < C; cc++) {
                for (int c = 0; c < C; c++) {
                    nc[cc] += n(c, cc);
                }
            }
            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    m0(c, cc) = 0;
                }
            }
            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    m0(c, cc) = (double) nr[c] * nc[cc] / _n;
                }
            }


            for (int c = 0; c < _vC[i1]; c++) {
                for (int cc = 0; cc < _vC[i2]; cc++) {
                    bva(i1, i2) += (double) pow((n(c, cc) - m0(c, cc)), 2) / m0(c, cc);
                }
            }
            bva(i1, i2) = (double) bva(i1, i2) / ((_vC[i1] - 1) * (_vC[i2] - 1));

        }
    }
    for (int d = 0; d < nlocdep; d++) {
        bva(locdep(d, 0) - 1, locdep(d, 1) - 1) = 0;
        bva(locdep(d, 1) - 1, locdep(d, 0) - 1) = 0;
    }


    for (int i = 0; i < _I; i++) {
        bvag[i] = 0;
    }

    for (int i1 = 0; i1 < _I; i1++) {

        for (int k = 0; k < _K; k++) {
            for (int cc = 0; cc < C; cc++) {
                n1(k, cc) = 0;
            }
        }
        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < C; cc++) {
                for (int j = 0; j < _n; j++) {
                    if (y(j, 1) - 1 == c && y(j, i1 + 3) == cc) {
                        n1(c, cc) += 1;
                    }
                }
            }
        }

        vector<int> nr1(_K);
        vector<int> nc1(C);
        for (int c = 0; c < _K; c++) {
            nr1[c] = 0;
        }
        for (int cc = 0; cc < C; cc++) {
            nc1[cc] = 0;
        }
        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < C; cc++) {
                nr1[c] += n1(c, cc);
            }
        }

        for (int cc = 0; cc < C; cc++) {
            for (int c = 0; c < _K; c++) {
                nc1[cc] += n1(c, cc);
            }
        }

        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < C; cc++) {
                m1(c, cc) = 0;
            }
        }
        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < C; cc++) {
                m1(c, cc) = (double) nr1[c] * nc1[cc] / _n;
            }
        }

        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                bvag[i1] += (double) pow((n1(c, cc) - m1(c, cc)), 2) / m1(c, cc);
            }
        }

        bvag[i1] = (double) bvag[i1] / ((_K - 1) * (_vC[i1] - 1));

    }


    for (int i = 0; i < _I; i++) {
        bvag1[i] = 0;
    }

    for (int i1 = 0; i1 < _I; i1++) {

        for (int q = 0; q < _Q; q++) {
            for (int cc = 0; cc < C; cc++) {
                n11(q, cc) = 0;
            }
        }
        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < C; cc++) {
                for (int j = 0; j < _n; j++) {
                    if (y(j, 2) - 1 == c && y(j, i1 + 3) == cc) {
                        n11(c, cc) += 1;
                    }
                }
            }
        }

        vector<int> nr11(_Q);
        vector<int> nc11(C);
        for (int c = 0; c < _Q; c++) {
            nr11[c] = 0;
        }
        for (int cc = 0; cc < C; cc++) {
            nc11[cc] = 0;
        }
        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < C; cc++) {
                nr11[c] += n11(c, cc);
            }
        }

        for (int cc = 0; cc < C; cc++) {
            for (int c = 0; c < _Q; c++) {
                nc11[cc] += n11(c, cc);
            }
        }

        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < C; cc++) {
                m11(c, cc) = 0;
            }
        }
        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < C; cc++) {
                m11(c, cc) = (double) nr11[c] * nc11[cc] / _n;
            }
        }

        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                bvag1[i1] += (double) pow((n11(c, cc) - m11(c, cc)), 2) / m11(c, cc);
            }
        }

        bvag1[i1] = (double) bvag1[i1] / ((_Q - 1) * (_vC[i1] - 1));

    }


    for (int i = 0; i < _I; i++) {
        bvap[i] = 0;
    }

    vector<int> vetcoppie(_K);
    for (int k = 0; k < _K; k++) {
        vetcoppie[k] = _nk[k] * (_nk[k] - 1) / 2;
    }
    int totcoppie = 0;
    for (int k = 0; k < _K; k++) {
        totcoppie += vetcoppie[k];
    }

    for (int i1 = 0; i1 < _I; i1++) {

        mat coppie(_K, _n);
        for (int k = 0; k < _K; k++) {
            for (int j = 0; j < _n; j++) {
                coppie(k, j) = -1;
            }
        }
        int co = 0;
        for (int k = 0; k < _K; k++) {
            for (int j = 0; j < _n; j++) {
                if (y(j, 1) - 1 == k) {
                    coppie(k, co) = y(j, i1 + 3);
                    co = co + 1;
                }
            }
            co = 0;
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                n2(c, cc) = 0;
            }
        }

        for (int k = 0; k < _K; k++) {

            for (int cc = 0; cc < _vC[i1]; cc++) {
                for (int c = 0; c < _vC[i1]; c++) {
                    for (int a = 0; a < _nk[k]; a++) {
                        if (coppie(k, a) == cc) {
                            for (int b = (a + 1); b < _nk[k]; b++) {
                                if (coppie(k, b) == c) {
                                    n2(cc, c) += 1;
                                }
                            }
                        }
                    }
                }
            }

        }

        vector<int> nr2(_vC[i1]);
        vector<int> nc2(_vC[i1]);
        for (int c = 0; c < _vC[i1]; c++) {
            nr2[c] = 0;
        }
        for (int cc = 0; cc < _vC[i1]; cc++) {
            nc2[cc] = 0;
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                nr2[c] += n2(c, cc);
            }
        }

        for (int cc = 0; cc < _vC[i1]; cc++) {
            for (int c = 0; c < _vC[i1]; c++) {
                nc2[cc] += n2(c, cc);
            }
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                m2(c, cc) = 0;
            }
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                m2(c, cc) = (double) nr2[c] * nc2[cc] / totcoppie;
            }
        }

        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = c; cc < _vC[i1]; cc++) {
                if (c == cc) {
                    bvap[i1] += (double) pow((n2(c, cc) - m2(c, cc)), 2) / m2(c, cc);
                }
                if (c != cc) {
                    bvap[i1] += (double) pow(((n2(c, cc) + n2(cc, c)) - (m2(c, cc) + m2(cc, c))), 2) /
                                (m2(c, cc) + m2(cc, c));
                }
            }
        }

        bvap[i1] = (double) (bvap[i1] * _K) / (_n * _vC[i1] * (_vC[i1] - 1) / 2);

    }


    for (int i = 0; i < _I; i++) {
        bvap1[i] = 0;
    }

    vector<int> vetcoppie1(_Q);
    for (int q = 0; q < _Q; q++) {
        vetcoppie1[q] = _nq[q] * (_nq[q] - 1) / 2;
    }
    int totcoppie1 = 0;
    for (int q = 0; q < _Q; q++) {
        totcoppie1 += vetcoppie1[q];
    }

    for (int i1 = 0; i1 < _I; i1++) {

        mat coppie1(_Q, _n);
        for (int q = 0; q < _Q; q++) {
            for (int j = 0; j < _n; j++) {
                coppie1(q, j) = -1;
            }
        }
        int co = 0;
        for (int q = 0; q < _Q; q++) {
            for (int j = 0; j < _n; j++) {
                if (y(j, 2) - 1 == q) {
                    coppie1(q, co) = y(j, i1 + 3);
                    co = co + 1;
                }
            }
            co = 0;
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                n22(c, cc) = 0;
            }
        }

        for (int q = 0; q < _Q; q++) {

            for (int cc = 0; cc < _vC[i1]; cc++) {
                for (int c = 0; c < _vC[i1]; c++) {
                    for (int a = 0; a < _nq[q]; a++) {
                        if (coppie1(q, a) == cc) {
                            for (int b = (a + 1); b < _nq[q]; b++) {
                                if (coppie1(q, b) == c) {
                                    n22(cc, c) += 1;
                                }
                            }
                        }
                    }
                }
            }

        }

        vector<int> nr22(_vC[i1]);
        vector<int> nc22(_vC[i1]);
        for (int c = 0; c < _vC[i1]; c++) {
            nr22[c] = 0;
        }
        for (int cc = 0; cc < _vC[i1]; cc++) {
            nc22[cc] = 0;
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                nr22[c] += n22(c, cc);
            }
        }

        for (int cc = 0; cc < _vC[i1]; cc++) {
            for (int c = 0; c < _vC[i1]; c++) {
                nc22[cc] += n22(c, cc);
            }
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                m22(c, cc) = 0;
            }
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                m22(c, cc) = (double) nr22[c] * nc22[cc] / totcoppie1;
            }
        }

        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = c; cc < _vC[i1]; cc++) {
                if (c == cc) {
                    bvap1[i1] += (double) pow((n22(c, cc) - m22(c, cc)), 2) / m22(c, cc);
                }
                if (c != cc) {
                    bvap1[i1] += (double) pow(((n22(c, cc) + n22(cc, c)) - (m22(c, cc) + m22(cc, c))), 2) /
                                 (m22(c, cc) + m22(cc, c));
                }
            }
        }

        bvap1[i1] = (double) (bvap1[i1] * _Q) / (_n * _vC[i1] * (_vC[i1] - 1) / 2);

    }


    mat y2(_n, _G - nlocdep);

    if (nlocdep > 0) {

        vector<int> ynew(_n);
        int vn;

        for (int j = 0; j < _n; j++) {
            y2(j, 0) = y(j, 0);
        }
        for (int j = 0; j < _n; j++) {
            y2(j, 1) = y(j, 1);
        }
        for (int j = 0; j < _n; j++) {
            y2(j, 2) = y(j, 2);
        }
        for (int d = 0; d < nlocdep; d++) {
            vn = 0;
            for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                    for (int nn = 0; nn < _n; nn++) {
                        if (y(nn, locdep(d, 0) + 2) == c && y(nn, locdep(d, 1) + 2) == cc) {
                            ynew[nn] = vn;
                        }
                    }
                    vn = vn + 1;
                }
            }
            for (int j = 0; j < _n; j++) {
                y2(j, d + 3) = ynew[j];
            }
        }

        vector<int> vd(_I - 2 * nlocdep);
        int tg = 0;
        int dd = 0;
        if ((_I - 2 * nlocdep) > 1) {
            for (int i = 1; i < (_I + 1); i++) {
                for (int d = 0; d < nlocdep; d++) {
                    if (locdep(d, 0) != i && locdep(d, 1) != i) {
                        dd = dd + 1;
                    }
                }
                if (dd == nlocdep) {
                    vd[tg] = i;
                    tg = tg + 1;
                }
                dd = 0;
            }

            _I = _G - 3 - nlocdep;
            for (int f = nlocdep; f < _I; f++) {
                for (int j = 0; j < _n; j++) {
                    y2(j, f + 3) = y(j, vd[f - nlocdep] + 2);
                }
            }

        }

        int _vd;
        if ((_I - 2 * nlocdep) == 1) {
            dd = 0;
            for (int i = 1; i < (_I + 1); i++) {
                for (int d = 0; d < nlocdep; d++) {
                    if (locdep(d, 0) != i && locdep(d, 1) != i) {
                        dd = dd + 1;
                    }
                }
                if (dd == nlocdep) {
                    _vd = i;
                }
                dd = 0;
            }

            _I = _G - 3 - nlocdep;
            for (int f = nlocdep; f < _I; f++) {
                for (int j = 0; j < _n; j++) {
                    y2(j, f + 3) = y(j, _vd + 2);
                }
            }

        }

    }


    if (nlocdep == 0) {
        for (int j = 0; j < _n; j++) {
            for (int i = 0; i < _G; i++) {
                y2(j, i) = y(j, i);
            }
        }
    }

    _I = _G - 3 - nlocdep;


    vector<int> vC(_I);
    //int ma;
    for (int i = 0; i < _I; i++) {
        ma = y2(0, i + 3);
        for (int j = 0; j < _n; j++) {
            if (y2(j, i + 3) > ma) {
                ma = y2(j, i + 3);
            }
        }
        vC[i] = ma + 1;
    }

    int _C;
    //int maa;
    maa = vC[0];
    for (int i = 0; i < _I; i++) {
        if (vC[i] > maa) {
            maa = vC[i];
        }
    }
    _C = maa;


    Tensor<double, 5> _num5(_I, _C, _L, _H, _R);
    //Tensor<double,5> _aaa(_I,_n,_L,_H,_R);
    Tensor<double, 4> _P(_starting, _I, _C, _L);
    Tensor<double, 6> _p_hat(_starting2, _I, _C, _L, _H, _R);
    for (int st2 = 0; st2 < _starting2; st2++) {
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _p_hat(st2, i, c, l, h, r) = 0;
                        }
                    }
                }
            }
        }
    }
    arma::cube _p1(_I, _C, _L, fill::randu);
    Tensor<double, 5> _p(_I, _C, _L, _H, _R);
    //_p.setRandom();
    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _p(i, c, l, h, r) = 0;
                    }
                }
            }
        }
    }


    for (int st = 0; st < _starting; st++) {

        arma::cube _pb(_I, _C, _L);
        double eqprob;
        vector<double> vec(_C * _L * _I);
        eqprob = (double) 1 / (_C * _L * _I);
        std::fill(vec.begin(), vec.end(), eqprob);
        std::random_device rdtest;
        std::mt19937 generator2(rdtest());
        boost::random::discrete_distribution<int> distribution3(vec.begin(), vec.end());
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    int sample = distribution3(generator2);
                    _pb(i, c, l) = sample + 1;
                }
            }
        }

        mat sumc(_I, _L);
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                sumc(i, l) = 0;
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int c = 0; c < vC[i]; c++) {
                    sumc(i, l) += _pb(i, c, l);
                }
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int c = 0; c < vC[i]; c++) {
                    _p1(i, c, l) = (double) _pb(i, c, l) / sumc(i, l);
                }
            }
        }

        vector<int> _z0(_L);
        double eqprob0;
        vector<double> vec0(_L);
        eqprob0 = (double) 1 / (_L);
        std::fill(vec0.begin(), vec0.end(), eqprob0);
        boost::random::discrete_distribution<int> distribution0(vec0.begin(), vec0.end());
        for (int l = 0; l < _L; l++) {
            int sample = distribution0(generator2);
            _z0[l] = sample + 1;
        }

        int sumc0 = 0;
        for (int l = 0; l < _L; l++) {
            sumc0 += _z0[l];
        }

        for (int l = 0; l < _L; l++) {
            _pl[l] = (double) _z0[l] / sumc0;
        }

        mat _f1(_n, _L);
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _f1(j, l) = 1;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int i = 0; i < _I; i++) {
                    _f1(j, l) = _f1(j, l) * _p1(i, y2(j, i + 3), l);
                }
            }
        }

        double L, L1;
        vector<double> _flog(_n);
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _flog[j] += _pl[l] * _f1(j, l);
            }
        }
        for (int j = 0; j < _n; j++) {
            _flog[j] = log(_flog[j]);
        }
        L1 = 0;
        for (int j = 0; j < _n; j++) {
            L1 += _flog[j];
        }

        for (int u = 0; u < _iter0; u++) {

            mat num0(_n, _L);
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = (double) log(_f1(j, l)) + log(_pl[l]);
                }
            }
            vector<double> _max0(_n);
            for (int j = 0; j < _n; j++) {
                _max0[j] = num0(j, 0);
                for (int l = 0; l < _L; l++) {
                    if (num0(j, l) > _max0[j]) {
                        _max0[j] = num0(j, l);
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = num0(j, l) - _max0[j];
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = exp(num0(j, l));
                }
            }

            vector<double> _den0(_n);
            for (int j = 0; j < _n; j++) {
                _den0[j] = 0;
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _den0[j] += num0(j, l);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _tau(j, l) = (double) num0(j, l) / _den0[j];
                }
            }

            for (int l = 0; l < _L; l++) {
                _pl[l] = 0;
            }

            for (int l = 0; l < _L; l++) {
                for (int j = 0; j < _n; j++) {
                    _pl[l] += _tau(j, l);
                }
            }

            for (int l = 0; l < _L; l++) {
                _pl[l] = _pl[l] / _n;
            }

            for (int l = 0; l < _L; l++) {
                _de2[l] = 0;
            }

            for (int l = 0; l < _L; l++) {
                for (int j = 0; j < _n; j++) {
                    _de2[l] += _tau(j, l);
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        _p1(i, c, l) = 0;
                    }
                }
            }

            if (nlocdep == 0) {

                for (int i = 0; i < _I; i++) {
                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i, y2(j, i + 3), l) += _tau(j, l);
                        }
                    }
                }

                for (int i = 0; i < _I; i++) {
                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            _p1(i, c, l) = (double) _p1(i, c, l) / _de2[l];
                        }
                    }
                }

            }

            if (nlocdep > 0) {

                for (int i = nlocdep; i < _I; i++) {
                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i, y2(j, i + 3), l) += _tau(j, l);
                        }
                    }
                }

                for (int i = nlocdep; i < _I; i++) {
                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            _p1(i, c, l) = (double) _p1(i, c, l) / _de2[l];
                        }
                    }
                }

                for (int d = 0; d < nlocdep; d++) {

                    mat nyy(_vC[locdep(d, 0) - 1], _vC[locdep(d, 1) - 1]);
                    for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                        for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                            nyy(c, cc) = 0;
                        }
                    }

                    for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                        for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                            for (int j = 0; j < _n; j++) {
                                if (y(j, locdep(d, 0) + 2) == c && y(j, locdep(d, 1) + 2) == cc) {
                                    nyy(c, cc) += 1;
                                }
                            }
                        }
                    }

                    mat nyx1(_vC[locdep(d, 0) - 1], _L);
                    for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                        for (int l = 0; l < _L; l++) {
                            nyx1(c, l) = 0;
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            nyx1(y(j, locdep(d, 0) + 2), l) += _tau(j, l);
                        }
                    }


                    mat nyx2(_vC[locdep(d, 1) - 1], _L);
                    for (int c = 0; c < _vC[locdep(d, 1) - 1]; c++) {
                        for (int l = 0; l < _L; l++) {
                            nyx2(c, l) = 0;
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            nyx2(y(j, locdep(d, 1) + 2), l) += _tau(j, l);
                        }
                    }

                    arma::cube myx(_vC[locdep(d, 0) - 1], _vC[locdep(d, 1) - 1], _L);
                    arma::cube myxp(_vC[locdep(d, 0) - 1], _vC[locdep(d, 1) - 1], _L);
                    arma::cube myxs(_vC[locdep(d, 0) - 1], _vC[locdep(d, 1) - 1], _L);
                    mat myx1(_vC[locdep(d, 0) - 1], _L);
                    mat myx2(_vC[locdep(d, 1) - 1], _L);
                    mat my(_vC[locdep(d, 0) - 1], _vC[locdep(d, 1) - 1]);
                    for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                        for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                            for (int l = 0; l < _L; l++) {
                                myx(c, cc, l) = 1;
                            }
                        }
                    }

                    //ipf
                    for (int uu = 0; uu < _iter2; uu++) {

                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int l = 0; l < _L; l++) {
                                myx1(c, l) = 0;
                            }
                        }
                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int l = 0; l < _L; l++) {
                                for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                    myx1(c, l) += myx(c, cc, l);
                                }
                            }
                        }

                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                for (int l = 0; l < _L; l++) {
                                    if (myx1(c, l) != 0) {
                                        myxp(c, cc, l) = (double) myx(c, cc, l) * nyx1(c, l) / myx1(c, l);
                                    }
                                    if (myx1(c, l) == 0) {
                                        myxp(c, cc, l) = 0;
                                    }
                                }
                            }
                        }

                        for (int c = 0; c < _vC[locdep(d, 1) - 1]; c++) {
                            for (int l = 0; l < _L; l++) {
                                myx2(c, l) = 0;
                            }
                        }

                        for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                            for (int l = 0; l < _L; l++) {
                                for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                                    myx2(cc, l) += myxp(c, cc, l);
                                }
                            }
                        }

                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                for (int l = 0; l < _L; l++) {
                                    if (myx2(cc, l) != 0) {
                                        myxs(c, cc, l) = (double) myxp(c, cc, l) * nyx2(cc, l) / myx2(cc, l);
                                    }
                                    if (myx2(cc, l) == 0) {
                                        myxs(c, cc, l) = 0;
                                    }
                                }
                            }
                        }

                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                my(c, cc) = 0;
                            }
                        }
                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                for (int l = 0; l < _L; l++) {
                                    my(c, cc) += myxs(c, cc, l);
                                }
                            }
                        }

                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                for (int l = 0; l < _L; l++) {
                                    if (my(c, cc) != 0) {
                                        myx(c, cc, l) = (double) myxs(c, cc, l) * nyy(c, cc) / my(c, cc);
                                    }
                                    if (my(c, cc) == 0) {
                                        myx(c, cc, l) = 0;
                                    }
                                }
                            }
                        }

                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                for (int l = 0; l < _L; l++) {
                                    if (abs(myx(c, cc, l) - myxp(c, cc, l)) < eps4 &&
                                        abs(myxp(c, cc, l) - myxs(c, cc, l)) < eps4) {
                                        uu = _iter2 - 1;
                                    }
                                }
                            }
                        }


                    }

                    vector<double> nx(_L);
                    for (int l = 0; l < _L; l++) {
                        nx[l] = 0;
                    }
                    for (int l = 0; l < _L; l++) {
                        for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                            for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                                nx[l] += myx(c, cc, l);
                            }
                        }
                    }

                    int ar = 0;
                    for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                        for (int cc = 0; cc < _vC[locdep(d, 1) - 1]; cc++) {
                            for (int l = 0; l < _L; l++) {
                                _p1(d, ar, l) = (double) myx(c, cc, l) / nx[l];
                            }
                            ar = ar + 1;
                        }
                    }

                }

            }

            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int c = 0; c < _C; c++) {
                        if (_p1(i, c, l) > 0.9999) {
                            _p1(i, c, l) = 0.999;
                        }
                    }
                }
            }
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int c = 0; c < _C; c++) {
                        if (_p1(i, c, l) < 0.00001) {
                            _p1(i, c, l) = 0.0001;
                        }
                    }
                }
            }


            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _f1(j, l) = 1;
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int i = 0; i < _I; i++) {
                        _f1(j, l) = _f1(j, l) * _p1(i, y2(j, i + 3), l);
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                _flog[j] = 0;
            }
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _flog[j] += _pl[l] * _f1(j, l);
                }
            }
            for (int j = 0; j < _n; j++) {
                _flog[j] = log(_flog[j]);
            }
            L = 0;
            for (int j = 0; j < _n; j++) {
                L += _flog[j];
            }
            if (abs(L1 - L) < eps0) {
                u = _iter0 - 1;
            }
            L1 = L;

            //cout<<u<<endl;
            //cout<<L1<<endl;

        }


        vector<double> tempp(_C);
        vector<double> tempp1(_n);
        for (int l = 0; l < _L - 1; l++) {
            min1 = l;
            for (int j = l + 1; j < _L; j++) {
                if (_pl[j] < _pl[min1]) {
                    min1 = j;
                }
            }
            temp1 = _pl[l];
            _pl[l] = _pl[min1];
            _pl[min1] = temp1;
            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    tempp[c] = _p1(i, c, l);
                    _p1(i, c, l) = _p1(i, c, min1);
                    _p1(i, c, min1) = tempp[c];
                }
            }
            for (int j = 0; j < _n; j++) {
                tempp1[j] = _tau(j, l);
                _tau(j, l) = _tau(j, min1);
                _tau(j, min1) = tempp1[j];
            }
            for (int j = 0; j < _n; j++) {
                tempp1[j] = _f1(j, l);
                _f1(j, l) = _f1(j, min1);
                _f1(j, min1) = tempp1[j];
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _Tau(st, j, l) = (double) _tau(j, l);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _F(st, j, l) = (double) _f1(j, l);
            }
        }

        //Seed[st] = rdtest();
        logLik[st] = L;
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    _P(st, i, c, l) = _p1(i, c, l);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            _PL(st, l) = (double) _pl[l];
        }

    }

    double _maxlogl;
    int indicestarting;
    _maxlogl = logLik[0];
    for (int st = 0; st < _starting; st++) {
        if (logLik[st] > _maxlogl) {
            _maxlogl = logLik[st];
            indicestarting = st;
        }
    }

    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _p(i, c, l, h, r) = _P(indicestarting, i, c, l);
                    }
                }
            }
        }
    }

    for (int st2 = 0; st2 < _starting2; st2++) {

        mu = 0;
        vp = 0;

        explogL2 = 0;
        Entropy2 = 0;
        logLapprox2 = 0;

        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _pi_l_hat(st2, l, h, r) = 0;
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _pi_lc_hat(l, h, r) = 0;
                }
            }
        }

        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _p_hat(st2, i, c, l, h, r) = 0;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < (_I + nlocdep); i++) {
            for (int c = 0; c < C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _p_finale(i, c, l, h, r) = 0;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < (_I + nlocdep); i++) {
            for (int c = 0; c < C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _p_finale2(st2, i, c, l, h, r) = 0;
                        }
                    }
                }
            }
        }

        for (int h = 0; h < _H; h++) {
            _p_h_hat(st2, h) = 0;
        }

        for (int r = 0; r < _R; r++) {
            _p_r_hat(st2, r) = 0;
        }

        for (int h = 0; h < _H; h++) {
            _p_hc_hat(st2, h) = 0;
        }

        for (int r = 0; r < _R; r++) {
            _p_rc_hat(st2, r) = 0;
        }
        std::fill(_p_hcc_hat.begin(), _p_hcc_hat.end(), 0);
        std::fill(_p_rcc_hat.begin(), _p_rcc_hat.end(), 0);


        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _Prob1_hat(k, h) = 0;
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _Prob2_hat(q, r) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _Prob3_hat(j, l) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _Prob32_hat(j, l, h, r) = 0;
                    }
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        empiricalpwzy1(k, q, h, r) = 0;
                    }
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        empiricalpwzy2(k, q, h, r) = 0;
                    }
                }
            }
        }

        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _p(i, c, l, h, r) = _P(indicestarting, i, c, l);
                        }
                    }
                }
            }
        }

        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _pi_l(l, h, r) = (double) 1 / _L;
                }
            }
        }


        vector<int> _cv(_H);
        int _pcv = 0;
        vector<int> _cv2(_R);
        int _pcv2 = 0;

        arma::vec dir_param_pr(_R);
        for (int r = 0; r < _R; r++) {
            dir_param_pr(r) = 4;
        }
        arma::vec dir_param_ph(_H);
        for (int h = 0; h < _H; h++) {
            dir_param_ph(h) = 4;
        }

        arma::cube dir_param_plhr(_L,_H,_R);
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    //dir_param_plhr(l,h, r) = _pi_l(l,h, r);
                    dir_param_plhr(l,h, r) = 4;
                }
            }
        }

        arma::cube dir_param_p(_I, _C, _L);
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            dir_param_p(i, c, l) =1;// _p(i, c, l, h, r);
                        }
                    }
                }
            }
        }


        for (int u1 = 0; u1 < _iter; u1++) {

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _f(j, l, h, r) = 1;
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int i = 0; i < _I; i++) {
                                _f(j, l, h, r) = _f(j, l, h, r) * _p(i, y2(j, i + 3), l, h, r);
                            }
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int l = 0; l < _L; l++) {
                            _pf(j, l, h, r) = (double) _pi_l(l, h, r) * _f(j, l, h, r);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _fykq(j, h, r) = 0;
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int l = 0; l < _L; l++) {
                            _fykq(j, h, r) += _pf(j, l, h, r);
                        }
                    }
                }
            }

            //for(int gi=0;gi<_iter2;gi++){

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    _ff(k, h) = 0;
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int q = 0; q < _Q; q++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _res(k, q, h, r) = 0;
                        }
                    }
                }
            }

            int m = 0;
            for (int k = 0; k < _K; k++) {
                for (int q = 0; q < _Q; q++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int j = m; j < _nkq(k, q) + m; j++) {
                                _res(k, q, h, r) += log(_fykq(j, h, r));
                            }
                        }
                    }
                    m = m + _nkq(k, q);
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int q = 0; q < _Q; q++) {
                        for (int r = 0; r < _R; r++) {
                            _res(k, q, h, r) = _z(q, r) * _res(k, q, h, r);
                        }
                    }
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int q = 0; q < _Q; q++) {
                            _ff(k, h) += _res(k, q, h, r);
                        }
                    }
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    _num2(k, h) = log(_p_h[h]) + _ff(k, h);
                }
            }

            vector<double> _maxh(_K);
            for (int k = 0; k < _K; k++) {
                _maxh[k] = _num2(k, 0);
                for (int h = 0; h < _H; h++) {
                    if (_num2(k, h) > _maxh[k]) {
                        _maxh[k] = _num2(k, h);
                    }
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    _num2(k, h) = _num2(k, h) - _maxh[k];
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    _num2(k, h) = exp(_num2(k, h));
                }
            }

            for (int k = 0; k < _K; k++) {
                _den2[k] = 0;
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    _den2[k] += _num2(k, h);
                }
            }


            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    _Prob1(k, h) = (double) _num2(k, h) / _den2[k];
                }
            }

            _pcv = 0;
            while (_pcv == 0) {
                _pcv = 1;
                for (int h = 0; h < _H; h++) {
                    _cv[h] = 0;
                }
                _w.zeros(_K, _H);
                //boost::mt19937 generator2(1234);
                for (int k = 0; k < _K; k++) {
                    boost::random::discrete_distribution<int> distribution4(_Prob1.row(k).begin(), _Prob1.row(k).end());
                    int sample = distribution4(generator);
                    _w(k, sample) = 1;
                }

                for (int h = 0; h < _H; h++) {
                    for (int k = 0; k < _K; k++) {
                        _cv[h] += _w(k, h);
                    }
                }
                for (int h = 0; h < _H; h++) {
                    _pcv = _pcv * _cv[h];
                }
            }

            std::fill(_numh.begin(), _numh.end(), 0);

            for (int h = 0; h < _H; h++) {
                for (int k = 0; k < _K; k++) {
                    _numh[h] += _w(k, h);
                }
            }

            vector<int> _nh(_H);
            for (int h = 0; h < _H; h++) {
                _nh[h]=0;
            }
            for (int h = 0; h < _H; h++) {
                for (int k = 0; k < _K; k++) {
                    _nh[h] += _w(k,h);
                }
            }

            for (int h = 0; h < _H; h++) {
                dir_param_ph(h) += _nh[h];
            }
            arma::vec sample_dir_ph = rdirichlet(dir_param_ph);
            for (int h = 0; h < _H; h++) {
                _p_h[h] = sample_dir_ph(h);
            }


            vector<double> _p_hc(_H);
            for (int h = 0; h < _H; h++) {
                _p_hc[h] = _p_h[h];
            }

            mat temp5(_L, _R);
            for (int i = 0; i < (_H - 1); i++) {
                min1 = i;
                for (int j = (i + 1); j < _H; j++) {
                    if (_p_h[j] < _p_h[min1]) {
                        min1 = j;
                    }
                }
                temp1 = _p_h[i];
                _p_h[i] = _p_h[min1];
                _p_h[min1] = temp1;

            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _ff2(q, r) = 0;
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _res2(q, k, h, r) = 0;
                        }
                    }
                }
            }

            int mm = 0;
            for (int q = 0; q < _Q; q++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int j = mm; j < _nkq(k, q) + mm; j++) {
                                _res2(q, k, h, r) += log(_fykq(_indvecqk[j], h, r));
                            }
                        }
                    }
                    mm = mm + _nkq(k, q);
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    for (int k = 0; k < _K; k++) {
                        for (int h = 0; h < _H; h++) {
                            _res2(q, k, h, r) = _w(k, h) * _res2(q, k, h, r);
                        }
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    for (int k = 0; k < _K; k++) {
                        for (int h = 0; h < _H; h++) {
                            _ff2(q, r) += _res2(q, k, h, r);
                        }
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _num3(q, r) = log(_p_r[r]) + _ff2(q, r);
                }
            }

            vector<double> _maxr(_Q);
            for (int q = 0; q < _Q; q++) {
                _maxr[q] = _num3(q, 0);
                for (int r = 0; r < _R; r++) {
                    if (_num3(q, r) > _maxr[q]) {
                        _maxr[q] = _num3(q, r);
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _num3(q, r) = _num3(q, r) - _maxr[q];
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _num3(q, r) = exp(_num3(q, r));
                }
            }

            for (int q = 0; q < _Q; q++) {
                _den3[q] = 0;
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _den3[q] += _num3(q, r);
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _Prob2(q, r) = (double) _num3(q, r) / _den3[q];
                }
            }

            _pcv2 = 0;
            while (_pcv2 == 0) {
                _pcv2 = 1;
                for (int r = 0; r < _R; r++) {
                    _cv2[r] = 0;
                }
                _z.zeros(_Q, _R);
                for (int q = 0; q < _Q; q++) {
                    boost::random::discrete_distribution<int> distribution5(_Prob2.row(q).begin(), _Prob2.row(q).end());
                    int sample = distribution5(generator);
                    _z(q, sample) = 1;
                }

                for (int r = 0; r < _R; r++) {
                    for (int q = 0; q < _Q; q++) {
                        _cv2[r] += _z(q, r);
                    }
                }
                for (int r = 0; r < _R; r++) {
                    _pcv2 = _pcv2 * _cv2[r];
                }
            }

            std::fill(_numr.begin(), _numr.end(), 0);

            for (int r = 0; r < _R; r++) {
                for (int q = 0; q < _Q; q++) {
                    _numr[r] += _z(q, r);
                }
            }

            vector<int> _nr(_R);
            for (int r = 0; r < _R; r++) {
                _nr[r]=0;
            }
            for (int r = 0; r < _R; r++) {
                for (int q = 0; q < _Q; q++) {
                    _nr[r] += _z(q,r);
                }
            }

            for (int r = 0; r < _R; r++) {
                dir_param_pr(r) += _nr[r];
            }
            arma::vec sample_dir_pr = rdirichlet(dir_param_pr);
            for (int r = 0; r < _R; r++) {
                _p_r[r] = sample_dir_pr(r);
            }


            vector<double> _p_rc(_R);
            for (int r = 0; r < _R; r++) {
                _p_rc[r] = _p_r[r];
            }

            double temp6;
            vector<int> temp7(_Q);
            mat temp8(_n, _H);
            arma::cube temp9(_n, _L, _H);
            mat temp10(_L, _H);
            for (int i = 0; i < (_R - 1); i++) {
                min2 = i;
                for (int j = (i + 1); j < _R; j++) {
                    if (_p_r[j] < _p_r[min2]) {
                        min2 = j;
                    }
                }
                temp6 = _p_r[i];
                _p_r[i] = _p_r[min2];
                _p_r[min2] = temp6;

            }

            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    _wjk(j, h) = 0;
                }
            }
            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    _wjk(j, h) = _w(y2(j, 1) - 1, h);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int r = 0; r < _R; r++) {
                    _zjq(j, r) = 0;
                }
            }
            for (int j = 0; j < _n; j++) {
                for (int r = 0; r < _R; r++) {
                    _zjq(j, r) = _z(y2(j, 2) - 1, r);
                }
            }


            for (int k = 0; k < _K; k++) {
                for (int q = 0; q < _Q; q++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            empiricalpwzy2(k, q, h, r) = 0;
                        }
                    }
                }
            }

            for (int gi = 0; gi < _iter4; gi++) {

                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _ff(k, h) = 0;
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int q = 0; q < _Q; q++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _res(k, q, h, r) = 0;
                            }
                        }
                    }
                }

                m = 0;
                for (int k = 0; k < _K; k++) {
                    for (int q = 0; q < _Q; q++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                for (int j = m; j < _nkq(k, q) + m; j++) {
                                    _res(k, q, h, r) += log(_fykq(j, h, r));
                                }
                            }
                        }
                        m = m + _nkq(k, q);
                    }
                }


                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        for (int q = 0; q < _Q; q++) {
                            for (int r = 0; r < _R; r++) {
                                _res(k, q, h, r) = _z(q, r) * _res(k, q, h, r);
                            }
                        }
                    }
                }


                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int q = 0; q < _Q; q++) {
                                _ff(k, h) += _res(k, q, h, r);
                            }
                        }
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _num2(k, h) = log(_p_hc[h]) + _ff(k, h);
                    }
                }

                //vector<double> _maxh(_K);
                for (int k = 0; k < _K; k++) {
                    _maxh[k] = _num2(k, 0);
                    for (int h = 0; h < _H; h++) {
                        if (_num2(k, h) > _maxh[k]) {
                            _maxh[k] = _num2(k, h);
                        }
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _num2(k, h) = _num2(k, h) - _maxh[k];
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _num2(k, h) = exp(_num2(k, h));
                    }
                }

                for (int k = 0; k < _K; k++) {
                    _den2[k] = 0;
                }

                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _den2[k] += _num2(k, h);
                    }
                }


                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _Prob1(k, h) = (double) _num2(k, h) / _den2[k];
                    }
                }


                _pcv = 0;
                while (_pcv == 0) {
                    _pcv = 1;
                    for (int h = 0; h < _H; h++) {
                        _cv[h] = 0;
                    }
                    _w.zeros(_K, _H);
//boost::mt19937 generator2(1234);
                    for (int k = 0; k < _K; k++) {
                        boost::random::discrete_distribution<int> distribution7(_Prob1.row(k).begin(),
                                                                                _Prob1.row(k).end());
                        int sample = distribution7(generator);
                        _w(k, sample) = 1;
                    }
                    for (int h = 0; h < _H; h++) {
                        for (int k = 0; k < _K; k++) {
                            _cv[h] += _w(k, h);
                        }
                    }
                    for (int h = 0; h < _H; h++) {
                        _pcv = _pcv * _cv[h];
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _ff2(q, r) = 0;
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int k = 0; k < _K; k++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _res2(q, k, h, r) = 0;
                            }
                        }
                    }
                }

                mm = 0;
                for (int q = 0; q < _Q; q++) {
                    for (int k = 0; k < _K; k++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                for (int j = mm; j < _nkq(k, q) + mm; j++) {
                                    _res2(q, k, h, r) += log(_fykq(_indvecqk[j], h, r));
                                }
                            }
                        }
                        mm = mm + _nkq(k, q);
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        for (int k = 0; k < _K; k++) {
                            for (int h = 0; h < _H; h++) {
                                _res2(q, k, h, r) = _w(k, h) * _res2(q, k, h, r);
                            }
                        }
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        for (int k = 0; k < _K; k++) {
                            for (int h = 0; h < _H; h++) {
                                _ff2(q, r) += _res2(q, k, h, r);
                            }
                        }
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _num3(q, r) = log(_p_rc[r]) + _ff2(q, r);
                    }
                }

                //vector<double> _maxr(_Q);
                for (int q = 0; q < _Q; q++) {
                    _maxr[q] = _num3(q, 0);
                    for (int r = 0; r < _R; r++) {
                        if (_num3(q, r) > _maxr[q]) {
                            _maxr[q] = _num3(q, r);
                        }
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _num3(q, r) = _num3(q, r) - _maxr[q];
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _num3(q, r) = exp(_num3(q, r));
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    _den3[q] = 0;
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _den3[q] += _num3(q, r);
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _Prob2(q, r) = (double) _num3(q, r) / _den3[q];
                    }
                }

                _pcv2 = 0;
                while (_pcv2 == 0) {
                    _pcv2 = 1;
                    for (int r = 0; r < _R; r++) {
                        _cv2[r] = 0;
                    }
                    _z.zeros(_Q, _R);
                    for (int q = 0; q < _Q; q++) {
                        boost::random::discrete_distribution<int> distribution8(_Prob2.row(q).begin(),
                                                                                _Prob2.row(q).end());
                        int sample = distribution8(generator);
                        _z(q, sample) = 1;
                    }

                    for (int r = 0; r < _R; r++) {
                        for (int q = 0; q < _Q; q++) {
                            _cv2[r] += _z(q, r);
                        }
                    }
                    for (int r = 0; r < _R; r++) {
                        _pcv2 = _pcv2 * _cv2[r];
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int q = 0; q < _Q; q++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                empiricalpwzy2(k, q, h, r) += _w(k, h) * _z(q, r);
                            }
                        }
                    }
                }

            }

            for (int k = 0; k < _K; k++) {
                for (int q = 0; q < _Q; q++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            empiricalpwzy2(k, q, h, r) = (double) empiricalpwzy2(k, q, h, r) / _iter4;
                        }
                    }
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    _Pwy2(k, h) = 0;
                }
            }
            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _Pzy2(q, r) = 0;
                }
            }
            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _Pwy2(k, h) += empiricalpwzy2(k, 0, h, r);
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    for (int h = 0; h < _H; h++) {
                        _Pzy2(q, r) += empiricalpwzy2(0, q, h, r);
                    }
                }
            }


            //}

            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _B(j, h, r) = 0;
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _B(j, h, r) = _wjk(j, h) * _zjq(j, r);
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Prob32(j, l, h, r) = log(_pi_l(l, h, r)) + log(_f(j, l, h, r));
                        }
                    }
                }
            }

            Tensor<double, 3> _maxl2(_n, _H, _R);
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int j = 0; j < _n; j++) {
                        _maxl2(j, h, r) = _Prob32(j, 0, h, r);
                        for (int l = 0; l < _L; l++) {
                            if (_Prob32(j, l, h, r) > _maxl2(j, h, r)) {
                                _maxl2(j, h, r) = _Prob32(j, l, h, r);
                            }
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Prob32(j, l, h, r) = _Prob32(j, l, h, r) - _maxl2(j, h, r);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Prob32(j, l, h, r) = exp(_Prob32(j, l, h, r));
                        }
                    }
                }
            }

            Tensor<double, 3> _denomi(_n, _H, _R);
            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _denomi(j, h, r) = 0;
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int l = 0; l < _L; l++) {
                            _denomi(j, h, r) += (double) _Prob32(j, l, h, r);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Prob32(j, l, h, r) = (double) _Prob32(j, l, h, r) / _denomi(j, h, r);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            if (_B(j, h, r) == 1) {
                                _Prob3(j, l) = log(_pi_l(l, h, r)) + log(_f(j, l, h, r));
                            }
                        }
                    }
                }
            }


            vector<double> _maxl(_n);
            for (int j = 0; j < _n; j++) {
                _maxl[j] = _Prob3(j, 0);
                for (int l = 0; l < _L; l++) {
                    if (_Prob3(j, l) > _maxl[j]) {
                        _maxl[j] = _Prob3(j, l);
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob3(j, l) = _Prob3(j, l) - _maxl[j];
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob3(j, l) = exp(_Prob3(j, l));
                }
            }

            for (int j = 0; j < _n; j++) {
                _denom[j] = 0;
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _denom[j] += _Prob3(j, l);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob3(j, l) = (double) _Prob3(j, l) / _denom[j];
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _xjkq(j, l) = 0;
                }
            }

            for (int j = 0; j < _n; j++) {
                boost::random::discrete_distribution<int> distribution6(_Prob3.row(j).begin(), _Prob3.row(j).end());
                int sample = distribution6(generator);
                _xjkq(j, sample) = 1;
            }

            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _den4(h, r) = 0;
                }
            }
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int j = 0; j < _n; j++) {
                        _den4(h, r) += _B(j, h, r);
                    }
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _num4(l, h, r) = 0;
                    }
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = 0; j < _n; j++) {
                            _num4(l, h, r) += _B(j, h, r) * _xjkq(j, l);
                        }
                    }
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        dir_param_plhr(l, h, r) += _num4(l, h, r);
                    }
                }
            }

            arma::vec dir_param_pl(_L);
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int l1 = 0; l1 < _L; l1++) {
                        dir_param_pl(l1) = dir_param_plhr(l1, h, r);
                    }
                    arma::vec sample_dir_pl = rdirichlet(dir_param_pl);
                    for (int l = 0; l < _L; l++) {
                        _pi_l(l, h, r) = sample_dir_pl(l);
                    }
                }
            }

            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int l = 0; l < _L; l++) {
                        if (_pi_l(l, h, r) == 0) {
                            _pi_l(l, h, r) = 0.0001;
                        }
                        if (_pi_l(l, h, r) == 1) {
                            _pi_l(l, h, r) = 0.999;
                        }
                    }
                }
            }

            //if (u1 >= _burn && u1 != it[mu]) {

            double explogL = 0;
            Tensor<double, 4> _Pxwzy(_n, _L, _H, _R);
            arma::cube _Pwzyest(_n, _H, _R);

            for (int j = 0; j < _n; j++) {
                for (int k = 0; k < _K; k++) {
                    for (int q = 0; q < _Q; q++) {
                        if (y2(j, 1) == (k + 1) && y2(j, 2) == (q + 1)) {
                            for (int h = 0; h < _H; h++) {
                                for (int r = 0; r < _R; r++) {
                                    _Pwzyest(j, h, r) = empiricalpwzy2(k, q, h, r);
                                }
                            }
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            //_Pxwzy(j, l, h, r) = (double) _Pwzyest(j, h, r) * _Prob3(j, l);
                            _Pxwzy(j, l, h, r) =(double) _Pwzyest(j, h, r) * _Prob32(j, l,h,r);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob4l(j, l) = 0;
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Prob4l(j, l) += (double) _Pxwzy(j, l, h, r);
                        }
                    }
                }
            }

            mat _Prob1est(_n, _H);
            for (int j = 0; j < _n; j++) {
                for (int k = 0; k < _K; k++) {
                    if (y(j, 1) == (k + 1)) {
                        for (int h = 0; h < _H; h++) {
                            _Prob1est(j, h) =(double) _Prob1(k, h);
                        }
                    }
                }
            }
            mat _Prob2est(_n, _R);
            for (int j = 0; j < _n; j++) {
                for (int q = 0; q < _Q; q++) {
                    if (y(j, 2) == (q + 1)) {
                        for (int r = 0; r < _R; r++) {
                            _Prob2est(j, r) =(double) _Prob2(q, r);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Pxwzyvar(j, l,h,r) = (double) _Prob1est(j,h)*_Prob2est(j,r)*_Prob32(j, l, h, r);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob4lvar(j, l) = 0;
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Prob4lvar(j, l) += (double) _Pxwzyvar(j, l, h, r);
                        }
                    }
                }
            }

            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    explogL += _Pwy2(k, h) * log(_p_hc[h]);
                }
            }
            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    explogL += _Pzy2(q, r) * log(_p_rc[r]);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            explogL += _Pxwzy(j, l, h, r) * log(_pi_l(l, h, r) * _f(j, l, h, r));
                        }
                    }
                }
            }

            double logLapprox;

            double Entropy = 0;
            for (int k = 0; k < _K; k++) {
                for (int q = 0; q < _Q; q++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            if (empiricalpwzy2(k, q, h, r) != 0) {
                                Entropy += -empiricalpwzy2(k, q, h, r) * log(empiricalpwzy2(k, q, h, r));
                            }
                        }
                    }
                }
            }

            double Entropy1 = 0;
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            if (_Prob32(j, l, h, r) != 0) {
                                Entropy1 += -_Pxwzy(j, l, h, r) * log(_Prob32(j, l, h, r));
                            }
                        }
                    }
                }
            }

            Entropy = Entropy + Entropy1;

            logLapprox = explogL + Entropy;

            explogL2 += explogL;
            Entropy2 += Entropy;
            logLapprox2 += logLapprox;

            //}


            vector<int> _dephr(_I);
            for (int i = 0; i < _I; i++) {
                _dephr[i] = 0;
            }
            //_dephr[0]=1;
            //_dephr[1]=3;
            //_dephr[2]=2;
            //_dephr[3]=2;
            //_dephr[4]=2;
            //_dephr[5]=3;
            //_dephr[6]=2;

            arma::cube _nicl(_I,_C,_L);
            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        _nicl(i, c, l) = 0;
                    }
                }
            }
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int j = 0; j < _n; j++) {
                        _nicl(i,y(j, i + 3),l) += _xjkq(j, l);
                    }
                }
            }


            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        dir_param_p(i, c, l) += _nicl(i,c,l);
                    }
                }
            }
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    arma::vec dir_param_pil(vC[i]);
                    for (int c = 0; c < vC[i]; c++) {
                        dir_param_pil(c) = dir_param_p(i, c, l);
                    }
                    arma::vec sample_dir_pil = rdirichlet(dir_param_pil);
                    for (int c = 0; c < vC[i]; c++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _p(i, c, l, h, r) = sample_dir_pil(c);
                            }
                        }
                    }
                }
            }


            for (int l = 0; l < _L; l++) {
                _den5[l] = 0;
            }
            for (int l = 0; l < _L; l++) {
                for (int j = 0; j < _n; j++) {
                    _den5[l] += _xjkq(j, l);
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _den6(l, r) = 0;
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    for (int j = 0; j < _n; j++) {
                        _den6(l,r) += _zjq(j,r) * _xjkq(j, l);
                    }
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    _den7(l, h) = 0;
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int j = 0; j < _n; j++) {
                        _den7(l, h) += _wjk(j,h) * _xjkq(j, l);
                    }
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int l = 0; l < _L; l++) {
                            for (int c = 0; c < _C; c++) {
                                if (_p(i, c, l, h, r) > 0.9999) {
                                    _p(i, c, l, h, r) = 0.999;
                                }
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < _I; i++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int l = 0; l < _L; l++) {
                            for (int c = 0; c < _C; c++) {
                                if (_p(i, c, l, h, r) < 0.00001) {
                                    _p(i, c, l, h, r) = 0.0001;
                                }
                            }
                        }
                    }
                }
            }



            if (u1 >= _burn && u1 == it[mu]) {
                mu=mu+1;
               for (int h = 0; h < _H; h++) {
                    _p_h_hat(st2, h) += _p_h[h];
                }

                for (int r = 0; r < _R; r++) {
                    _p_r_hat(st2, r) += _p_r[r];
                }

                for (int h = 0; h < _H; h++) {
                    _p_hc_hat(st2,h) += _p_hc[h];
                }

                for (int r = 0; r < _R; r++) {
                    _p_rc_hat(st2,r) += _p_rc[r];
                }

                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _pi_l_hat(st2, l, h, r) += _pi_l(l, h, r);
                        }
                    }
                }

                for (int i = 0; i < _I; i++) {
                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            for (int h = 0; h < _H; h++) {
                                for (int r = 0; r < _R; r++) {
                                    _p_hat(st2, i, c, l, h, r) += _p(i, c, l, h, r);
                                }
                            }
                        }
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _Prob1_hat(k, h) += (double) _Prob1(k, h);
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _Prob2_hat(q, r) += (double) _Prob2(q, r);
                    }
                }

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _Prob3_hat(j, l) += (double) _Prob3(j, l);
                    }
                }

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _Prob32_hat(j, l, h, r) += (double) _Prob32(j, l, h, r);
                            }
                        }
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int q = 0; q < _Q; q++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                empiricalpwzy1(k, q, h, r) += _w(k, h) * _z(q, r);
                            }
                        }
                    }
                }

                for (int i = 0; i < (_I + nlocdep); i++) {
                    for (int c = 0; c < C; c++) {
                        for (int l = 0; l < _L; l++) {
                            for (int h = 0; h < _H; h++) {
                                for (int r = 0; r < _R; r++) {
                                    _p_finale(i, c, l, h, r) = 0;
                                }
                            }
                        }
                    }
                }

                if (nlocdep > 0) {
                    int ar = 0;
                    int ar1 = 0;
                    for (int d = 0; d < nlocdep; d++) {
                        for (int l = 0; l < _L; l++) {
                            for (int h = 0; h < _H; h++) {
                                for (int r = 0; r < _R; r++) {
                                    for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                                        for (int cc = ar; cc < (ar + _vC[locdep(d, 1) - 1]); cc++) {
                                            _p_finale(locdep(d, 0) - 1, c, l, h, r) += _p(d, cc, l, h, r);
                                        }
                                        ar = ar + _vC[locdep(d, 1) - 1];
                                    }
                                    ar = 0;
                                    for (int ccc = 0; ccc < _vC[locdep(d, 1) - 1]; ccc++) {
                                        for (int cccc = ar1; cccc < (_vC[locdep(d, 0) - 1] * _vC[locdep(d, 1) - 1]); cccc++) {
                                            _p_finale(locdep(d, 1) - 1, ccc, l, h, r) += _p(d, cccc, l, h, r);
                                            cccc = cccc + _vC[locdep(d, 1) - 1] - 1;
                                            //cccc = cccc + _vC[locdep(d, 0) - 1];
                                        }
                                        //ar1 = ar1 + _vC[locdep(d, 1) - 1] - 1;
                                        ar1 = ar1 + 1;
                                    }
                                    ar1 = 0;
                                }
                            }
                        }
                    }

                    int dd = 0;
                    int ddd = nlocdep;
                    for (int i = 0; i < (_I + nlocdep); i++) {
                        for (int d = 0; d < nlocdep; d++) {
                            if (locdep(d, 0) != (i + 1) && locdep(d, 1) != (i + 1)) {
                                dd = dd + 1;
                            }
                        }
                        if (dd == nlocdep) {
                            for (int c = 0; c < C; c++) {
                                for (int l = 0; l < _L; l++) {
                                    for (int h = 0; h < _H; h++) {
                                        for (int r = 0; r < _R; r++) {
                                            _p_finale(i, c, l, h, r) = _p(ddd, c, l, h, r);
                                        }
                                    }
                                }
                            }
                            ddd = ddd + 1;
                        }
                        dd = 0;
                    }
                }

                if (nlocdep == 0) {
                    for (int i = 0; i < (_I + nlocdep); i++) {
                        for (int c = 0; c < C; c++) {
                            for (int l = 0; l < _L; l++) {
                                for (int h = 0; h < _H; h++) {
                                    for (int r = 0; r < _R; r++) {
                                        _p_finale(i, c, l, h, r) = _p(i, c, l, h, r);
                                    }
                                }
                            }
                        }
                    }
                }

                //bva
                /*for (int j = 0; j < _n; j++) {
                    yy(j, 0) = y(j, 0);
                }
                for (int j = 0; j < _n; j++) {
                    yy(j, 1) = y(j, 1);
                }
                for (int j = 0; j < _n; j++) {
                    yy(j, 2) = y(j, 2);
                }

                vector<int> samples21(_K);
                //for(int pr=0;pr<50;pr++) {
                std::random_device rdtest;
                std::mt19937 generatorr(rdtest());
                boost::random::discrete_distribution<int> distributionn(_p_hc.begin(), _p_hc.end());
                for (int k = 0; k < _K; k++) {
                    samples21[k] = distributionn(generatorr);
                }

                vector<int> _zj1(_n);
                for (int j = 0; j < _n; j++) {
                    _zj1[j] = samples21[yy(j, 1) - 1];
                }

                //int sumsamples=0;
                vector<int> samples2(_Q);
                //for(int pr=0;pr<50;pr++) {
                //std::random_device rdtest;
                std::mt19937 generatorr2(rdtest());
                boost::random::discrete_distribution<int> distributionn2(_p_rc.begin(), _p_rc.end());
                for (int q = 0; q < _Q; q++) {
                    samples2[q] = distributionn2(generatorr2);
                }

                vector<int> _zj(_n);
                for (int j = 0; j < _n; j++) {
                    _zj[j] = samples2[yy(j, 2) - 1];
                }

                vector<double> Pll(_L);
                vector<int> samples(_n);
                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        Pll[l] = _pi_l(l, _zj1[j], _zj[j]);
                    }
                    std::random_device rdtest;
                    std::mt19937 generator(rdtest());
                    boost::random::discrete_distribution<int> distributionn1(Pll.begin(), Pll.end());
                    //for (int j = 0; j < _n; j++) {
                    samples[j] = distributionn1(generator);
                }

                //vector<double> Pi(C);
                int samp;
                for (int in = 0; in < (_I + nlocdep); in++) {
                    vector<double> Pi(_vC[in]);
                    for (int j = 0; j < _n; j++) {
                        for (int c = 0; c < _vC[in]; c++) {
                            Pi[c] = _p_finale(in, c, samples[j], _zj1[j], _zj[j]);
                        }
                        std::random_device rdtest;
                        std::mt19937 generator(rdtest());
                        boost::random::discrete_distribution<int> distribution2(Pi.begin(), Pi.end());
                        samp = distribution2(generator);
                        yy(j, in + 3) = samp;
                    }
                }

                //BVA
                for (int i1 = 0; i1 < (_I + nlocdep); i1++) {
                    for (int i2 = 0; i2 < (_I + nlocdep); i2++) {
                        BVA(st2, vp, i1, i2) = 0;
                    }
                }
                for (int i1 = 0; i1 < (_I + nlocdep); i1++) {
                    for (int i2 = i1 + 1; i2 < (_I + nlocdep); i2++) {

                        for (int c = 0; c < C; c++) {
                            for (int cc = 0; cc < C; cc++) {
                                n(c, cc) = 0;
                            }
                        }
                        for (int c = 0; c < C; c++) {
                            for (int cc = 0; cc < C; cc++) {
                                for (int j = 0; j < _n; j++) {
                                    if (yy(j, i1 + 3) == c && yy(j, i2 + 3) == cc) {
                                        n(c, cc) += 1;
                                    }
                                }
                            }
                        }

                        vector<int> nr(C);
                        vector<int> nc(C);
                        for (int c = 0; c < C; c++) {
                            nr[c] = 0;
                        }
                        for (int cc = 0; cc < C; cc++) {
                            nc[cc] = 0;
                        }
                        for (int c = 0; c < C; c++) {
                            for (int cc = 0; cc < C; cc++) {
                                nr[c] += n(c, cc);
                            }
                        }
                        for (int cc = 0; cc < C; cc++) {
                            for (int c = 0; c < C; c++) {
                                nc[cc] += n(c, cc);
                            }
                        }
                        for (int c = 0; c < C; c++) {
                            for (int cc = 0; cc < C; cc++) {
                                m0(c, cc) = 0;
                            }
                        }
                        for (int c = 0; c < C; c++) {
                            for (int cc = 0; cc < C; cc++) {
                                m0(c, cc) = (double) nr[c] * nc[cc] / _n;
                            }
                        }

                        for (int c = 0; c < _vC[i1]; c++) {
                            for (int cc = 0; cc < _vC[i2]; cc++) {
                                BVA(st2, vp, i1, i2) += (double) pow((n(c, cc) - m0(c, cc)), 2) / m0(c, cc);
                            }
                        }

                        BVA(st2, vp, i1, i2) = (double) BVA(st2, vp, i1, i2) / ((_vC[i1] - 1) * (_vC[i2] - 1));

                    }
                }

                for (int i = 0; i < (_I + nlocdep); i++) {
                    BVAG(st2, vp, i) = 0;
                }

                for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

                    for (int k = 0; k < _K; k++) {
                        for (int cc = 0; cc < C; cc++) {
                            n1(k, cc) = 0;
                        }
                    }
                    for (int c = 0; c < _K; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            for (int j = 0; j < _n; j++) {
                                if (yy(j, 1) - 1 == c && yy(j, i1 + 3) == cc) {
                                    n1(c, cc) += 1;
                                }
                            }
                        }
                    }

                    vector<int> nr1(_K);
                    vector<int> nc1(C);
                    for (int c = 0; c < _K; c++) {
                        nr1[c] = 0;
                    }
                    for (int cc = 0; cc < C; cc++) {
                        nc1[cc] = 0;
                    }
                    for (int c = 0; c < _K; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            nr1[c] += n1(c, cc);
                        }
                    }

                    for (int cc = 0; cc < C; cc++) {
                        for (int c = 0; c < _K; c++) {
                            nc1[cc] += n1(c, cc);
                        }
                    }

                    for (int c = 0; c < _K; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            m1(c, cc) = 0;
                        }
                    }
                    for (int c = 0; c < _K; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            m1(c, cc) = (double) nr1[c] * nc1[cc] / _n;
                        }
                    }

                    for (int c = 0; c < _K; c++) {
                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            BVAG(st2, vp, i1) += (double) pow((n1(c, cc) - m1(c, cc)), 2) / m1(c, cc);
                        }
                    }

                    BVAG(st2, vp, i1) = (double) BVAG(st2, vp, i1) / ((_K - 1) * (_vC[i1] - 1));

                }

                for (int i = 0; i < (_I + nlocdep); i++) {
                    BVAG1(st2, vp, i) = 0;
                }

                for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

                    for (int q = 0; q < _Q; q++) {
                        for (int cc = 0; cc < C; cc++) {
                            n11(q, cc) = 0;
                        }
                    }
                    for (int c = 0; c < _Q; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            for (int j = 0; j < _n; j++) {
                                if (yy(j, 2) - 1 == c && yy(j, i1 + 3) == cc) {
                                    n11(c, cc) += 1;
                                }
                            }
                        }
                    }

                    vector<int> nr11(_Q);
                    vector<int> nc11(C);
                    for (int c = 0; c < _Q; c++) {
                        nr11[c] = 0;
                    }
                    for (int cc = 0; cc < C; cc++) {
                        nc11[cc] = 0;
                    }
                    for (int c = 0; c < _Q; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            nr11[c] += n11(c, cc);
                        }
                    }

                    for (int cc = 0; cc < C; cc++) {
                        for (int c = 0; c < _Q; c++) {
                            nc11[cc] += n11(c, cc);
                        }
                    }

                    for (int c = 0; c < _Q; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            m11(c, cc) = 0;
                        }
                    }
                    for (int c = 0; c < _Q; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            m11(c, cc) = (double) nr11[c] * nc11[cc] / _n;
                        }
                    }

                    for (int c = 0; c < _Q; c++) {
                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            BVAG1(st2, vp, i1) += (double) pow((n11(c, cc) - m11(c, cc)), 2) / m11(c, cc);
                        }
                    }

                    BVAG1(st2, vp, i1) = (double) BVAG1(st2, vp, i1) / ((_Q - 1) * (_vC[i1] - 1));

                }

                for (int i = 0; i < (_I + nlocdep); i++) {
                    BVAP(st2, vp, i) = 0;
                }

                vector<int> vetcoppie(_K);
                for (int k = 0; k < _K; k++) {
                    vetcoppie[k] = _nk[k] * (_nk[k] - 1) / 2;
                }
                int totcoppie = 0;
                for (int k = 0; k < _K; k++) {
                    totcoppie += vetcoppie[k];
                }

                for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

                    mat coppie(_K, _n);
                    for (int k = 0; k < _K; k++) {
                        for (int j = 0; j < _n; j++) {
                            coppie(k, j) = -1;
                        }
                    }
                    int co = 0;
                    for (int k = 0; k < _K; k++) {
                        for (int j = 0; j < _n; j++) {
                            if (yy(j, 1) - 1 == k) {
                                coppie(k, co) = yy(j, i1 + 3);
                                co = co + 1;
                            }
                        }
                        co = 0;
                    }

                    for (int c = 0; c < C; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            n2(c, cc) = 0;
                        }
                    }

                    for (int k = 0; k < _K; k++) {

                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            for (int c = 0; c < _vC[i1]; c++) {
                                for (int a = 0; a < _nk[k]; a++) {
                                    if (coppie(k, a) == cc) {
                                        for (int b = (a + 1); b < _nk[k]; b++) {
                                            if (coppie(k, b) == c) {
                                                n2(cc, c) += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }

                    vector<int> nr2(_vC[i1]);
                    vector<int> nc2(_vC[i1]);
                    for (int c = 0; c < _vC[i1]; c++) {
                        nr2[c] = 0;
                    }
                    for (int cc = 0; cc < _vC[i1]; cc++) {
                        nc2[cc] = 0;
                    }
                    for (int c = 0; c < _vC[i1]; c++) {
                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            nr2[c] += n2(c, cc);
                        }
                    }

                    for (int cc = 0; cc < _vC[i1]; cc++) {
                        for (int c = 0; c < _vC[i1]; c++) {
                            nc2[cc] += n2(c, cc);
                        }
                    }

                    for (int c = 0; c < C; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            m2(c, cc) = 0;
                        }
                    }
                    for (int c = 0; c < _vC[i1]; c++) {
                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            m2(c, cc) = (double) nr2[c] * nc2[cc] / totcoppie;
                        }
                    }

                    for (int c = 0; c < _vC[i1]; c++) {
                        for (int cc = c; cc < _vC[i1]; cc++) {
                            if (c == cc) {
                                BVAP(st2, vp, i1) += (double) pow((n2(c, cc) - m2(c, cc)), 2) / m2(c, cc);
                            }
                            if (c != cc) {
                                BVAP(st2, vp, i1) +=
                                        (double) pow(((n2(c, cc) + n2(cc, c)) - (m2(c, cc) + m2(cc, c))), 2) /
                                        (m2(c, cc) + m2(cc, c));
                            }
                        }
                    }

                    BVAP(st2, vp, i1) = (double) (BVAP(st2, vp, i1) * _K) / (_n * _vC[i1] * (_vC[i1] - 1) / 2);

                }

                for (int i = 0; i < (_I + nlocdep); i++) {
                    BVAP1(st2, vp, i) = 0;
                }

                vector<int> vetcoppie1(_Q);
                for (int q = 0; q < _Q; q++) {
                    vetcoppie1[q] = _nq[q] * (_nq[q] - 1) / 2;
                }
                int totcoppie1 = 0;
                for (int q = 0; q < _Q; q++) {
                    totcoppie1 += vetcoppie1[q];
                }

                for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

                    mat coppie1(_Q, _n);
                    for (int q = 0; q < _Q; q++) {
                        for (int j = 0; j < _n; j++) {
                            coppie1(q, j) = -1;
                        }
                    }
                    int co = 0;
                    for (int q = 0; q < _Q; q++) {
                        for (int j = 0; j < _n; j++) {
                            if (yy(j, 2) - 1 == q) {
                                coppie1(q, co) = yy(j, i1 + 3);
                                co = co + 1;
                            }
                        }
                        co = 0;
                    }

                    for (int c = 0; c < C; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            n22(c, cc) = 0;
                        }
                    }

                    for (int q = 0; q < _Q; q++) {

                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            for (int c = 0; c < _vC[i1]; c++) {
                                for (int a = 0; a < _nq[q]; a++) {
                                    if (coppie1(q, a) == cc) {
                                        for (int b = (a + 1); b < _nq[q]; b++) {
                                            if (coppie1(q, b) == c) {
                                                n22(cc, c) += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }

                    vector<int> nr22(_vC[i1]);
                    vector<int> nc22(_vC[i1]);
                    for (int c = 0; c < _vC[i1]; c++) {
                        nr22[c] = 0;
                    }
                    for (int cc = 0; cc < _vC[i1]; cc++) {
                        nc22[cc] = 0;
                    }
                    for (int c = 0; c < _vC[i1]; c++) {
                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            nr22[c] += n22(c, cc);
                        }
                    }

                    for (int cc = 0; cc < _vC[i1]; cc++) {
                        for (int c = 0; c < _vC[i1]; c++) {
                            nc22[cc] += n22(c, cc);
                        }
                    }

                    for (int c = 0; c < C; c++) {
                        for (int cc = 0; cc < C; cc++) {
                            m22(c, cc) = 0;
                        }
                    }
                    for (int c = 0; c < _vC[i1]; c++) {
                        for (int cc = 0; cc < _vC[i1]; cc++) {
                            m22(c, cc) = (double) nr22[c] * nc22[cc] / totcoppie1;
                        }
                    }

                    for (int c = 0; c < _vC[i1]; c++) {
                        for (int cc = c; cc < _vC[i1]; cc++) {
                            if (c == cc) {
                                BVAP1(st2, vp, i1) += (double) pow((n22(c, cc) - m22(c, cc)), 2) / m22(c, cc);
                            }
                            if (c != cc) {
                                BVAP1(st2, vp, i1) +=
                                        (double) pow(((n22(c, cc) + n22(cc, c)) - (m22(c, cc) + m22(cc, c))), 2) /
                                        (m22(c, cc) + m22(cc, c));
                            }
                        }
                    }

                    BVAP1(st2, vp, i1) = (double) (BVAP1(st2, vp, i1) * _Q) / (_n * _vC[i1] * (_vC[i1] - 1) / 2);

                }*/

                vp = vp + 1;

            }


            if (u1 == _iter - 1) {

                //_p_h_hat
                for (int h = 0; h < _H; h++) {
                    _p_h_hat(st2, h) = (double) _p_h_hat(st2, h) / (nn);
                }

                for (int h = 0; h < _H; h++) {
                    _p_hc_hat(st2,h) = (double) _p_hc_hat(st2,h) / (nn);
                }
                for (int h = 0; h < _H; h++) {
                    _p_hcc_hat[h] = (double) _p_hc_hat(st2,h);
                }

                //_p_r_hat
                for (int r = 0; r < _R; r++) {
                    _p_r_hat(st2, r) = (double) _p_r_hat(st2, r) / (nn);
                }

                for (int r = 0; r < _R; r++) {
                    _p_rc_hat(st2,r) = (double) _p_rc_hat(st2,r) / (nn);
                }
                for (int r = 0; r < _R; r++) {
                    _p_rcc_hat[r] = (double) _p_rc_hat(st2,r);
                }

                //_pi_l_hat
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _pi_l_hat(st2, l, h, r) = (double) _pi_l_hat(st2, l, h, r) / (nn);
                        }
                    }
                }

                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _pi_lc_hat(l, h, r) = (double) _pi_l_hat(st2, l, h, r);
                        }
                    }
                }


                /*for (int i = 0; i < (_H - 1); i++) {
                    min1 = i;
                    for (int j = (i + 1); j < _H; j++) {
                        if (_p_hcc_hat[j] < _p_hcc_hat[min1]) {
                            min1 = j;
                        }
                    }

                    temp1 = _p_hcc_hat[i];
                    _p_hcc_hat[i] = _p_hcc_hat[min1];
                    _p_hcc_hat[min1] = temp1;

                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            temp5(l, r) = (double) _pi_l_hat(st2, l, i, r);
                            _pi_l_hat(st2, l, i, r) =(double) _pi_l_hat(st2, l, min1, r);
                            _pi_l_hat(st2, l, min1, r) = (double) temp5(l, r);
                        }
                    }

                }

                for (int i = 0; i < (_R - 1); i++) {
                    min2 = i;
                    for (int j = (i + 1); j < _R; j++) {
                        if (_p_rcc_hat[j] < _p_rcc_hat[min2]) {
                            min2 = j;
                        }
                    }

                    temp6 = _p_rcc_hat[i];
                    _p_rcc_hat[i] = _p_rcc_hat[min2];
                    _p_rcc_hat[min2] = temp6;

                    for (int l = 0; l < _L; l++) {
                        for (int h = 0; h < _H; h++) {
                            temp10(l, h) = (double) _pi_l_hat(st2,l, h, i);
                            _pi_l_hat(st2, l, h, i) = _pi_l_hat(st2, l, h, min2);
                            _pi_l_hat(st2, l, h, min2) = temp10(l, h);
                        }
                    }

                }*/


                //_p_hat
                for (int i = 0; i < _I; i++) {
                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            for (int h = 0; h < _H; h++) {
                                for (int r = 0; r < _R; r++) {
                                    _p_hat(st2, i, c, l, h, r) =
                                            (double) _p_hat(st2, i, c, l, h, r) / (nn);
                                }
                            }
                        }
                    }
                }

                if (nlocdep > 0) {
                    int ar = 0;
                    int ar1 = 0;
                    for (int d = 0; d < nlocdep; d++) {
                        for (int l = 0; l < _L; l++) {
                            for (int h = 0; h < _H; h++) {
                                for (int r = 0; r < _R; r++) {
                                    for (int c = 0; c < _vC[locdep(d, 0) - 1]; c++) {
                                        for (int cc = ar; cc < (ar + _vC[locdep(d, 1) - 1]); cc++) {
                                            _p_finale2(st2, locdep(d, 0) - 1, c, l, h, r) += _p_hat(st2, d, cc, l, h,
                                                                                                    r);
                                        }
                                        ar = ar + _vC[locdep(d, 1) - 1];
                                    }
                                    ar = 0;
                                    for (int ccc = 0; ccc < _vC[locdep(d, 1) - 1]; ccc++) {
                                        for (int cccc = ar1;
                                             cccc < (_vC[locdep(d, 0) - 1] * _vC[locdep(d, 1) - 1]); cccc++) {
                                            _p_finale2(st2, locdep(d, 1) - 1, ccc, l, h, r) += _p_hat(st2, d, cccc, l,
                                                                                                      h, r);
                                            cccc = cccc + _vC[locdep(d, 1) - 1] - 1;
                                            //cccc = cccc + _vC[locdep(d, 0) - 1];
                                        }
                                        //ar1 = ar1 + _vC[locdep(d, 1) - 1] - 1;
                                        ar1 = ar1 + 1;
                                    }
                                    ar1 = 0;
                                }
                            }
                        }
                    }

                    int dd = 0;
                    int ddd = nlocdep;
                    for (int i = 0; i < (_I + nlocdep); i++) {
                        for (int d = 0; d < nlocdep; d++) {
                            if (locdep(d, 0) != (i + 1) && locdep(d, 1) != (i + 1)) {
                                dd = dd + 1;
                            }
                        }
                        if (dd == nlocdep) {
                            for (int c = 0; c < C; c++) {
                                for (int l = 0; l < _L; l++) {
                                    for (int h = 0; h < _H; h++) {
                                        for (int r = 0; r < _R; r++) {
                                            _p_finale2(st2, i, c, l, h, r) = _p_hat(st2, ddd, c, l, h, r);
                                        }
                                    }
                                }
                            }
                            ddd = ddd + 1;
                        }
                        dd = 0;
                    }
                }

                if (nlocdep == 0) {
                    for (int i = 0; i < (_I + nlocdep); i++) {
                        for (int c = 0; c < C; c++) {
                            for (int l = 0; l < _L; l++) {
                                for (int h = 0; h < _H; h++) {
                                    for (int r = 0; r < _R; r++) {
                                        _p_finale2(st2, i, c, l, h, r) = _p_hat(st2, i, c, l, h, r);
                                    }
                                }
                            }
                        }
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _Prob1_hat(k, h) = (double) _Prob1_hat(k, h) / (nn);
                    }
                }

                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _Prob2_hat(q, r) = (double) _Prob2_hat(q, r) / (nn);
                    }
                }

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _Prob3_hat(j, l) = (double) _Prob3_hat(j, l) / (nn);
                    }
                }

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _Prob32_hat(j, l, h, r) = (double) _Prob32_hat(j, l, h, r) / (nn);
                            }
                        }
                    }
                }

                for (int k = 0; k < _K; k++) {
                    for (int q = 0; q < _Q; q++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                empiricalpwzy1(k, q, h, r) = (double) empiricalpwzy1(k, q, h, r) / (nn);
                            }
                        }
                    }
                }

                explogL2 = (double) explogL2 / (nn);
                Entropy2 = (double) Entropy2 / (nn);
                logLapprox2 = (double) logLapprox2 / (nn);


            }
            for (int h = 0; h < _H; h++) {
                _p_h[h] = _p_hc[h];
            }
            for (int r = 0; r < _R; r++) {
                _p_r[r] = _p_rc[r];
            }

            //cout << u1 + 1 << endl;
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _Pwy2(k, h) = 0;
            }
        }
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _Pzy2(q, r) = 0;
            }
        }
        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _Pwy2(k, h) += empiricalpwzy1(k, 0, h, r);
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int h = 0; h < _H; h++) {
                    _Pzy2(q, r) += empiricalpwzy1(0, q, h, r);
                }
            }
        }

        double Entropy = 0;
        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        if (empiricalpwzy1(k, q, h, r) != 0) {
                            Entropy += -empiricalpwzy1(k, q, h, r) * log(empiricalpwzy1(k, q, h, r));
                        }
                    }
                }
            }
        }

        Tensor<double, 4> _Pxwzy(_n, _L, _H, _R);
        arma::cube _Pwzyest(_n, _H, _R);

        for (int j = 0; j < _n; j++) {
            for (int k = 0; k < _K; k++) {
                for (int q = 0; q < _Q; q++) {
                    if (y(j, 1) == (k + 1) && y(j, 2) == (q + 1)) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _Pwzyest(j, h, r) = empiricalpwzy1(k, q, h, r);
                            }
                        }
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _Pxwzy(j, l, h, r) = _Pwzyest(j, h, r) * _Prob3_hat(j, l);
                    }
                }
            }
        }

        double Entropy1 = 0;
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        if (_Prob32_hat(j, l, h, r) != 0) {
                            Entropy1 += -_Pxwzy(j, l, h, r) * log(_Prob32_hat(j, l, h, r));
                        }
                    }
                }
            }
        }

        Entropy = Entropy + Entropy1;

        double explogL = 0;
        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                explogL += _Pwy2(k, h) * log(_p_hc_hat[h]);
            }
        }
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                explogL += _Pzy2(q, r) * log(_p_rc_hat[r]);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        explogL += _Pxwzy(j, l, h, r) * log(_pi_lc_hat(l, h, r) * _f(j, l, h, r));
                    }
                }
            }
        }

        logLapprox2 = explogL + Entropy;

        logLtot[st2] = logLapprox2;
        explogLtot[st2] = explogL; //2

        //cout<<st2+1<<endl;
    }

    _maxlogl = logLtot[0];
    indicestarting = 0;
    for (int st2 = 0; st2 < _starting2; st2++) {
        if (logLtot[st2] > _maxlogl) {
            _maxlogl = logLtot[st2];
            indicestarting = st2;
        }
    }

    if(_starting2==1){
        indicestarting=0;
    }

    logLapprox2 = logLtot[indicestarting];
    explogL2 = explogLtot[indicestarting];

    vector<double> phf(_H);
    vector<double> prf(_R);
    vector<double> phcf(_H);
    vector<double> prcf(_R);
    arma::cube plf(_L, _H, _R);
    Tensor<double, 5> pif(_I + nlocdep,C, _L, _H, _R);
    for (int h = 0; h < _H; h++) {
        phf[h] = (double) _p_h_hat(indicestarting, h);
    }
    for (int r = 0; r < _R; r++) {
        prf[r] = (double) _p_r_hat(indicestarting, r);
    }
    for (int h = 0; h < _H; h++) {
        phcf[h] = (double) _p_hc_hat(indicestarting, h);
    }
    for (int r = 0; r < _R; r++) {
        prcf[r] = (double) _p_rc_hat(indicestarting, r);
    }
    for (int l = 0; l < _L; l++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                plf(l, h, r) = (double) _pi_l_hat(indicestarting, l, h, r);
            }
        }
    }
    for(int i=0;i<(_I+nlocdep);i++) {
        for (int c = 0; c < C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        pif(i, c, l, h, r) = (double) _p_finale2(indicestarting, i, c, l, h, r);
                    }
                }
            }
        }
    }


    cout<<"ph"<<endl;
    for (int h = 0; h < _H; h++) {
        cout<<phcf[h]<<",";
    }
    cout<<"\n";
    cout<<"pr"<<endl;
    for (int r = 0; r < _R; r++) {
        cout<<prcf[r]<<",";
    }
    cout<<"\n";
    cout<<"phc"<<endl;
    for (int h = 0; h < _H; h++) {
        cout<<phcf[h]<<",";
    }
    cout<<"\n";
    cout<<"prc"<<endl;
    for (int r = 0; r < _R; r++) {
        cout<<prcf[r]<<",";
    }
    cout<<"\n";
    cout<<"pl(l,0,0)"<<endl;
    for (int l = 0; l < _L; l++) {
        cout<<plf(l,0,0)<<",";
    }
    cout<<"\n";
    cout<<"pl(l,0,1)"<<endl;
    for (int l = 0; l < _L; l++) {
        cout<<plf(l,0,1)<<",";
    }
    cout<<"\n";
    cout<<"pl(l,1,0)"<<endl;
    for (int l = 0; l < _L; l++) {
        cout<<plf(l,1,0)<<",";
    }
    cout<<"\n";
    cout<<"pl(l,1,1)"<<endl;
    for (int l = 0; l < _L; l++) {
        cout<<plf(l,1,1)<<",";
    }
    cout<<"\n";
    cout<<"pl(l,2,0)"<<endl;
    for (int l = 0; l < _L; l++) {
        cout<<plf(l,2,0)<<",";
    }
    cout<<"\n";
    cout<<"pl(l,2,1)"<<endl;
    for (int l = 0; l < _L; l++) {
        cout<<plf(l,2,1)<<",";
    }
    cout<<"\n";
    cout<<"p1"<<endl;
    for (int c = 0; c < _vC[0]; c++) {
        for (int l = 0; l < _L; l++) {
            cout<<pif(0, c, l, 0, 0)<<",";
        }
        cout<<"\n";
    }
    cout<<"\n";
    cout<<"p2"<<endl;
    for (int c = 0; c < _vC[1]; c++) {
        for (int l = 0; l < _L; l++) {
            cout<<pif(1, c, l, 0, 0)<<",";
        }
        cout<<"\n";
    }
    cout<<"\n";
    cout<<"p3"<<endl;
    for (int c = 0; c < _vC[2]; c++) {
        for (int l = 0; l < _L; l++) {
            cout<<pif(2, c, l, 0, 0)<<",";
        }
        cout<<"\n";
    }
    cout<<"\n";
    cout<<"p4"<<endl;
    for (int c = 0; c < _vC[3]; c++) {
        for (int l = 0; l < _L; l++) {
            cout<<pif(3, c, l, 0, 0)<<",";
        }
        cout<<"\n";
    }
    cout<<"\n";
    cout<<"p5"<<endl;
    for (int c = 0; c < _vC[4]; c++) {
        for (int l = 0; l < _L; l++) {
            cout<<pif(4, c, l, 0, 0)<<",";
        }
        cout<<"\n";
    }
    cout<<"\n";
    cout<<"p6"<<endl;
    for (int c = 0; c < _vC[5]; c++) {
        for (int l = 0; l < _L; l++) {
            cout<<pif(5, c, l, 0, 0)<<",";
        }
        cout<<"\n";
    }
    cout<<"\n";

    mat pbva((_I+nlocdep), (_I+nlocdep));
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        for (int i2 = 0; i2 < (_I+nlocdep); i2++) {
            pbva(i1, i2) = 0;
        }
    }
    for (int in = 0; in < (_iter - _burn - nn); in++) {
        for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
            for (int i2 = 0; i2 < (_I+nlocdep); i2++) {
                if (BVA(indicestarting,in, i1, i2) > bva(i1, i2)) {
                    pbva(i1, i2) += 1;
                }
            }
        }
    }
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        for (int i2 = 0; i2 < (_I+nlocdep); i2++) {
            pbva(i1, i2) = (double) pbva(i1, i2) / (_iter - _burn - nn);
        }
    }

    vector<double> pbvag((_I+nlocdep));
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvag[i1] = 0;
    }
    for (int in = 0; in < (_iter - _burn - nn); in++) {
        for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
            if (BVAG(indicestarting,in, i1) > bvag[i1]) {
                pbvag[i1] += 1;
            }
        }
    }
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvag[i1] = (double) pbvag[i1] / (_iter - _burn - nn);
    }



    vector<double> pbvag1((_I+nlocdep));
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvag1[i1] = 0;
    }
    for (int in = 0; in < (_iter - _burn - nn); in++) {
        for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
            if (BVAG1(indicestarting,in, i1) > bvag1[i1]) {
                pbvag1[i1] += 1;
            }
        }
    }
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvag1[i1] = (double) pbvag1[i1] / (_iter - _burn - nn);
    }

    vector<double> pbvap((_I+nlocdep));
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvap[i1] = 0;
    }
    for (int in = 0; in < (_iter - _burn - nn); in++) {
        for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
            if (BVAP(indicestarting,in, i1) > bvap[i1]) {
                pbvap[i1] += 1;
            }
        }
    }
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvap[i1] = (double) pbvap[i1] / (_iter - _burn - nn);
    }


    vector<double> pbvap1((_I+nlocdep));
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvap1[i1] = 0;
    }
    for (int in = 0; in < (_iter - _burn - nn); in++) {
        for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
            if (BVAP1(indicestarting,in, i1) > bvap1[i1]) {
                pbvap1[i1] += 1;
            }
        }
    }
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        pbvap1[i1] = (double) pbvap1[i1] / (_iter - _burn - nn);
    }

    /*cout << "pbva" << endl;
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        for (int i2 = 0; i2 < (_I+nlocdep); i2++) {
            cout << pbva(i1, i2) << ",";
        }
        cout << "\n";
    }
    cout << "\n";
    cout << "pbvag" << endl;
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        cout << pbvag[i1] << ",";
    }
    cout << "\n";
    cout << "pbvag1" << endl;
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        cout << pbvag1[i1] << ",";
    }
    cout << "\n";
    cout << "pbvap" << endl;
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        cout << pbvap[i1] << ",";
    }
    cout << "\n";
    cout << "pbvap1" << endl;
    for (int i1 = 0; i1 < (_I+nlocdep); i1++) {
        cout << pbvap1[i1] << ",";
    }
    cout << "\n";*/


    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _f(j, l, h, r) = 1;
                }
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int i = 0; i < _I; i++) {
                        _f(j, l, h, r) = _f(j, l, h, r) * _p_hat(indicestarting,i, y2(j, i + 3), l, h, r);
                    }
                }
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _pf(j, l, h, r) = (double) plf(l, h, r) * _f(j, l, h, r);
                }
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _fykq(j, h, r) = 0;
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _fykq(j, h, r) += _pf(j, l, h, r);
                }
            }
        }
    }

    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _res(k, q, h, r) = 1;
                }
            }
        }
    }


    int ml = 0;
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int j = ml; j < _nkq(k, q) + ml; j++) {
                        _res(k, q, h, r) = _res(k, q, h, r) * _fykq(j, h, r);
                    }
                }
            }
            ml = ml + _nkq(k, q);
        }
    }

    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            _tt(k, q) = 0;
        }
    }

    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    if (_res(k, q, h, r) != 1) {
                        _tt(k, q) += phf[h] * prf[r] * _res(k, q, h, r);
                    }
                }
            }
        }
    }

    _maxL = 0;
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            if (_tt(k, q) != 0) {
                _maxL += log(_tt(k, q));
            }
        }
    }

    //Entropia
    double entropyhigh;
    double entropylow;
    double numentr1 = 0;
    double numentr2 = 0;
    double denentr1 = 0;
    double denentr2 = 0;

    //-------------------------------------------------------------------

    // classification statistics

    //entropy
    double centropy = 0;
    double centropy2 = 0;
    double centropy3 = 0;
    int _pcv2;
    vector<int> _cv2(_R);
    int _pcv;
    vector<int> _cv(_H);
    //classificazione
    while (_pcv2 == 0) {
        for (int r = 0; r < _R; r++) {
            _cv2[r] = 0;
        }
        _pcv2 = 1;
        _z.zeros(_Q, _R);
        boost::random::discrete_distribution<int> distribution2(prf.begin(), prf.end());
        for (int q = 0; q < _Q; q++) {
            int sample = distribution2(generator);
            _z(q, sample) = 1;
        }

        for (int r = 0; r < _R; r++) {
            for (int q = 0; q < _Q; q++) {
                _cv2[r] += _z(q, r);
            }
        }
        for (int r = 0; r < _R; r++) {
            _pcv2 = _pcv2 * _cv2[r];
        }
    }

    Tensor<double, 4> empiricalpwzy(_K, _Q, _H, _R);
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    empiricalpwzy(k, q, h, r) = 0;
                }
            }
        }
    }

    mat _Prob1_hat1(_K, _H);
    mat _Prob2_hat1(_Q, _R);
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            _Prob1_hat1(k, h) = 0;
        }
    }
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _Prob2_hat1(q, r) = 0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _f(j, l, h, r) = 1;
                }
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int i = 0; i < _I; i++) {
                        _f(j, l, h, r) = _f(j, l, h, r) * _p_hat(indicestarting,i, y2(j, i + 3), l, h, r);
                    }
                }
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _pf(j, l, h, r) = (double) plf(l, h, r) * _f(j, l, h, r);
                }
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _fykq(j, h, r) = 0;
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _fykq(j, h, r) += _pf(j, l, h, r);
                }
            }
        }
    }

    _w_hat.zeros(_K, _H);
    _z_hat.zeros(_Q, _R);
    _x_hat.zeros(_n, _L);

    for (int gi = 0; gi < _iter3; gi++) {

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _ff(k, h) = 0;
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _res(k, q, h, r) = 0;
                    }
                }
            }
        }


        int m = 0;
        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = m; j < _nkq(k, q) + m; j++) {
                            _res(k, q, h, r) += log(_fykq(j, h, r));
                        }
                    }
                }
                m = m + _nkq(k, q);
            }
        }


        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _res(k, q, h, r) = _z(q, r) * _res(k, q, h, r);
                    }
                }
            }
        }


        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int q = 0; q < _Q; q++) {
                        _ff(k, h) += _res(k, q, h, r);
                    }
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _num2(k, h) = log(phf[h]) + _ff(k, h);
            }
        }

        vector<double> _maxh(_K);
        for (int k = 0; k < _K; k++) {
            _maxh[k] = _num2(k, 0);
            for (int h = 0; h < _H; h++) {
                if (_num2(k, h) > _maxh[k]) {
                    _maxh[k] = _num2(k, h);
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _num2(k, h) = _num2(k, h) - _maxh[k];
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _num2(k, h) = exp(_num2(k, h));
            }
        }


        for (int k = 0; k < _K; k++) {
            _den2[k] = 0;
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _den2[k] += _num2(k, h);
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _Prob1(k, h) = (double) _num2(k, h) / _den2[k];
            }
        }

        _pcv = 0;
        while (_pcv == 0) {
            _pcv = 1;
            for (int h = 0; h < _H; h++) {
                _cv[h] = 0;
            }
            _w.zeros(_K, _H);
//boost::mt19937 generator2(1234);
            for (int k = 0; k < _K; k++) {
                boost::random::discrete_distribution<int> distribution7(_Prob1.row(k).begin(), _Prob1.row(k).end());
                int sample = distribution7(generator);
                _w(k, sample) = 1;
            }
            for (int h = 0; h < _H; h++) {
                for (int k = 0; k < _K; k++) {
                    _cv[h] += _w(k, h);
                }
            }
            for (int h = 0; h < _H; h++) {
                _pcv = _pcv * _cv[h];
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _w_hat(k, h) += _w(k, h);
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 0;
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _res2(q, k, h, r) = 0;
                    }
                }
            }
        }

        int mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = mm; j < _nkq(k, q) + mm; j++) {
                            _res2(q, k, h, r) += log(_fykq(_indvecqk[j], h, r));
                        }
                    }
                }
                mm = mm + _nkq(k, q);
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _res2(q, k, h, r) = _w(k, h) * _res2(q, k, h, r);
                    }
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _ff2(q, r) += _res2(q, k, h, r);
                    }
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = log(prf[r]) + _ff2(q, r);
            }
        }

        vector<double> _maxr(_Q);
        for (int q = 0; q < _Q; q++) {
            _maxr[q] = _num3(q, 0);
            for (int r = 0; r < _R; r++) {
                if (_num3(q, r) > _maxr[q]) {
                    _maxr[q] = _num3(q, r);
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = _num3(q, r) - _maxr[q];
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = exp(_num3(q, r));
            }
        }

        for (int q = 0; q < _Q; q++) {
            _den3[q] = 0;
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _den3[q] += _num3(q, r);
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _Prob2(q, r) = (double) _num3(q, r) / _den3[q];
            }
        }

        _pcv2 = 0;
        while (_pcv2 == 0) {
            _pcv2 = 1;
            for (int r = 0; r < _R; r++) {
                _cv2[r] = 0;
            }
            _z.zeros(_Q, _R);
            for (int q = 0; q < _Q; q++) {
                boost::random::discrete_distribution<int> distribution8(_Prob2.row(q).begin(), _Prob2.row(q).end());
                int sample = distribution8(generator);
                _z(q, sample) = 1;
            }

            for (int r = 0; r < _R; r++) {
                for (int q = 0; q < _Q; q++) {
                    _cv2[r] += _z(q, r);
                }
            }
            for (int r = 0; r < _R; r++) {
                _pcv2 = _pcv2 * _cv2[r];
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _z_hat(q, r) += _z(q, r);
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        empiricalpwzy(k, q, h, r) += _w(k, h) * _z(q, r);
                    }
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _Prob1_hat1(k, h) += _Prob1(k, h);
            }
        }
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _Prob2_hat1(q, r) += _Prob2(q, r);
            }
        }

    }

    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            _Prob1_hat1(k, h) = _Prob1_hat1(k, h) / _iter3;
        }
    }
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _Prob2_hat1(q, r) = _Prob2_hat1(q, r) / _iter3;
        }
    }

    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    empiricalpwzy(k, q, h, r) = (double) empiricalpwzy(k, q, h, r) / _iter3;
                }
            }
        }
    }

    denentr2=0;
    for (int h = 0; h < _H; h++) {
        for (int r = 0; r < _R; r++) {
            denentr2 += -1 * phf[h] * prf[r] * (log(prf[r]) + log(phf[h]));
        }
    }
    denentr2 = denentr2 * ckq;

    numentr2 = 0;
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    if (_tt(k, q) != 0 && empiricalpwzy(k, q, h, r) != 0) {
                        //if(empiricalpwzy(k, q, h, r)!=0) {
                        numentr2 += (double) -1 * empiricalpwzy(k, q, h, r) * log(empiricalpwzy(k, q, h, r));
                    }
                }
            }
        }
    }
    entropyhigh = 1 - (double) numentr2 / (denentr2);


    numentr2 = 0;
    for (int k = 0; k < _K; k++) {
        for (int q = 0; q < _Q; q++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    if (_tt(k, q) != 0 && empiricalpwzy1(k, q, h, r) != 0) {
                        //if(empiricalpwzy(k, q, h, r)!=0) {
                        numentr2 += (double) -1 * empiricalpwzy1(k, q, h, r) * log(empiricalpwzy1(k, q, h, r));
                    }
                }
            }
        }
    }
    entropyhigh = 1 - (double) numentr2 / (denentr2);


    mat _Pwy(_K, _H);
    mat _Pwy1(_K, _H);
    mat _Pzy(_Q, _R);
    mat _Pzy1(_Q, _R);
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            _Pwy(k, h) = 0;
        }
    }
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _Pzy(q, r) = 0;
        }
    }
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            _Pwy1(k, h) = 0;
        }
    }
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _Pzy1(q, r) = 0;
        }
    }
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _Pwy(k, h) += empiricalpwzy(k, 0, h, r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int h = 0; h < _H; h++) {
                _Pzy(q, r) += empiricalpwzy(0, q, h, r);
            }
        }
    }

    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _Pwy1(k, h) += empiricalpwzy1(k, 0, h, r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int h = 0; h < _H; h++) {
                _Pzy1(q, r) += empiricalpwzy1(0, q, h, r);
            }
        }
    }

    centropy2 = 0;
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            if (_Pwy(k, h) != 0) {
                centropy2 += -_Pwy(k, h) * log(_Pwy(k, h));
            }
        }
    }

    double denentr = 0;
    for (int h = 0; h < _H; h++) {
        denentr += -phf[h] * log(phf[h]);
    }
    denentr=denentr*_K;
    centropy2 = 1 - centropy2 / denentr;

    centropy3 = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            if (_Pzy(q, r) != 0) {
                centropy3 += -_Pzy(q, r) * log(_Pzy(q, r));
            }
        }
    }

    denentr = 0;
    for (int r = 0; r < _R; r++) {
        denentr += -prf[r] * log(prf[r]);
    }
    denentr=denentr*_Q;
    centropy3 = 1 - centropy3 / denentr;


    centropy2 = 0;
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            if (_Pwy1(k, h) != 0) {
                centropy2 += -_Pwy1(k, h) * log(_Pwy1(k, h));
            }
        }
    }

    denentr = 0;
    for (int h = 0; h < _H; h++) {
        denentr += -phf[h] * log(phf[h]);
    }
    denentr=denentr*_K;

    centropy2 = 1 - centropy2 / denentr;

    centropy3 = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            if (_Pzy1(q, r) != 0) {
                centropy3 += -_Pzy1(q, r) * log(_Pzy1(q, r));
            }
        }
    }

    denentr = 0;
    for (int r = 0; r < _R; r++) {
        denentr += -prf[r] * log(prf[r]);
    }
    denentr=denentr*_Q;

    centropy3 = 1 - centropy3 / denentr;

    //BVR-------------------------------------------------------------------------------------
    arma::cube empiricalpwzyest(_n, _H, _R);

    /*for (int j = 0; j < _n; j++) {
        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                if (y(j, 1) == (k + 1) && y(j, 2) == (q + 1)) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            empiricalpwzyest(j, h, r) = empiricalpwzy(k, q, h, r);
                        }
                    }
                }
            }
        }
    }

    for (int i1 = 0; i1 < (_I + nlocdep); i1++) {
        for (int i2 = 0; i2 < (_I + nlocdep); i2++) {
            BVR(i1, i2) = 0;
        }
    }
    for (int i1 = 0; i1 < (_I + nlocdep); i1++) {
        for (int i2 = i1 + 1; i2 < (_I + nlocdep); i2++) {

            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    n(c, cc) = 0;
                }
            }
            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    for (int j = 0; j < _n; j++) {
                        if (y(j, i1 + 3) == c && y(j, i2 + 3) == cc) {
                            n(c, cc) += 1;
                        }
                    }
                }
            }

            vector<int> nr(C);
            vector<int> nc(C);
            for (int c = 0; c < C; c++) {
                nr[c] = 0;
            }
            for (int cc = 0; cc < C; cc++) {
                nc[cc] = 0;
            }
            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    nr[c] += n(c, cc);
                }
            }
            for (int cc = 0; cc < C; cc++) {
                for (int c = 0; c < C; c++) {
                    nc[cc] += n(c, cc);
                }
            }
            for (int c = 0; c < C; c++) {
                for (int cc = 0; cc < C; cc++) {
                    m0(c, cc) = 0;
                }
            }

            for (int c = 0; c < _vC[i1]; c++) {
                for (int cc = 0; cc < _vC[i2]; cc++) {
                    for (int j = 0; j < _n; j++) {
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                for (int l = 0; l < _L; l++) {
                                    m0(c, cc) += (double) empiricalpwzyest(j, h, r)* plf(l,h,r) * pif(i1, c, l, h, r) * pif(i2, cc, l, h, r);
                                }
                            }
                        }
                    }
                }
            }

            for (int c = 0; c < _vC[i1]; c++) {
                for (int cc = 0; cc < _vC[i2]; cc++) {
                    BVR(i1, i2) += (double) pow((n(c, cc) - m0(c, cc)), 2) / m0(c, cc);
                }
            }

            BVR(i1, i2) = (double) BVR(i1, i2) / ((_vC[i1] - 1) * (_vC[i2] - 1));

        }
    }


    for (int i = 0; i < (_I + nlocdep); i++) {
        BVRG[i] = 0;
    }

    for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

        for (int k = 0; k < _K; k++) {
            for (int cc = 0; cc < C; cc++) {
                n1(k, cc) = 0;
            }
        }
        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < C; cc++) {
                for (int j = 0; j < _n; j++) {
                    if (y(j, 1) - 1 == c && y(j, i1 + 3) == cc) {
                        n1(c, cc) += 1;
                    }
                }
            }
        }

        vector<int> nr1(_K);
        vector<int> nc1(C);
        for (int c = 0; c < _K; c++) {
            nr1[c] = 0;
        }
        for (int cc = 0; cc < C; cc++) {
            nc1[cc] = 0;
        }
        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < C; cc++) {
                nr1[c] += n1(c, cc);
            }
        }

        for (int cc = 0; cc < C; cc++) {
            for (int c = 0; c < _K; c++) {
                nc1[cc] += n1(c, cc);
            }
        }

        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < C; cc++) {
                m1(c, cc) = 0;
            }
        }
        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int l = 0; l < _L; l++) {
                            for (int q = 0; q < _Q; q++) {
                                m1(c, cc) += (double) empiricalpwzy(c, q, h, r) * plf(l, h, r) *
                                             pif(i1, cc, l, h, r);
                            }
                        }
                    }
                }
            }
        }

        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                m1(c,cc)=(double) m1(c,cc)*_nk[c];
            }
        }

        for (int c = 0; c < _K; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                BVRG[i1] += (double) pow((n1(c, cc) - m1(c, cc)), 2) / m1(c, cc);
            }
        }

        BVRG[i1] = (double) BVRG[i1] / ((_K - 1) * (_vC[i1] - 1));

    }

    for (int i = 0; i < (_I + nlocdep); i++) {
        BVRG1[i] = 0;
    }

    for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

        for (int q = 0; q < _Q; q++) {
            for (int cc = 0; cc < C; cc++) {
                n11(q, cc) = 0;
            }
        }
        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < C; cc++) {
                for (int j = 0; j < _n; j++) {
                    if (y(j, 2) - 1 == c && y(j, i1 + 3) == cc) {
                        n11(c, cc) += 1;
                    }
                }
            }
        }

        vector<int> nr11(_Q);
        vector<int> nc11(C);
        for (int c = 0; c < _Q; c++) {
            nr11[c] = 0;
        }
        for (int cc = 0; cc < C; cc++) {
            nc11[cc] = 0;
        }
        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < C; cc++) {
                nr11[c] += n11(c, cc);
            }
        }

        for (int cc = 0; cc < C; cc++) {
            for (int c = 0; c < _Q; c++) {
                nc11[cc] += n11(c, cc);
            }
        }

        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < C; cc++) {
                m11(c, cc) = 0;
            }
        }

        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int l = 0; l < _L; l++) {
                            for (int k = 0; k < _K; k++) {
                                m11(c, cc) += (double) empiricalpwzy(k, c, h, r) * plf(l, h, r) *
                                              pif(i1, cc, l, h, r);
                            }
                        }
                    }
                }
            }
        }

        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                m1(c,cc)=(double) m1(c,cc)*_nq[c];
            }
        }

        for (int c = 0; c < _Q; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                BVRG1[i1] += (double) pow((n11(c, cc) - m11(c, cc)), 2) / m11(c, cc);
            }
        }

        BVRG1[i1] = (double) BVRG1[i1] / ((_Q - 1) * (_vC[i1] - 1));

    }

    for (int i = 0; i < (_I + nlocdep); i++) {
        BVRP[i] = 0;
    }

    //vector<int> vetcoppie(_K);
    for (int k = 0; k < _K; k++) {
        vetcoppie[k] = _nk[k] * (_nk[k] - 1) / 2;
    }
    //int totcoppie = 0;
    for (int k = 0; k < _K; k++) {
        totcoppie += vetcoppie[k];
    }

    for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

        mat coppie(_K, _n);
        for (int k = 0; k < _K; k++) {
            for (int j = 0; j < _n; j++) {
                coppie(k, j) = -1;
            }
        }
        int co = 0;
        for (int k = 0; k < _K; k++) {
            for (int j = 0; j < _n; j++) {
                if (y(j, 1) - 1 == k) {
                    coppie(k, co) = y(j, i1 + 3);
                    co = co + 1;
                }
            }
            co = 0;
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                n2(c, cc) = 0;
            }
        }

        for (int k = 0; k < _K; k++) {

            for (int cc = 0; cc < _vC[i1]; cc++) {
                for (int c = 0; c < _vC[i1]; c++) {
                    for (int a = 0; a < _nk[k]; a++) {
                        if (coppie(k, a) == cc) {
                            for (int b = (a + 1); b < _nk[k]; b++) {
                                if (coppie(k, b) == c) {
                                    n2(cc, c) += 1;
                                }
                            }
                        }
                    }
                }
            }

        }

        vector<int> nr2(_vC[i1]);
        vector<int> nc2(_vC[i1]);
        for (int c = 0; c < _vC[i1]; c++) {
            nr2[c] = 0;
        }
        for (int cc = 0; cc < _vC[i1]; cc++) {
            nc2[cc] = 0;
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                nr2[c] += n2(c, cc);
            }
        }

        for (int cc = 0; cc < _vC[i1]; cc++) {
            for (int c = 0; c < _vC[i1]; c++) {
                nc2[cc] += n2(c, cc);
            }
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                m2(c, cc) = 0;
            }
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                for(int k=0;k<_K;k++) {

                    mat as(_H,_R);
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            as(h,r) = 0;
                        }
                    }
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int l = 0; l < _L; l++) {
                                as(h, r) += pif(i1, c, l, h, r) * plf(l, h, r);
                            }
                        }
                    }
                    mat as1(_H,_R);
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            as1(h,r) = 0;
                        }
                    }
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int l = 0; l < _L; l++) {
                                as1(h,r) += pif(i1, cc, l,h, r) * plf(l,h, r);
                            }
                        }
                    }

                    double as2=0;
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int q = 0; q < _Q; q++) {
                                as2 += as(h, r) * as1(h, r) * empiricalpwzy(k, q, h, r);
                            }
                        }
                    }

                    m2(c, cc) += (double) vetcoppie[k]*as2;

                }
            }
        }

        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = c; cc < _vC[i1]; cc++) {
                if (c == cc) {
                    BVRP[i1] += (double) pow((n2(c, cc) - m2(c, cc)), 2) / m2(c, cc);
                }
                if (c != cc) {
                    BVRP[i1] +=
                            (double) pow(((n2(c, cc) + n2(cc, c)) - (m2(c, cc) + m2(cc, c))), 2) /
                            (m2(c, cc) + m2(cc, c));
                }
            }
        }

        BVRP[i1] = (double) (BVRP[i1] * _K) / (_n * _vC[i1] * (_vC[i1] - 1) / 2);

    }

    for (int i = 0; i < (_I + nlocdep); i++) {
        BVRP1[i] = 0;
    }

    //vector<int> vetcoppie1(_Q);
    for (int q = 0; q < _Q; q++) {
        vetcoppie1[q] = _nq[q] * (_nq[q] - 1) / 2;
    }
    //int totcoppie1 = 0;
    for (int q = 0; q < _Q; q++) {
        totcoppie1 += vetcoppie1[q];
    }

    for (int i1 = 0; i1 < (_I + nlocdep); i1++) {

        mat coppie1(_Q, _n);
        for (int q = 0; q < _Q; q++) {
            for (int j = 0; j < _n; j++) {
                coppie1(q, j) = -1;
            }
        }
        int co = 0;
        for (int q = 0; q < _Q; q++) {
            for (int j = 0; j < _n; j++) {
                if (y(j, 2) - 1 == q) {
                    coppie1(q, co) = y(j, i1 + 3);
                    co = co + 1;
                }
            }
            co = 0;
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                n22(c, cc) = 0;
            }
        }

        for (int q = 0; q < _Q; q++) {

            for (int cc = 0; cc < _vC[i1]; cc++) {
                for (int c = 0; c < _vC[i1]; c++) {
                    for (int a = 0; a < _nq[q]; a++) {
                        if (coppie1(q, a) == cc) {
                            for (int b = (a + 1); b < _nq[q]; b++) {
                                if (coppie1(q, b) == c) {
                                    n22(cc, c) += 1;
                                }
                            }
                        }
                    }
                }
            }

        }

        vector<int> nr22(_vC[i1]);
        vector<int> nc22(_vC[i1]);
        for (int c = 0; c < _vC[i1]; c++) {
            nr22[c] = 0;
        }
        for (int cc = 0; cc < _vC[i1]; cc++) {
            nc22[cc] = 0;
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                nr22[c] += n22(c, cc);
            }
        }

        for (int cc = 0; cc < _vC[i1]; cc++) {
            for (int c = 0; c < _vC[i1]; c++) {
                nc22[cc] += n22(c, cc);
            }
        }

        for (int c = 0; c < C; c++) {
            for (int cc = 0; cc < C; cc++) {
                m22(c, cc) = 0;
            }
        }
        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = 0; cc < _vC[i1]; cc++) {
                for(int q=0;q<_Q;q++) {

                    mat as(_H,_R);
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            as(h,r) = 0;
                        }
                    }
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int l = 0; l < _L; l++) {
                                as(h, r) += pif(i1, c, l, h, r) * plf(l, h, r);
                            }
                        }
                    }
                    mat as1(_H,_R);
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            as1(h,r) = 0;
                        }
                    }
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int l = 0; l < _L; l++) {
                                as1(h,r) += pif(i1, cc, l,h, r) * plf(l,h, r);
                            }
                        }
                    }

                    double as2=0;
                    for(int h=0;h<_H;h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int k = 0; k < _K; k++) {
                                as2 += as(h, r) * as1(h, r) * empiricalpwzy(k, q, h, r);
                            }
                        }
                    }

                    m22(c, cc) += (double) vetcoppie1[q]*as2;

                }
            }
        }


        for (int c = 0; c < _vC[i1]; c++) {
            for (int cc = c; cc < _vC[i1]; cc++) {
                if (c == cc) {
                    BVRP1[i1] += (double) pow((n22(c, cc) - m22(c, cc)), 2) / m22(c, cc);
                }
                if (c != cc) {
                    BVRP1[i1] +=
                            (double) pow(((n22(c, cc) + n22(cc, c)) - (m22(c, cc) + m22(cc, c))), 2) /
                            (m22(c, cc) + m22(cc, c));
                }
            }
        }

        BVRP1[i1] = (double) (BVRP1[i1] * _Q) / (_n * _vC[i1] * (_vC[i1] - 1) / 2);

    }*/

//-------------------------------------------------------------------------------

    std::fill(_indh.begin(),_indh.end(),0);
    for(int k=0;k<_K;k++){
        int _max=_w_hat(k,0);
        for(int h=0;h<_H;h++){
            if(_w_hat(k,h)>_max){
                _max=_w_hat(k,h);
                _indh[k]=h;
            }
        }
    }


    std::fill(_indr.begin(),_indr.end(),0);
    for(int q=0;q<_Q;q++){
        int _max2=_z_hat(q,0);
        for(int r=0;r<_R;r++){
            if(_z_hat(q,r)>_max2){
                _max2=_z_hat(q,r);
                _indr[q]=r;
            }
        }
    }


    _wjk.zeros(_n,_H);

    for(int j=0;j<_n;j++){
        for(int h=0;h<_H;h++){
            if(_indh[y(j,1)-1]==h){
                _wjk(j,h)=1;
            }
        }
    }

    _zjq.zeros(_n,_R);

    for(int j=0;j<_n;j++){
        for(int r=0;r<_R;r++){
            if(_indr[y(j,2)-1]==r){
                _zjq(j,r)=1;
            }
        }
    }

    for(int j=0;j<_n;j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _B(j,h,r)=0;
            }
        }
    }

    for(int j=0;j<_n;j++){
        for(int h=0;h<_H;h++){
            for(int r=0;r<_R;r++){
                _B(j,h,r)=_wjk(j,h)*_zjq(j,r);
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    if (_B(j, h, r) == 1) {
                        _Prob3(j, l) = log(plf(l, h, r)) + log(_f(j, l,h,r));
                    }
                }
            }
        }
    }

    vector<double> _maxl(_n);
    for (int j = 0; j < _n; j++) {
        _maxl[j] = _Prob3(j, 0);
        for (int l = 0; l < _L; l++) {
            if (_Prob3(j, l) > _maxl[j]) {
                _maxl[j] = _Prob3(j, l);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3(j, l) = _Prob3(j, l) - _maxl[j];
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3(j, l) = exp(_Prob3(j, l));
        }
    }

    for (int j = 0; j < _n; j++) {
        _denom[j] = 0;
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _denom[j] += _Prob3(j, l);
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3(j, l) = (double) _Prob3(j, l) / _denom[j];
        }
    }

//statistics
    vector<int> _indL(_n);
    for (int j = 0; j < _n; j++) {
        double _maxl = _Prob3_hat(j, 0);
        for (int l = 0; l < _L; l++) {
            if (_Prob3_hat(j, l) > _maxl) {
                _maxl = _Prob3_hat(j, l);
                _indL[j] = l;
            }
        }
    }

    vector<int> _indH(_K);
    std::fill(_indH.begin(), _indH.end(), 0);
    for (int k = 0; k < _K; k++) {
        double _max2 = _Prob1_hat1(k, 0);
        for (int h = 0; h < _H; h++) {
            if (_Prob1_hat1(k, h) > _max2) {
                _max2 = _Prob1_hat1(k, h);
                _indH[k] = h;
            }
        }
    }

    vector<int> _indR(_Q);
    std::fill(_indR.begin(), _indR.end(), 0);
    for (int q = 0; q < _Q; q++) {
        double _max2 = _Prob2_hat1(q, 0);
        for (int r = 0; r < _R; r++) {
            if (_Prob2_hat1(q, r) > _max2) {
                _max2 = _Prob2_hat1(q, r);
                _indR[q] = r;
            }
        }
    }

    //classification error
    double ce = 0;
    double ce2 = 0;
    double ce3 = 0;
    for (int j = 0; j < _n; j++) {
        ce += 1 - _Prob3_hat(j, _indL[j]);
        //ce += 1 - _Prob3(j, _indL[j]);
    }
    for (int k = 0; k < _K; k++) {
        //ce2 += 1 - _Prob1_hat(k, _indH[k]);
        ce2 += 1 - _Prob1_hat1(k, _indH[k]);
        //ce2 += 1 - _Pwy1(k, _indH[k]);
    }
    for (int q = 0; q < _Q; q++) {
        //ce3 += 1 - _Prob2_hat(q, _indR[q]);
        ce3 += 1 - _Prob2_hat1(q, _indR[q]);
        //ce3 += 1 - _Pzy1(q, _indR[q]);
    }
    ce = ce / _n;
    ce2 = ce2 / _K;
    ce3 = ce3 / _Q;

    //lambda
    vector<double> _Plmean(_L);
    for (int l = 0; l < _L; l++) {
        _Plmean[l] = 0;
    }
    for (int l = 0; l < _L; l++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _Plmean[l] += plf(l, h, r) * phf[h] * prf[r];
            }
        }
    }

    double _maxx;
    for (int l = 0; l < _L; l++) {
        _maxx = _Plmean[0];
        if (_Plmean[l] > _maxx) {
            _maxx = _Plmean[l];
        }
    }

    double lambda, lambda2, lambda3;
    lambda = 1 - ce / (1 - _maxx);

    for (int h = 0; h < _H; h++) {
        _maxx = _p_h_hat[0];
        if (_p_h_hat[h] > _maxx) {
            _maxx = phf[h];
        }
    }

    lambda2 = 1 - ce2 / (1 - _maxx);
    for (int r = 0; r < _R; r++) {
        _maxx = _p_r_hat[0];
        if (_p_r_hat[r] > _maxx) {
            _maxx = prf[r];
        }
    }

    lambda3 = 1 - ce3 / (1 - _maxx);

    mat Pxy(_n,_L);
    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            Pxy(j,l)=_Prob3_hat(j,l);
        }
    }

    mat _Pxy(_n,_L);
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Pxy(j, l) = (double) _f(j, l,0,0) * _Plmean[l];
        }
    }

    vector<double> denpxy(_n);
    for (int j = 0; j < _n; j++) {
        denpxy[j] = 0;
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            denpxy[j] += _Pxy(j, l);
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Pxy(j, l) = (double) _Pxy(j, l) / denpxy[j];
        }
    }
    numentr1=0;
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            numentr1 += (double) -1 * _Pxy(j, l) * log(_Pxy(j, l));
        }
    }

    numentr1=0;
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            numentr1 += (double) -1 * Pxy(j, l) * log(Pxy(j, l));
        }
    }


    denentr1=0;
    for (int l = 0; l < _L; l++) {
        denentr1 += -1 * _Plmean[l] * log(_Plmean[l]);
    }
    denentr1 = denentr1 * _n;

    entropylow = 1 - (double) numentr1 / (denentr1);

    //double explogL = 0;
    Tensor<double, 4> _Pxwzy(_n, _L, _H, _R);
    arma::cube _Pwzyest(_n, _H, _R);

    for (int j = 0; j < _n; j++) {
        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                if (y(j, 1) == (k+1) && y(j, 2) == (q+1)) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _Pwzyest(j, h, r) = empiricalpwzy1(k, q, h, r);
                        }
                    }
                }
            }
        }
    }

    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _Pxwzy(j,l,h,r)=_Pwzyest(j,h,r)*Pxy(j,l);
                }
            }
        }
    }

    mat Pxy1(_n,_L);
    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            Pxy1(j,l)=0;
        }
    }
    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    Pxy1(j,l)+=_Pxwzy(j,l,h,r);
                }
            }
        }
    }
    numentr1=0;
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            numentr1 += (double) -1 * Pxy1(j, l) * log(Pxy1(j, l));
        }
    }

    for (int j = 0; j < _n; j++) {
        double _maxl = Pxy1(j, 0);
        for (int l = 0; l < _L; l++) {
            if (Pxy1(j, l) > _maxl) {
                _maxl = Pxy1(j, l);
                _indL[j] = l;
            }
        }
    }
    ce=0;
    for (int j = 0; j < _n; j++) {
        ce += 1 - Pxy1(j, _indL[j]);
    }
    ce=ce/_n;


    int _b=1;
    for(int n=0;n<_n;n++){
        if(y(n,0)==_b){
            _classif2(_b-1,0)=_b;
            for(int k=0;k<_K;k++){
                if(y(n,1)==k+1){
                    _classif2(_b-1,1)=_indh[k]+1;
                }
            }
            for(int q=0;q<_Q;q++){
                if(y(n,2)==q+1){
                    _classif2(_b-1,2)=_indr[q]+1;
                }
            }
            _b=_b+1;
            n=-1;
        }
    }

    mat _x_hat_ord(_n,_L);
    int _bb=1;
    for(int j=0;j<_n;j++){
        if(y(j,0)==_bb){
            for(int l=0;l<_L;l++){
                //_x_hat_ord(_bb-1,l)=_x_hat(j,l);
                _x_hat_ord(_bb-1,l)=_Prob3(j,l);
            }
            _bb=_bb+1;
            j=-1;
        }
    }

    std::fill(_indl.begin(),_indl.end(),0);
    for(int j=0;j<_n;j++){
        double _max3 = _x_hat_ord(j,0);
        for(int l=0;l<_L;l++){
            if(_x_hat_ord(j,l) > _max3){
                _max3 = _x_hat_ord(j,l);
                _indl[j]=l;
            }
        }
    }

 for(int j=0;j<_n;j++){
        _classif(j,0)=j+1;
        _classif(j,1)=_indl[j]+1;
        _classif(j,2)=_classif2(j,1);
        _classif(j,3)=_classif2(j,2);
    }

    cout<<"\n";

    double AIC3, AIC, BIC, SABIC, CAIC,AIC3e,AICe,BICg,CAICg,BICe,SABICe,SABICg,CAICe,BICge,SABICge,CAICge;
    int npar=0;
    for(int i=0;i<_I;i++) {
        npar += _L * (_vC[i]-1);
    }
    //ICL_BIC=-2*logLapprox2+((_R-1)*log(_Q))+((_H-1)*log(_K))+(_H*_R*(_L-1)*log(_n))+(npar*log(_n));
    //ICL_BICg=-2*logLapprox2+((_R-1)*log(_Q))+((_H-1)*log(_K))+(_H*_R*(_L-1)*log(_n))+(npar*log(ckq));
    //ICL_BIC2=-2*logLapprox2+((_R-1)+(_H-1)+(_H*_R*(_L-1))+npar)*log(_n);
    //ICL_BIC2g=-2*logLapprox2+((_R-1)+(_H-1)+(_H*_R*(_L-1))+npar)*log(ckq);
    AIC3=-2*logLapprox2+ 3 * (npar+_H+_R-2+(_H*_R*(_L-1)));
    AIC=-2*logLapprox2+ 2 * (npar+_H+_R-2+(_H*_R*(_L-1)));
    //ICL_BICe=-2*explogL2+((_R-1)*log(_Q))+((_H-1)*log(_K))+(_H*_R*(_L-1)*log(_n))+(npar*log(_n));
    //ICL_BICge=-2*explogL2+((_R-1)*log(_Q))+((_H-1)*log(_K))+(_H*_R*(_L-1)*log(_n))+(npar*log(ckq));
    //ICL_BIC2e=-2*explogL2+((_R-1)+(_H-1)+(_H*_R*(_L-1))+npar)*log(_n);
    //ICL_BIC2ge=-2*explogL2+((_R-1)+(_H-1)+(_H*_R*(_L-1))+npar)*log(ckq);
    AIC3e=-2*explogL2+ 3 * (npar+_H+_R-2+(_H*_R*(_L-1)));
    AICe=-2*explogL2+ 2 * (npar+_H+_R-2+(_H*_R*(_L-1)));
    BICe = -2 * explogL2 + log(_n) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    BICge = -2 * explogL2 + log(ckq) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    SABICe = -2 * explogL2 + log((_n+2)/24) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    SABICge = -2 * explogL2 + log((ckq+2)/24) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    CAICe = -2 * explogL2 + (log(_n)+1) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    CAICge = -2 * explogL2 + (log(ckq)+1) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    BIC = -2 * logLapprox2 + log(_n) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    BICg = -2 * logLapprox2 + log(ckq) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    SABIC = -2 * logLapprox2 + log((_n+2)/24) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    SABICg = -2 * logLapprox2 + log((ckq+2)/24) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    CAIC = -2 * logLapprox2 + (log(_n)+1) * (npar+_H+_R-2+(_H*_R*(_L-1)));
    CAICg = -2 * logLapprox2 + (log(ckq)+1) * (npar+_H+_R-2+(_H*_R*(_L-1)));

    auto end = chrono::steady_clock::now();
    int elapsed_time = chrono::duration_cast<chrono::seconds>(end - start).count();
    std::cout << "Time: " << elapsed_time << " sec" << std::endl;

}


