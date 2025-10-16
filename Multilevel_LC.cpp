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

    int _L,_H,_R,_K,_Q,_iter,_iter2,_iter3,_starting,_burn,_thinn;
    _L=4;
    _R=3;
    _iter=2000;
    _iter2=500;
    _iter3=500;
    _starting=5;
    _burn=150;
    _thinn=5;


    //thinning
    int nn=0;
    int nnn=_iter2-_burn;
    while(_thinn<=nnn){
        nn=nn+1;
        nnn=nnn-_thinn;
    }
    vector<int> it(nn);
    it[0]=_burn+_thinn-1;
    for(int j=1;j<nn;j++){
        it[j]=it[j-1]+_thinn;
    }
    int mu=0;


    mat y;
    y=vector_to_matrix(student_vector);

    _Q = y(0, 1);
    for (int j = 0; j < _n; j++) {
        if (y(j, 1) > _Q) {
            _Q = y(j, 1);
        }
    }
    cout<<_Q<<endl;

    vector<int> _vettoreqk(_n);
    int xx = 0;
    for (int q = 0; q < _Q; q++) {
        for (int j = 0; j < _n; j++) {
            if (y(j, 1) == q + 1) {
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

    vector<int> _nq(_Q);

    for (int q = 0; q < _Q; q++) {
        _nq[q] = 0;
    }
    for (int q = 0; q < _Q; q++) {
        for (int j = 0; j < _n; j++) {
            if (y(j, 1) == q + 1) {
                _nq[q] += 1;
            }
        }
    }

    mat _pi_l(_L, _R);
    for (int r = 0; r < _R; r++) {
        for (int l = 0; l < _L; l++) {
            _pi_l(l, r) = (double) 1 / _L;
        }
    }

int _I = _G - 2;
    vector<int> item(_I);
    item[0]=0;
    item[1]=0;
    item[2]=0;
    item[3]=0;
    item[4]=0;
    item[5]=0;

    vector<int> vC(_I);
    for(int i=0;i<_I;i++){
        vC[i]=2;  }
    int ma;
    for(int i=0;i<_I;i++){
        if(item[i]==0) {
            ma = y(0, i + 2);
            for (int j = 0; j < _n; j++) {
                if (y(j, i + 2) > ma) {
                    ma = y(j, i + 2);
                }
            }
            vC[i] = ma + 1;
        }
    }


    int _C;
    int maa;
    maa = vC[0];
    for (int i = 0; i < _I; i++) {
        if (vC[i] > maa) {
            maa = vC[i];
        }
    }
    _C=maa;

    mat _Prob2(_Q, _R);
    //boost::mt19937 generator2(1234);


    mat _z(_Q, _R);
    _z.zeros(_Q, _R);
    vector<double> vec(_R);
    double eqprob = (double) 1 / _R;
    std::fill(vec.begin(), vec.end(), eqprob);
    std::random_device rdtest;
    std::mt19937 generator(rdtest());
    boost::random::discrete_distribution<int> distribution2(vec.begin(), vec.end());
    for (int q = 0; q < _Q; q++) {
        int sample = distribution2(generator);
        _z(q, sample) = 1;
    }

    vector<int> _numr(_R);
    vector<double> _p_r(_R);

    std::fill(_numr.begin(), _numr.end(), 0);

    for (int r = 0; r < _R; r++) {
        for (int q = 0; q < _Q; q++) {
            _numr[r] += _z(q, r);
        }
    }


    for (int r = 0; r < _R; r++) {
        _p_r[r] = (double) _numr[r] / _Q;
    }

    mat _z_hat(_Q, _R);

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _z_hat(q, r) = 0;
        }
    }


    arma::cube _f(_n, _L, _R);
    mat _num(_n, _L);
    vector<double> _den(_n);
    vector<double> _denom(_n);
    mat _Prob3(_n, _L);
    mat _ff2(_Q, _R);
    mat _res2(_Q, _R);
    mat _num3(_Q, _R);
    vector<double> _den3(_Q);
    //mat _wjk(_n, _H);
    mat _zjq(_n, _R);
    vector<double> _den4(_R);
    //Tensor<double,4> _aa(_n,_L,_H,_R);
    mat _num4(_L, _R);
    vector<double> _den5(_L);
    mat _den6(_L,_R);
    Tensor<double,4> _num5(_I, _C, _L, _R);
    arma::cube _pf(_n, _L, _R);
    mat _fykq(_n,_R);
    arma::cube _res3(_n, _L,_R);
    mat _xjkq(_n, _L);
    mat _x_hat(_n, _L);
    _x_hat.zeros(_n, _L);
    mat _x_hat_int(_n, _L);
    _x_hat_int.zeros(_n, _L);
    double _logL1;
    arma::cube _Prob4(_n, _L, _R);
    arma::cube _Prob3i(_n, _L, _R);

    vector<double> _p_rc_hat(_R);
    vector<double> _p_rcc_hat(_R);
    std::fill(_p_rc_hat.begin(), _p_rc_hat.end(), 0);

     vector<int> _indr(_Q);
    mat _classif2(_n, 3);

    vector<int> _indl(_n);

    mat _classif(_n, 3);

    double _maxL;

    mat nxz(_L,_R);
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            nxz(l,r)=0;
        }
    }

    mat empiricalpzy1(_Q, _R);
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            empiricalpzy1(q, r) = 0;
        }
    }
    arma::cube empiricalpxzy1(_n,_L, _R);
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            empiricalpxzy1(j,l, r) = 0;
        }
    }
    }
    arma::cube empiricalpxzy2(_n,_L, _R);
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                empiricalpxzy2(j,l, r) = 0;
            }
        }
    }

    double LogL=0;

    mat _Prob3_hat(_n,_L);
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3_hat(j, l) = 0;
        }
    }

    arma::cube _Prob32(_n,_L,_R);
    arma::cube _Prob32_hat(_n,_L,_R);
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _Prob32_hat(j, l,r) = 0;
            }
        }
    }

    mat _Prob2_hat(_Q,_R);
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _Prob2_hat(q, r) = 0;
        }
    }

    mat _pi_l_hat(_L, _R);
    mat _pi_lc_hat(_L, _R);
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            _pi_l_hat(l,r) = 0;
        }
    }

    Tensor<double, 4> _p_hat(_I, _C, _L, _R);
    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _p_hat(i, c, l, r) = 0;
                }
            }
        }
    }

    vector<double> _p_r_hat(_R);
    std::fill(_p_r_hat.begin(), _p_r_hat.end(), 0);

    //arma::cube _p_h_all(_starting, _iter, _H);
    arma::cube _p_r_all(_starting, _iter, _R);
    //vector<vector<vector<double>>> _p_h_all;
    mat _logL_all(_starting, _iter);

    double explogL2=0;
    double Entropy2=0;
    double logLapprox2=0;

    //arma::cube _de(_n,_H,_R);
    //arma::cube _de2(_L,_H,_R);
    //Tensor<double,4> _tau(_n, _L,_H,_R);

    //arma_rng::set_seed_random();
    arma::cube _p1(_I, _C, _L);
    Tensor<double, 4> _p(_I, _C, _L, _R);
    //_p.setRandom();
    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _p(i, c, l, r) = 0;
                }
            }
        }
    }

    vector<double> _pl(_L);
    for (int l = 0; l < _L; l++) {
        _pl[l] = (double) 1 / _L;
    }

    for (int r = 0; r < _R; r++) {
        for (int l = 0; l < _L; l++) {
            _pi_l(l, r) = (double) 1 / _L;
        }
    }


    vector<double> _de(_n);
    vector<double> _de2(_L);
    mat _tau(_n, _L);

    double eps=0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000001;
    double eps2=0.000000000000000000000000000000000000000000000000000000000000000000000000000000000001;
    double eps3=0.000000000000000000000000000000000000000000000000000000000000000000000000000000000001;

    Tensor<double,4> _P(_starting,_I,_C,_L);
    Tensor<double,2> _PL(_starting,_L);
    vector<double> logLik(_starting);

    for (int st = 0; st < _starting; st++) {

        arma::cube _p1(_I, _C, _L);
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    _p1(i, c, l) = 0;
                }
            }
        }

        //---
        vector<int> _z0(_L);
        double eqprob0;
        vector<double> vec0(_L);
        eqprob0 = (double) 1 / (_L);
        std::fill(vec0.begin(), vec0.end(), eqprob0);
        boost::random::discrete_distribution<int> distribution0(vec0.begin(), vec0.end());
        for (int l = 0; l < _L; l++) {
            int sample = distribution0(generator);
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

        for (int i = 0; i < _I; i++) {

            if (item[i] == 0) {

                mat _pb(_C, _L);
                double eqprob;
                vector<double> vec(_C * _L);
                eqprob = (double) 1 / (_C * _L);
                std::fill(vec.begin(), vec.end(), eqprob);
                std::random_device rdtest;
                std::mt19937 generator(rdtest());
                boost::random::discrete_distribution<int> distribution2(vec.begin(), vec.end());
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        int sample = distribution2(generator);
                        _pb(c, l) = sample + 1;
                    }
                }

                vector<int> sumc(_L);
                for (int l = 0; l < _L; l++) {
                    sumc[l] = 0;
                }

                for (int l = 0; l < _L; l++) {
                    for (int c = 0; c < vC[i]; c++) {
                        sumc[l] += _pb(c, l);
                    }
                }

                for (int l = 0; l < _L; l++) {
                    for (int c = 0; c < vC[i]; c++) {
                        _p1(i, c, l) = (double) _pb(c, l) / sumc[l];
                    }
                }

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 2), l);
                    }
                }


            }

            if (item[i] == 1) {

                arma_rng::set_seed_random();
                mat _pu(_I, _L, fill::randu);
                arma::vec y1(_n);

                for (int j = 0; j < _n; j++) {
                    y1[j] = y(j, i + 2);
                }

                double maxy = y1.max();
                double miny = y1.min();

                for (int l = 0; l < _L; l++) {
                    _p1(i, 0, l) = _pu(i, l) * (maxy - miny) + miny;
                }

                for (int l = 0; l < _L; l++) {
                    _p1(i, 1, l) = var(y1);
                }

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _f1(j, l) = (double) _f1(j, l) * pow(2 * M_PI * _p1(i, 1, l), -0.5) *
                                    exp(-(0.5 * pow(_p1(i, 1, l), -1)) * pow((y(j, i + 2) - _p1(i, 0, l)), 2));
                    }
                }

            }

            if (item[i] == 2) {

                for (int l = 0; l < _L; l++) {
                    for (int j = 0; j < _n; j++) {
                        _p1(i, 0, l) += y(j, i + 2);
                    }
                }
                for (int l = 0; l < _L; l++) {
                    _p1(i, 0, l) = (double) _p1(i, 0, l) / _n;
                }

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _f1(j, l) = (double) _f1(j, l) * exp(-_p1(i, 0, l)) * pow(_p1(i, 0, l), (y(j, i + 2))) /
                                    fattoriale(y(j, i + 2));
                    }
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


        for (int u = 0; u < _iter; u++) {

            mat num0(_n,_L);
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) =(double) log(_f1(j,l))+log(_pl[l]);
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
                    _tau(j, l) = (double) num0(j,l) / _den0[j];
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
                        _p1(i,c, l) = 0;
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _f1(j, l) = 1;
                }
            }

            for (int i = 0; i < _I; i++) {

                if (item[i] == 0) {

                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i, y(j, i + 2), l) += _tau(j, l);
                        }
                    }

                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            _p1(i, c, l) = (double) _p1(i, c, l) / _de2[l];
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        for (int c = 0; c < _C; c++) {
                            if (_p1(i, c, l) > 0.9999) {
                                _p1(i, c, l) = 0.999;
                            }
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        for (int c = 0; c < _C; c++) {
                            if (_p1(i, c, l) < 0.0001) {
                                _p1(i, c, l) = 0.001;
                            }
                        }
                    }

                    for (int j = 0; j < _n; j++) {
                        for (int l = 0; l < _L; l++) {
                            _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 2), l);
                        }
                    }

                }

                if(item[i]==1){

                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i,0, l) += _tau(j, l) * y(j, i + 2);
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        _p1(i,0, l) = (double) _p1(i,0, l) / _de2[l];
                    }

                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i,1, l) += _tau(j, l) * pow((y(j, i + 2) - _p1(i,0, l)), 2);
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        _p1(i,1, l) = (double) _p1(i,1, l) / _de2[l];
                    }

                    for (int j = 0; j < _n; j++) {
                        for (int l = 0; l < _L; l++) {
                            _f1(j, l) = (double)_f1(j, l) * pow(2 * M_PI * _p1(i,1, l), -0.5) * exp(-(0.5 * pow(_p1(i,1, l), -1)) * pow((y(j, i + 2) - _p1(i,0, l)), 2));
                        }
                    }

                }

                if(item[i]==2){

                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i,0,l) += _tau(j, l) * y(j,i+2);
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        _p1(i,0,l) = (double) _p1(i,0,l) / _de2[l];
                    }

                    for (int j = 0; j < _n; j++) {
                        for (int l = 0; l < _L; l++) {
                            _f1(j, l) = (double)_f1(j, l) * exp(-_p1(i,0,l))*pow(_p1(i,0,l),(y(j,i+2)))/fattoriale(y(j,i+2));
                        }
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
            if (abs(L1 - L) < eps) {
                u = _iter - 1;
            }
            L1 = L;
            //cout<<u<<endl;
        }

        logLik[st] = L;
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    _P(st,i, c, l) = _p1(i, c, l);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            _PL(st,l)=_pl[l];
        }

    }

    double _maxlogl;
    int indicestarting=0;
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
                for (int r = 0; r < _R; r++) {
                    _p(i, c, l, r) = _P(indicestarting,i, c, l);
                }
            }
        }
    }

    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            _pi_l(l,r)=_PL(indicestarting,l);
        }
    }

    mat _z1(_Q, _R);
    _z1.zeros(_Q, _R);
    double eqprob2;
    eqprob2 = (double) 1 / _R;
    vector<double> vec3(_R);
    std::fill(vec3.begin(), vec3.end(), eqprob2);
    //std::random_device rdtest;
    std::mt19937 generator2(rdtest());
    boost::random::discrete_distribution<int> distribution3(vec3.begin(), vec3.end());
    for (int q = 0; q < _Q; q++) {
        int sample = distribution3(generator2);
        _z1(q, sample) = 1;
    }

    std::fill(_numr.begin(), _numr.end(), 0);

    for (int r = 0; r < _R; r++) {
        for (int q = 0; q < _Q; q++) {
            _numr[r] += _z1(q, r);
        }
    }

    //aggiornamento p_r

    for (int r = 0; r < _R; r++) {
        _p_r[r] = (double) _numr[r] / _Q;
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _f(j, l, r) = 1;
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int i = 0; i < _I; i++) {
                    _f(j, l, r) =(double) _f(j, l, r) * _p(i, y(j, i + 2), l, r);
                }
            }
        }
    }
    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pf(j, l, r) = (double) _pi_l(l, r) * _f(j, l, r);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            _fykq(j, r) = 0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _fykq(j, r) +=(double) _pf(j, l, r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _ff2(q, r) = 1;
        }
    }

    int mm = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int j = mm; j < _nq[q] + mm; j++) {
                _ff2(q, r) =(double) _ff2(q, r) * _fykq(_indvecqk[j], r);
            }
        }
        mm = mm + _nq[q];
    }

    double _logL = 0;
    vector<double> _tt(_Q);
    for (int q = 0; q < _Q; q++) {
        _tt[q] = 0;
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _tt[q] += _p_r[r] * _ff2(q, r);
        }
    }

    for (int q = 0; q < _Q; q++) {
        _logL += log(_tt[q]);
    }


    for (int u1 = 0; u1 < _iter3; u1++) {

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 0;
            }
        }

        mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int j = mm; j < _nq[q] + mm; j++) {
                    _ff2(q, r) +=(double) log(_fykq(_indvecqk[j], r));
                }
            }
            mm = mm + _nq[q];
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) =(double) log(_p_r[r]) + _ff2(q, r);
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

        for (int r = 0; r < _R; r++) {
            _numr[r] = 0;
        }
        for (int r = 0; r < _R; r++) {
            for (int q = 0; q < _Q; q++) {
                _numr[r] += _Prob2(q, r);
            }
        }

        for (int r = 0; r < _R; r++) {
            _p_r[r] = (double) _numr[r] / _Q;
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3i(j, l, r) = (double) log(_pi_l(l, r)) + log(_f(j, l, r));
                }
            }
        }

        mat _maxl(_n, _R);
        for (int r = 0; r < _R; r++) {
            for (int j = 0; j < _n; j++) {
                _maxl(j, r) = _Prob3i(j, 0, r);
                for (int l = 0; l < _L; l++) {
                    if (_Prob3i(j, l, r) > _maxl(j, r)) {
                        _maxl(j, r) = _Prob3i(j, l, r);
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3i(j, l, r) = _Prob3i(j, l, r) - _maxl(j, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3i(j, l, r) = exp(_Prob3i(j, l, r));
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            _denom[j] = 0;
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _denom[j] += _Prob3i(j, l, r);
                }
            }
        }

        mat _denomi(_n, _R);
        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _denomi(j, r) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _denomi(j, r) +=(double) _Prob3i(j, l, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3i(j, l, r) = (double) _Prob3i(j, l, r) / _denomi(j, r);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int j = 0; j < _n; j++) {
                    _Prob4(j, l, r) = _Prob2(y(j, 1) - 1, r) * _Prob3i(j, l, r);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _num4(l, r) = 0;
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int j = 0; j < _n; j++) {
                    _num4(l, r) += _Prob4(j, l, r);
                }
            }
        }

        for (int r = 0; r < _R; r++) {
            _den4[r] = 0;
        }

        for (int r = 0; r < _R; r++) {
            //for (int l = 0; l < _L; l++) {
            for (int j = 0; j < _n; j++) {
                _den4[r] += _Prob2(y(j, 1) - 1, r);
            }
            //}
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                if (_den4[r] != 0) {
                    _pi_l(l, r) = (double) _num4(l, r) / _den4[r];
                }
            }
        }

        int vvv = 0;
        vector<int> _vt3(_L);
        int _indv;
        for (int b = 0; b < _L; b++) {
            _vt3[b] = 0;
        }
        int v = 0;

        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                if (_pi_l(l, r) == 0) {
                    _pi_l(l, r) = 0.01;
                    _vt3[l] = 1;
                    v = v + 1;
                }
            }

            vector<int> _vt4(_L);
            for (int b = 0; b < _L; b++) {
                _vt4[b] = 0;
            }

            for (int u = 0; u < _L; u++) {
                if (_pi_l(u, r) <= ((v * 0.01) / (_L - v))) {
                    vvv = vvv + 1;
                    _vt4[u] = 1;
                }
            }

            if (vvv == 0 && v > 0) {
                for (int u = 0; u < _L; u++) {
                    if (_vt3[u] == 0) {
                        _pi_l(u, r) = _pi_l(u, r) - double((v * 0.01) / (_L -
                                                                         v));//se ho 0.04 0 0.96 uno verrà negativo e problema anche in p che 0 può essere vero 0
                    }
                }
            }

            if (vvv >= 1 && v > 0) {
                //for(int u=0;u<_L;u++){
                //if(_vt4[u]==0 &&_vt3[u]==0){
                //_pi_l(u,r)=_pi_l(u,r)-vvv*0.01;//se ho 0.04 0 0.96 uno verrà negativo e problema anche in p che 0 può essere vero 0
                //u=_L;  //ma non è detto che gli posso togliere 0.02
                //}
                //}
                int _maxx = _pi_l(0, r);
                for (int l = 0; l < _L; l++) {
                    if (_pi_l(l, r) > _maxx) {
                        _maxx = _pi_l(l, r);
                        _indv = l;
                    }
                }
                _pi_l(_indv, r) = _pi_l(_indv, r) - 0.01;
            }

            v = 0;
            vvv = 0;
            for (int b = 0; b < _L; b++) {
                _vt3[b] = 0;
            }
        }

        vector<int> _dephr(_I);
        for (int i = 0; i < _I; i++) {
            _dephr[i] = 0;
        }
        //_dephr[0]=2;
        //_dephr[1]=2;
        //_dephr[2]=2;
        //_dephr[3]=1;
        //_dephr[4]=2;
        //_dephr[5]=2;
        //_dephr[6]=2;
        //_dephr[7]=3;
        mat _num7(_n, _L);
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _num7(j, l) = 0;
            }
        }
        for (int l = 0; l < _L; l++) {
            for (int j = 0; j < _n; j++) {
                for (int r = 0; r < _R; r++) {
                    _num7(j, l) += _Prob4(j, l, r);
                }
            }
        }

        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _num5(i, c, l, r) = 0;
                    }
                }
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = 0; j < _n; j++) {
                            if (_dephr[i] == 1) {
                                if (y(j, i + 2) == c) {
                                    _num5(i, c, l, r) += _Prob4(j, l, r);
                                }
                            } else if (_dephr[i] == 0) {
                                if (y(j, i + 2) == c) {
                                    _num5(i, c, l, r) += _num7(j, l);
                                }
                            }
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
                for (int r = 0; r < _R; r++) {
                    _den5[l] += _Prob4(j, l, r);
                }
            }
        }

        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if (_dephr[i] == 1) {
                            if (_num4(l, r) != 0) {
                                _p(i, c, l, r) = (double) _num5(i, c, l, r) / _num4(l, r);
                            }
                        } else if (_dephr[i] == 0) {
                            if (_den5[l] != 0) {
                                _p(i, c, l, r) = (double) _num5(i, c, l, r) / _den5[l];
                            }
                        }
                    }
                }
            }
        }

        for (int r = 0; r < _R; r++) {
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int c = 0; c < vC[i]; c++) {
                        if (_p(i, c, l, r) == 1) {
                            _p(i, c, l, r) = 0.99;
                            if(c==0){
                                _p(i, 1, l, r) = 0.01;
                            }
                            if(c==1){
                                _p(i, 0, l, r) = 0.01;
                            }
                        }
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _f(j, l, r) = 1;
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    for (int i = 0; i < _I; i++) {
                        _f(j, l, r) =(double) _f(j, l, r) * _p(i, y(j, i + 2), l, r);
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _pf(j, l, r) = (double) _pi_l(l, r) * _f(j, l, r);
                }
            }
        }


        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _fykq(j, r) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _fykq(j, r) += _pf(j, l, r);
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 1;
            }
        }

        mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int j = mm; j < _nq[q] + mm; j++) {
                    _ff2(q, r) =(double) _ff2(q, r) * _fykq(_indvecqk[j], r);
                }
            }
            mm = mm + _nq[q];
        }

        _logL1 = 0;

        for (int q = 0; q < _Q; q++) {
            _tt[q] = 0;
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _tt[q] += _p_r[r] * _ff2(q, r);
            }
        }

        for (int q = 0; q < _Q; q++) {
            _logL1 += log(_tt[q]);
        }

        if (abs(_logL - _logL1) < eps3) {
            u1 = _iter3 - 1;
        }
        _logL = _logL1;
        //cout<<_logL<<endl;

    }

//--------------------------stochastic EM--------------------------------------

    _pi_l(0,0)=0.40;
    _pi_l(0,1)=0.15;
    _pi_l(0,2)=0.20;
    _pi_l(1,0)=0.30;
    _pi_l(1,1)=0.10;
    _pi_l(1,2)=0.40;
    _pi_l(2,0)=0.20;
    _pi_l(2,1)=0.35;
    _pi_l(2,2)=0.15;
    _pi_l(3,0)=0.10;
    _pi_l(3,1)=0.40;
    _pi_l(3,2)=0.25;
    _p_r[0]=0.5;
    _p_r[1]=0.3;
    _p_r[2]=0.2;
    _p(0,0,0,0)=0.8;
    _p(0,1,0,0)=0.2;
    _p(0,2,0,0)=0;
    _p(0,3,0,0)=0;
    _p(0,0,1,0)=0.9;
    _p(0,1,1,0)=0.1;
    _p(0,2,1,0)=0;
    _p(0,3,1,0)=0;
    _p(0,0,2,0)=0.3;
    _p(0,1,2,0)=0.7;
    _p(0,2,2,0)=0;
    _p(0,3,2,0)=0;
    _p(0,0,3,0)=0.1;
    _p(0,1,3,0)=0.9;
    _p(0,2,3,0)=0;
    _p(0,3,3,0)=0;
    _p(1,0,0,0)=0.7;
    _p(1,1,0,0)=0.3;
    _p(1,2,0,0)=0;
    _p(1,3,0,0)=0;
    _p(1,0,1,0)=0.8;
    _p(1,1,1,0)=0.2;
    _p(1,2,1,0)=0;
    _p(1,3,1,0)=0;
    _p(1,0,2,0)=0.2;
    _p(1,1,2,0)=0.8;
    _p(1,2,2,0)=0;
    _p(1,3,2,0)=0;
    _p(1,0,3,0)=0.3;
    _p(1,1,3,0)=0.7;
    _p(1,2,3,0)=0;
    _p(1,3,3,0)=0;
    _p(2,0,0,0)=0.7;
    _p(2,1,0,0)=0.1;
    _p(2,2,0,0)=0.2;
    _p(2,3,0,0)=0;
    _p(2,0,1,0)=0.8;
    _p(2,1,1,0)=0.05;
    _p(2,2,1,0)=0.15;
    _p(2,3,1,0)=0;
    _p(2,0,2,0)=0.1;
    _p(2,1,2,0)=0.7;
    _p(2,2,2,0)=0.2;
    _p(2,3,2,0)=0;
    _p(2,0,3,0)=0.15;
    _p(2,1,3,0)=0.65;
    _p(2,2,3,0)=0.2;
    _p(2,3,3,0)=0;
    _p(3,0,0,0)=0.8;
    _p(3,1,0,0)=0.05;
    _p(3,2,0,0)=0.15;
    _p(3,3,0,0)=0;
    _p(3,0,1,0)=0.1;
    _p(3,1,1,0)=0.7;
    _p(3,2,1,0)=0.2;
    _p(3,3,1,0)=0;
    _p(3,0,2,0)=0.7;
    _p(3,1,2,0)=0.1;
    _p(3,2,2,0)=0.2;
    _p(3,3,2,0)=0;
    _p(3,0,3,0)=0.2;
    _p(3,1,3,0)=0.7;
    _p(3,2,3,0)=0.1;
    _p(3,3,3,0)=0;
    _p(4,0,0,0)=0.65;
    _p(4,1,0,0)=0.1;
    _p(4,2,0,0)=0.1;
    _p(4,3,0,0)=0.15;
    _p(4,0,1,0)=0.1;
    _p(4,1,1,0)=0.65;
    _p(4,2,1,0)=0.1;
    _p(4,3,1,0)=0.15;
    _p(4,0,2,0)=0.6;
    _p(4,1,2,0)=0.1;
    _p(4,2,2,0)=0.2;
    _p(4,3,2,0)=0.1;
    _p(4,0,3,0)=0.2;
    _p(4,1,3,0)=0.6;
    _p(4,2,3,0)=0.1;
    _p(4,3,3,0)=0.1;
    _p(5,0,0,0)=0.75;
    _p(5,1,0,0)=0.05;
    _p(5,2,0,0)=0.05;
    _p(5,3,0,0)=0.15;
    _p(5,0,1,0)=0.2;
    _p(5,1,1,0)=0.6;
    _p(5,2,1,0)=0.05;
    _p(5,3,1,0)=0.15;
    _p(5,0,2,0)=0.7;
    _p(5,1,2,0)=0.1;
    _p(5,2,2,0)=0.1;
    _p(5,3,2,0)=0.1;
    _p(5,0,3,0)=0.1;
    _p(5,1,3,0)=0.7;
    _p(5,2,3,0)=0.15;
    _p(5,3,3,0)=0.05;

    for(int i=0;i<_I;i++) {
        for (int c = 0; c<_C; c++) {
            for (int l = 0; l<_L; l++) {
                for (int r = 0; r<_R; r++) {
                    _p(i, c, l, r) = _p(i, c, l, 0);
                }
            }
        }
    }

    double temp3;
    vector<int> temp4(_Q);
    int min2;


    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _f(j, l, r) = 1;
            }
        }
    }

    for (int i = 0; i < _I; i++) {

        if (item[i] == 0) {

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = _f(j, l,r) * _p(i, y(j, i + 2), l,r);
                    }
                }
            }

        }

        if(item[i]==1){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * pow(2 * M_PI * _p(i, 1, l,r), -0.5) *
                                     exp(-(0.5 * pow(_p(i, 1, l,r), -1)) * pow((y(j, i + 2) - _p(i, 0, l,r)), 2));
                    }
                }
            }

        }

        if(item[i]==2){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * exp(-_p(i, 0, l,r)) * pow(_p(i, 0, l,r), (y(j, i + 2))) /
                                     fattoriale(y(j, i + 2));
                    }
                }
            }

        }

    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pf(j, l, r) = (double) _pi_l(l, r) * _f(j, l, r);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            _fykq(j, r) = 0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _fykq(j, r) += _pf(j, l, r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _ff2(q, r) = 1;
        }
    }

    mm = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int j = mm; j < _nq[q] + mm; j++) {
                _ff2(q, r) = _ff2(q, r) * _fykq(_indvecqk[j], r);
            }
        }
        mm = mm + _nq[q];
    }

    _logL = 0;
    //vector<double> _tt(_Q);
    for (int q = 0; q < _Q; q++) {
        _tt[q] = 0;
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _tt[q] += _p_r[r] * _ff2(q, r);
        }
    }

    for (int q = 0; q < _Q; q++) {
        _logL += log(_tt[q]);
    }

    for(int u1=0;u1<_iter2;u1++) {

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 0;
            }
        }

        mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int j = mm; j < _nq[q] + mm; j++) {
                    _ff2(q, r) += log(_fykq(_indvecqk[j], r));
                }
            }
            mm = mm + _nq[q];
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

        int _pcv2=0;
        vector<int> _cv2(_R);
        while (_pcv2 == 0) {
            for (int r = 0; r < _R; r++) {
                _cv2[r] = 0;
            }
            _pcv2 = 1;
            _z.zeros(_Q,_R);
            for(int q=0;q<_Q;q++){
                std::random_device rdtest;
                std::mt19937 generator2(rdtest());
                boost::random::discrete_distribution<int> distribution (_Prob2.row(q).begin(),_Prob2.row(q).end());
                int sample = distribution(generator2);
                _z(q,sample) = 1;
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

        std::fill(_numr.begin(),_numr.end(),0);

        for(int r=0;r<_R;r++){
            for(int q=0;q<_Q;q++){
                _numr[r]+=_z(q,r);
            }
        }


        for(int r=0;r<_R;r++){
            _p_r[r]=(double)_numr[r]/_Q;
        }


        vector<double> _p_rc(_R);
        for (int r = 0; r < _R; r++) {
            _p_rc[r]=_p_r[r];
        }

        double temp6;
        vector<int> temp7(_Q);
        //mat temp8(_n,_H);
        //arma::cube temp9(_n,_L,_H);
        vector<double> temp10(_L);
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

        for(int j=0;j<_n;j++){
            for(int r=0;r<_R;r++){
                _zjq(j,r)=0;
            }
        }

        for(int j=0;j<_n;j++){
            for(int r=0;r<_R;r++){
                _zjq(j,r)=_z(y(j,1)-1,r);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j,l,r)= log(_pi_l(l,r))+log(_f(j,l,r));
                }
            }
        }

        mat _maxl2(_n, _R);
        for (int r = 0; r < _R; r++) {
            for (int j = 0; j < _n; j++) {
                _maxl2(j, r) = _Prob32(j, 0, r);
                for (int l = 0; l < _L; l++) {
                    if (_Prob32(j, l, r) > _maxl2(j, r)) {
                        _maxl2(j, r) = _Prob32(j, l, r);
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j, l, r) = _Prob32(j, l, r) - _maxl2(j, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j, l, r) = exp(_Prob32(j, l, r));
                }
            }
        }

        mat _denomi(_n, _R);
        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _denomi(j, r) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _denomi(j, r) +=(double) _Prob32(j, l, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j, l, r) = (double) _Prob32(j, l, r) / _denomi(j, r);
                }
            }
        }


        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    if(_zjq(j,r)==1){
                        _Prob3(j,l)= log(_pi_l(l,r))+log(_f(j,l,r));
                    }
                }
            }
        }

        vector<double>_maxl(_n);
        for(int j=0;j<_n;j++){
            _maxl[j]=_Prob3(j,0);
            for(int l=0;l<_L;l++){
                if(_Prob3(j,l)>_maxl[j]){
                    _maxl[j]=_Prob3(j,l);
                }
            }
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _Prob3(j,l)=_Prob3(j,l)-_maxl[j];
            }
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _Prob3(j,l)=exp(_Prob3(j,l));
            }
        }

        for(int j=0;j<_n;j++){
            _denom[j]=0;
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _denom[j]+=_Prob3(j,l);
            }
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _Prob3(j,l)=(double) _Prob3(j,l)/_denom[j];
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++){
                _xjkq(j,l)=0;
            }
        }

        for(int j=0;j<_n;j++){
            std::random_device rdtest;
            std::mt19937 generator3(rdtest());
            boost::random::discrete_distribution<int> distribution (_Prob3.row(j).begin(),_Prob3.row(j).end());
            int sample = distribution(generator3);
            _xjkq(j,sample) = 1;
        }

        for(int r=0;r<_R;r++){
            _den4[r]=0;
        }

        for(int r=0;r<_R;r++){
            for(int j=0;j<_n;j++){
                _den4[r]+=_zjq(j,r);
            }
        }

        for(int l=0;l<_L;l++){
            for(int r=0;r<_R;r++){
                _num4(l,r)=0;
            }
        }

        for(int l=0;l<_L;l++){
            for(int r=0;r<_R;r++){
                for(int j=0;j<_n;j++){
                    _num4(l,r)+=_zjq(j,r) * _xjkq(j,l);
                }
            }
        }

        for(int l=0;l<_L;l++){
            for(int r=0;r<_R;r++){
                if(_den4[r]!=0){
                    _pi_l(l,r)=(double) _num4(l,r) / _den4[r];
                }
            }
        }

        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                if (_pi_l(l, r) == 0) {
                    _pi_l(l, r) = 0.001;
                }
                if (_pi_l(l, r) == 1) {
                    _pi_l(l, r) = 0.999;
                }
            }
        }



        vector<int> _dephr(_I);
        for (int i = 0; i < _I; i++) {
            _dephr[i]=0;
        }
        //_dephr[0]=2;
        //_dephr[1]=2;
        //_dephr[2]=2;
        //_dephr[3]=1;
        //_dephr[4]=2;
        //_dephr[5]=2;
        //_dephr[6]=2;
        //_dephr[7]=3;

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
                    if (_zjq(j, r) == 1) {
                        _den6(l, r) += _xjkq(j, l);
                    }
                }
            }
        }

        for (int i = 0; i < _I; i++) {

            if (item[i] == 0) {

                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            _num5(i, c, l, r) = 0;
                        }
                    }
                }

                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            for (int j = 0; j < _n; j++) {
                                if (_dephr[i] == 1) {
                                    if (y(j, i + 2) == c && _xjkq(j, l) == 1 && _zjq(j, r) == 1) {
                                        _num5(i, c, l, r) += 1;
                                    }
                                } else if (_dephr[i] == 0) {
                                    if (y(j, i + 2) == c && _xjkq(j, l) == 1) {
                                        _num5(i, c, l, r) += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            if (_dephr[i] == 1) {
                                if (_den6(l, r) != 0) {
                                    _p(i, c, l, r) = (double) _num5(i, c, l, r) / _den6(l, r);
                                }
                            } else if (_dephr[i] == 0) {
                                if (_den5[l] != 0) {
                                    _p(i, c, l, r) = (double) _num5(i, c, l, r) / _den5[l];
                                }
                            }
                        }
                    }
                }


                for (int c = 0; c < vC[i]; c++) {
                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            if (_p(i, c, l, r) > 0.9999) {
                                _p(i, c, l, r) = 0.999;
                            }
                            if (_p(i, c, l, r) < 0.00001) {
                                _p(i, c, l, r) = 0.001;
                            }
                        }
                    }
                }

            }

            if (item[i] == 1) {

                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            _num5(i, 0, l, r) = 0;
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            for (int j = 0; j < _n; j++) {
                                if (_dephr[i] == 1) {
                                    if (_xjkq(j, l) == 1 && _zjq(j, r) == 1) {
                                        _num5(i, 0, l, r) += y(j, i + 2);
                                    }
                                } else if (_dephr[i] == 0) {
                                    if (_xjkq(j, l) == 1) {
                                        _num5(i, 0, l, r) += y(j, i + 2);
                                    }
                                }
                            }
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            if (_dephr[i] == 1) {
                                if (_den6(l, r) != 0) {
                                    _p(i, 0, l, r) = (double) _num5(i, 0, l, r) / _den6(l, r);
                                }
                            } else if (_dephr[i] == 0) {
                                if (_den5[l] != 0) {
                                    _p(i, 0, l, r) = (double) _num5(i, 0, l, r) / _den5[l];
                                }
                            }
                        }
                    }

                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _num5(i, 1, l, r) = 0;
                    }
                }

                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {// sample from dirichlet distribution
                        for (int j = 0; j < _n; j++) {
                            if (_dephr[i] == 1) {
                                if (_xjkq(j, l) == 1 && _zjq(j, r) == 1) {
                                    _num5(i, 1, l, r) += pow((y(j, i + 2) - _p(i,0, l,r)), 2);
                                }
                            } else if (_dephr[i] == 0) {
                                if (_xjkq(j, l) == 1) {
                                    _num5(i, 1, l, r) += pow((y(j, i + 2) - _p(i,0, l,r)), 2);
                                }
                            }
                        }
                    }
                }

                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if (_dephr[i] == 1) {
                            if (_den6(l, r) != 0) {
                                _p(i, 1, l, r) = (double) _num5(i, 1, l, r) / _den6(l, r);
                            }
                        } else if (_dephr[i] == 0) {
                            if (_den5[l] != 0) {
                                _p(i, 1, l, r) = (double) _num5(i, 1, l, r) / _den5[l];
                            }
                        }
                    }
                }


            }

            if (item[i] == 2) {

                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _num5(i, 0, l, r) = 0;
                    }
                }

                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = 0; j < _n; j++) {
                            if (_dephr[i] == 1) {
                                if (_xjkq(j, l) == 1 && _zjq(j, r) == 1) {
                                    _num5(i, 0, l, r) += y(j, i + 2);
                                }
                            } else if (_dephr[i] == 0) {
                                if (_xjkq(j, l) == 1) {
                                    _num5(i, 0, l, r) += y(j, i + 2);
                                }
                            }
                        }
                    }
                }

                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if (_dephr[i] == 1) {
                            if (_den6(l, r) != 0) {
                                _p(i, 0, l, r) = (double) _num5(i, 0, l, r) / _den6(l, r);
                            }
                        } else if (_dephr[i] == 0) {
                            if (_den5[l] != 0) {
                                _p(i, 0, l, r) = (double) _num5(i, 0, l, r) / _den5[l];
                            }
                        }
                    }
                }

            }

            }


for (int j = 0; j < _n; j++) {
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            _f(j, l, r) = 1;
        }
    }
}

        for (int i = 0; i < _I; i++) {

            if (item[i] == 0) {

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                            for (int r = 0; r < _R; r++) {
                                _f(j, l,r) = _f(j, l,r) * _p(i, y(j, i + 2), l,r);
                            }
                    }
                }

            }

            if(item[i]==1){

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                            for (int r = 0; r < _R; r++) {
                                _f(j, l,r) = (double) _f(j, l,r) * pow(2 * M_PI * _p(i, 1, l,r), -0.5) *
                                               exp(-(0.5 * pow(_p(i, 1, l,r), -1)) * pow((y(j, i + 2) - _p(i, 0, l,r)), 2));
                            }
                    }
                }

            }

            if(item[i]==2){

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                            for (int r = 0; r < _R; r++) {
                                _f(j, l,r) = (double) _f(j, l,r) * exp(-_p(i, 0, l,r)) * pow(_p(i, 0, l,r), (y(j, i + 2))) /
                                               fattoriale(y(j, i + 2));
                            }
                    }
                }

            }

        }


for (int j = 0; j < _n; j++) {
    for (int r = 0; r < _R; r++) {
        for (int l = 0; l < _L; l++) {
            _pf(j,l,r)= (double)_pi_l(l,r) * _f(j,l,r);
        }
    }
}


for (int j = 0; j < _n; j++) {
    for (int r = 0; r < _R; r++) {
        _fykq(j,r)=0;
    }
}

for (int j = 0; j < _n; j++) {
    for (int r = 0; r < _R; r++) {
        for (int l = 0; l < _L; l++) {
            _fykq(j,r) += _pf(j,l,r);
        }
    }
}

for (int q = 0; q < _Q; q++) {
    for (int r = 0; r < _R; r++) {
        _ff2(q, r) = 1;
    }
}

mm=0;
for (int q = 0; q < _Q; q++) {
    for (int r = 0; r < _R; r++) {
        for (int j = mm; j < _nq[q]+mm; j++) {
            _ff2(q,r) =_ff2(q,r)* _fykq(_indvecqk[j],r);
        }
    }
    mm=mm+_nq[q];
}


double _logL=0;
vector<double> _tt(_Q);
for(int q=0;q<_Q;q++){
    _tt[q]=0;
}

for(int q=0;q<_Q;q++){
    for(int r=0;r<_R;r++){
        _tt[q]+=_p_r[r]*_ff2(q,r);
    }
}

for(int q=0;q<_Q;q++){
    _logL+=log(_tt[q]);
}

//_logL_all(v,u1)=_logL;


//media parametri
if(u1>it[mu]){
    mu=mu+1;
}

if(u1>=_burn && u1!=it[mu]){

    LogL+=_logL;

    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            nxz(l,r)+=_den6(l,r);
        }
    }

    for(int r=0;r<_R;r++){
        _p_r_hat[r] +=_p_r[r];
    }

    for (int r = 0; r < _R; r++) {
        _p_rc_hat[r] += _p_rc[r];
    }


    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            _pi_l_hat(l,r) += _pi_l(l,r);
        }
    }


    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for(int l=0;l<_L;l++){
                for (int r = 0; r < _R; r++) {
                    _p_hat(i, c, l, r) += _p(i, c, l, r);
                }
            }
        }
    }

            double explogL = 0;
            Tensor<double, 3> _Pxzy(_n, _L, _R);
            mat _Pzyest(_n, _R);

            for (int j = 0; j < _n; j++) {
                for (int q = 0; q < _Q; q++) {
                    if (y(j, 1) == (q+1)) {
                        for (int r = 0; r < _R; r++) {
                            _Pzyest(j, r) = _Prob2(q, r);
                        }
                    }
                }
            }

            for(int j=0;j<_n;j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _Pxzy(j,l,r)=_Pzyest(j,r)*_Prob32(j,l,r);
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    explogL+=_Prob2(q,r)*log(_p_rc[r]);
                }
            }


            for(int j=0;j<_n;j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if(_pi_l(l, r) * _f(j, l, r)!=0) {
                            explogL += _Pxzy(j, l, r) * log(_pi_l(l, r) * _f(j, l, r));
                            //explogL += _Pxzy(j, l, r) * (log(_pi_l(l, r)) + log(_f(j, l, r)));
                        }
                    }
                }
            }

            explogL2+=explogL;

            double Entropy = 0;
            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    if(_Prob2(q,r)!=0) {
                        Entropy += -_Prob2(q, r) * log(_Prob2(q, r));
                    }
                }
            }

            double Entropy1=0;
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if (_Prob32(j, l,r) != 0) {
                            Entropy1 += -_Pxzy(j, l, r) * log(_Prob32(j, l,r));
                        }
                    }
                }
            }

            Entropy=Entropy+Entropy1;
            Entropy2+=Entropy;

    double logLapprox;
            logLapprox=explogL+Entropy;
            logLapprox2+=logLapprox;

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob3_hat(j, l) += (double) _Prob3(j, l);
                }
            }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _Prob32_hat(j, l, r) += (double) _Prob32(j, l, r);
            }
        }
    }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _Prob2_hat(q, r) += (double) _Prob2(q, r);
                }
            }


        }


        if(u1==_iter2-1){

            LogL=LogL/(_iter2-_burn-nn);

            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    nxz(l,r)=(double)nxz(l,r)/(_iter2-_burn-nn);
                }
            }

            //_p_r_hat
            for(int r=0;r<_R;r++){
                _p_r_hat[r]=(double)_p_r_hat[r]/(_iter2-_burn-nn);
            }

            for (int r = 0; r < _R; r++) {
                _p_rc_hat[r] = (double) _p_rc_hat[r] / (_iter2-_burn-nn);
            }
            for (int r = 0; r < _R; r++) {
                _p_rcc_hat[r] = (double) _p_rc_hat[r];
            }


            //_pi_l_hat
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _pi_l_hat(l,r) = (double)_pi_l_hat(l,r)/(_iter2-_burn-nn);
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _pi_lc_hat(l,r) = (double)_pi_l_hat(l,r);
                }
            }


            /*for (int i = 0; i < (_R - 1); i++) {
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
                    temp10[l] =(double) _pi_l_hat(l,i);
                    _pi_l_hat(l,i) = _pi_l_hat(l,min2);
                    _pi_l_hat(l,min2) = temp10[l];
                }

            }*/


            //_p_hat
            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for(int l=0;l<_L;l++){
                        for (int r = 0; r < _R; r++) {
                            _p_hat(i, c, l, r) = (double) _p_hat(i, c, l, r) / (_iter2 - _burn - nn);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob3_hat(j, l) = (double) _Prob3_hat(j, l) / (_iter2 - _burn - nn);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _Prob32_hat(j, l,r) = (double) _Prob32_hat(j, l,r) / (_iter2 - _burn - nn);
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _Prob2_hat(q, r) = (double) _Prob2_hat(q, r) / (_iter2 - _burn - nn);
                }
            }


            explogL2=(double)explogL2 / (_iter2 - _burn - nn);
            Entropy2=(double)Entropy2 / (_iter2 - _burn - nn);
            logLapprox2=(double)logLapprox2 / (_iter2 - _burn - nn);

        }

        for (int r = 0; r < _R; r++) {
            _p_r[r]=_p_rc[r];
        }

        //cout<<u1+2<<endl;
    }
    /*cout<<"LogL"<<endl;
    cout<<LogL<<endl;

    cout<<"explogL2"<<endl;
    cout<<explogL2<<endl;

    cout<<"Entropy2"<<endl;
    cout<<Entropy2<<endl;*/


    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _f(j, l, r) = 1;
            }
        }
    }

    for (int i = 0; i < _I; i++) {

        if (item[i] == 0) {

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = _f(j, l,r) * _p_hat(i, y(j, i + 2), l,r);
                    }
                }
            }

        }

        if(item[i]==1){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * pow(2 * M_PI * _p_hat(i, 1, l,r), -0.5) *
                                     exp(-(0.5 * pow(_p_hat(i, 1, l,r), -1)) * pow((y(j, i + 2) - _p_hat(i, 0, l,r)), 2));
                    }
                }
            }

        }

        if(item[i]==2){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * exp(-_p_hat(i, 0, l,r)) * pow(_p_hat(i, 0, l,r), (y(j, i + 2))) /
                                     fattoriale(y(j, i + 2));
                    }
                }
            }

        }

    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pf(j,l,r)= (double)_pi_lc_hat(l,r) * _f(j,l,r);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            _fykq(j,r)=0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _fykq(j,r) += _pf(j,l,r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _ff2(q, r) = 1;
        }
    }

    mm=0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int j = mm; j < _nq[q]+mm; j++) {
                _ff2(q,r) =_ff2(q,r)* _fykq(_indvecqk[j],r);
            }
        }
        mm=mm+_nq[q];
    }


    //double _logL=0;
    //vector<double> _tt(_Q);
    for(int q=0;q<_Q;q++){
        _tt[q]=0;
    }

    for(int q=0;q<_Q;q++){
        for(int r=0;r<_R;r++){
            _tt[q]+=_p_rc_hat[r]*_ff2(q,r);
        }
    }

    _logL=0;
    for(int q=0;q<_Q;q++){
        _logL+=log(_tt[q]);
    }

    double explogL = 0;
    Tensor<double, 3> _Pxzy(_n, _L, _R);
    mat _Pzyest(_n, _R);

    for (int j = 0; j < _n; j++) {
        for (int q = 0; q < _Q; q++) {
            if (y(j, 1) == (q+1)) {
                for (int r = 0; r < _R; r++) {
                    _Pzyest(j, r) = _Prob2_hat(q, r);
                }
            }
        }
    }

    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _Pxzy(j,l,r)=_Pzyest(j,r)*_Prob32_hat(j,l,r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            explogL+=_Prob2_hat(q,r)*log(_p_rc_hat[r]);
        }
    }

    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                if(_pi_l_hat(l, r) * _f(j, l, r)!=0) {
                    explogL += _Pxzy(j, l, r) * log(_pi_lc_hat(l, r) * _f(j, l, r));
                }
            }
        }
    }

    double Entropy = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            if(_Prob2_hat(q,r)!=0) {
                Entropy += -_Prob2_hat(q, r) * log(_Prob2_hat(q, r));
            }
        }
    }

    double Entropy1=0;
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                if (_Prob32_hat(j, l,r) != 0) {
                    Entropy1 += -_Pxzy(j, l, r) * log(_Prob32_hat(j, l,r));
                }
            }
        }
    }

    Entropy=Entropy+Entropy1;

    double logLapprox;
    logLapprox=explogL+Entropy;

    cout<<"LogL"<<endl;
    cout<<_logL<<endl;

    /*cout<<"explogL"<<endl;
    cout<<explogL<<endl;

    cout<<"Entropy"<<endl;
    cout<<Entropy<<endl;*/

    cout<<"logLapprox Stochastic EM (a)"<<endl;
    cout<<logLapprox<<endl;

    cout<<"logLapprox Stochastic EM (b)"<<endl;
    cout<<logLapprox2<<endl;


    //------------------------bayesian version-------------------------------

    _pi_l(0,0)=0.40;
    _pi_l(0,1)=0.15;
    _pi_l(0,2)=0.20;
    _pi_l(1,0)=0.30;
    _pi_l(1,1)=0.10;
    _pi_l(1,2)=0.40;
    _pi_l(2,0)=0.20;
    _pi_l(2,1)=0.35;
    _pi_l(2,2)=0.15;
    _pi_l(3,0)=0.10;
    _pi_l(3,1)=0.40;
    _pi_l(3,2)=0.25;
    _p_r[0]=0.5;
    _p_r[1]=0.3;
    _p_r[2]=0.2;
    _p(0,0,0,0)=0.8;
    _p(0,1,0,0)=0.2;
    _p(0,2,0,0)=0;
    _p(0,3,0,0)=0;
    _p(0,0,1,0)=0.9;
    _p(0,1,1,0)=0.1;
    _p(0,2,1,0)=0;
    _p(0,3,1,0)=0;
    _p(0,0,2,0)=0.3;
    _p(0,1,2,0)=0.7;
    _p(0,2,2,0)=0;
    _p(0,3,2,0)=0;
    _p(0,0,3,0)=0.1;
    _p(0,1,3,0)=0.9;
    _p(0,2,3,0)=0;
    _p(0,3,3,0)=0;
    _p(1,0,0,0)=0.7;
    _p(1,1,0,0)=0.3;
    _p(1,2,0,0)=0;
    _p(1,3,0,0)=0;
    _p(1,0,1,0)=0.8;
    _p(1,1,1,0)=0.2;
    _p(1,2,1,0)=0;
    _p(1,3,1,0)=0;
    _p(1,0,2,0)=0.2;
    _p(1,1,2,0)=0.8;
    _p(1,2,2,0)=0;
    _p(1,3,2,0)=0;
    _p(1,0,3,0)=0.3;
    _p(1,1,3,0)=0.7;
    _p(1,2,3,0)=0;
    _p(1,3,3,0)=0;
    _p(2,0,0,0)=0.7;
    _p(2,1,0,0)=0.1;
    _p(2,2,0,0)=0.2;
    _p(2,3,0,0)=0;
    _p(2,0,1,0)=0.8;
    _p(2,1,1,0)=0.05;
    _p(2,2,1,0)=0.15;
    _p(2,3,1,0)=0;
    _p(2,0,2,0)=0.1;
    _p(2,1,2,0)=0.7;
    _p(2,2,2,0)=0.2;
    _p(2,3,2,0)=0;
    _p(2,0,3,0)=0.15;
    _p(2,1,3,0)=0.65;
    _p(2,2,3,0)=0.2;
    _p(2,3,3,0)=0;
    _p(3,0,0,0)=0.8;
    _p(3,1,0,0)=0.05;
    _p(3,2,0,0)=0.15;
    _p(3,3,0,0)=0;
    _p(3,0,1,0)=0.1;
    _p(3,1,1,0)=0.7;
    _p(3,2,1,0)=0.2;
    _p(3,3,1,0)=0;
    _p(3,0,2,0)=0.7;
    _p(3,1,2,0)=0.1;
    _p(3,2,2,0)=0.2;
    _p(3,3,2,0)=0;
    _p(3,0,3,0)=0.2;
    _p(3,1,3,0)=0.7;
    _p(3,2,3,0)=0.1;
    _p(3,3,3,0)=0;
    _p(4,0,0,0)=0.65;
    _p(4,1,0,0)=0.1;
    _p(4,2,0,0)=0.1;
    _p(4,3,0,0)=0.15;
    _p(4,0,1,0)=0.1;
    _p(4,1,1,0)=0.65;
    _p(4,2,1,0)=0.1;
    _p(4,3,1,0)=0.15;
    _p(4,0,2,0)=0.6;
    _p(4,1,2,0)=0.1;
    _p(4,2,2,0)=0.2;
    _p(4,3,2,0)=0.1;
    _p(4,0,3,0)=0.2;
    _p(4,1,3,0)=0.6;
    _p(4,2,3,0)=0.1;
    _p(4,3,3,0)=0.1;
    _p(5,0,0,0)=0.75;
    _p(5,1,0,0)=0.05;
    _p(5,2,0,0)=0.05;
    _p(5,3,0,0)=0.15;
    _p(5,0,1,0)=0.2;
    _p(5,1,1,0)=0.6;
    _p(5,2,1,0)=0.05;
    _p(5,3,1,0)=0.15;
    _p(5,0,2,0)=0.7;
    _p(5,1,2,0)=0.1;
    _p(5,2,2,0)=0.1;
    _p(5,3,2,0)=0.1;
    _p(5,0,3,0)=0.1;
    _p(5,1,3,0)=0.7;
    _p(5,2,3,0)=0.15;
    _p(5,3,3,0)=0.05;

    for(int i=0;i<_I;i++) {
        for (int c = 0; c<_C; c++) {
            for (int l = 0; l<_L; l++) {
                for (int r = 0; r<_R; r++) {
                    _p(i, c, l, r) = _p(i, c, l, 0);
                }
            }
        }
    }

    arma::vec dir_param_pr(_R);
        for (int r = 0; r < _R; r++) {
            dir_param_pr(r) = _p_r[r];
        }

    arma::mat dir_param_plr(_L,_R);
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            dir_param_plr(l,r) = _pi_l(l,r);
        }
    }
    arma::cube dir_param_p(_I, _C, _L);
    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                dir_param_p(i, c, l) = _p(i,c,l,0);
            }
        }
    }

    explogL2=0;
    Entropy2=0;
    logLapprox2=0;
    LogL=0;
    mu=0;

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _f(j, l, r) = 1;
            }
        }
    }

    for (int i = 0; i < _I; i++) {

        if (item[i] == 0) {

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = _f(j, l,r) * _p(i, y(j, i + 2), l,r);
                    }
                }
            }

        }

        if(item[i]==1){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * pow(2 * M_PI * _p(i, 1, l,r), -0.5) *
                                     exp(-(0.5 * pow(_p(i, 1, l,r), -1)) * pow((y(j, i + 2) - _p(i, 0, l,r)), 2));
                    }
                }
            }

        }

        if(item[i]==2){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * exp(-_p(i, 0, l,r)) * pow(_p(i, 0, l,r), (y(j, i + 2))) /
                                     fattoriale(y(j, i + 2));
                    }
                }
            }

        }

    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pf(j, l, r) = (double) _pi_l(l, r) * _f(j, l, r);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            _fykq(j, r) = 0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _fykq(j, r) += _pf(j, l, r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _ff2(q, r) = 1;
        }
    }

    mm = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int j = mm; j < _nq[q] + mm; j++) {
                _ff2(q, r) = _ff2(q, r) * _fykq(_indvecqk[j], r);
            }
        }
        mm = mm + _nq[q];
    }

    _logL = 0;
    //vector<double> _tt(_Q);
    for (int q = 0; q < _Q; q++) {
        _tt[q] = 0;
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _tt[q] += _p_r[r] * _ff2(q, r);
        }
    }

    for (int q = 0; q < _Q; q++) {
        _logL += log(_tt[q]);
    }

    for(int u1=0;u1<_iter2;u1++) {

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 0;
            }
        }

        mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int j = mm; j < _nq[q] + mm; j++) {
                    _ff2(q, r) += log(_fykq(_indvecqk[j], r));
                }
            }
            mm = mm + _nq[q];
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

        int _pcv2=0;
        vector<int> _cv2(_R);
        while (_pcv2 == 0) {
            for (int r = 0; r < _R; r++) {
                _cv2[r] = 0;
            }
            _pcv2 = 1;
            _z.zeros(_Q,_R);
            for(int q=0;q<_Q;q++){
                std::random_device rdtest;
                std::mt19937 generator2(rdtest());
                boost::random::discrete_distribution<int> distribution (_Prob2.row(q).begin(),_Prob2.row(q).end());
                int sample = distribution(generator2);
                _z(q,sample) = 1;
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

        std::fill(_numr.begin(),_numr.end(),0);

        for(int r=0;r<_R;r++){
            for(int q=0;q<_Q;q++){
                _numr[r]+=_z(q,r);
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
            _p_rc[r]=_p_r[r];
        }

        double temp6;
        vector<int> temp7(_Q);
        //mat temp8(_n,_H);
        //arma::cube temp9(_n,_L,_H);
        vector<double> temp10(_L);
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

        for(int j=0;j<_n;j++){
            for(int r=0;r<_R;r++){
                _zjq(j,r)=0;
            }
        }

        for(int j=0;j<_n;j++){
            for(int r=0;r<_R;r++){
                _zjq(j,r)=_z(y(j,1)-1,r);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j,l,r)= log(_pi_l(l,r))+log(_f(j,l,r));
                }
            }
        }

        mat _maxl2(_n, _R);
        for (int r = 0; r < _R; r++) {
            for (int j = 0; j < _n; j++) {
                _maxl2(j, r) = _Prob32(j, 0, r);
                for (int l = 0; l < _L; l++) {
                    if (_Prob32(j, l, r) > _maxl2(j, r)) {
                        _maxl2(j, r) = _Prob32(j, l, r);
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j, l, r) = _Prob32(j, l, r) - _maxl2(j, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j, l, r) = exp(_Prob32(j, l, r));
                }
            }
        }

        mat _denomi(_n, _R);
        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _denomi(j, r) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _denomi(j, r) +=(double) _Prob32(j, l, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob32(j, l, r) = (double) _Prob32(j, l, r) / _denomi(j, r);
                }
            }
        }


        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    if(_zjq(j,r)==1){
                        _Prob3(j,l)= log(_pi_l(l,r))+log(_f(j,l,r));
                    }
                }
            }
        }

        vector<double>_maxl(_n);
        for(int j=0;j<_n;j++){
            _maxl[j]=_Prob3(j,0);
            for(int l=0;l<_L;l++){
                if(_Prob3(j,l)>_maxl[j]){
                    _maxl[j]=_Prob3(j,l);
                }
            }
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _Prob3(j,l)=_Prob3(j,l)-_maxl[j];
            }
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _Prob3(j,l)=exp(_Prob3(j,l));
            }
        }

        for(int j=0;j<_n;j++){
            _denom[j]=0;
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _denom[j]+=_Prob3(j,l);
            }
        }

        for(int j=0;j<_n;j++){
            for(int l=0;l<_L;l++){
                _Prob3(j,l)=(double) _Prob3(j,l)/_denom[j];
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++){
                _xjkq(j,l)=0;
            }
        }

        for(int j=0;j<_n;j++){
            std::random_device rdtest;
            std::mt19937 generator3(rdtest());
            boost::random::discrete_distribution<int> distribution (_Prob3.row(j).begin(),_Prob3.row(j).end());
            int sample = distribution(generator3);
            _xjkq(j,sample) = 1;
        }

        for(int r=0;r<_R;r++){
            _den4[r]=0;
        }

        for(int r=0;r<_R;r++){
            for(int j=0;j<_n;j++){
                _den4[r]+=_zjq(j,r);
            }
        }

        for(int l=0;l<_L;l++){
            for(int r=0;r<_R;r++){
                _num4(l,r)=0;
            }
        }

        for(int l=0;l<_L;l++){
            for(int r=0;r<_R;r++){
                for(int j=0;j<_n;j++){
                    _num4(l,r)+=_zjq(j,r) * _xjkq(j,l);
                }
            }
        }


        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                dir_param_plr(l, r) += _num4(l, r);
            }
        }

        arma::vec dir_param_pl(_L);
        for(int r=0;r<_R;r++) {
            for (int l1 = 0; l1 < _L; l1++) {
                dir_param_pl(l1) = dir_param_plr(l1, r);
            }
            arma::vec sample_dir_pl = rdirichlet(dir_param_pl);
            for (int l = 0; l < _L; l++) {
                _pi_l(l,r) = sample_dir_pl(l);
            }
        }


        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                if (_pi_l(l, r) == 0) {
                    _pi_l(l, r) = 0.001;
                }
                if (_pi_l(l, r) == 1) {
                    _pi_l(l, r) = 0.999;
                }
            }
        }


        vector<int> _dephr(_I);
        for (int i = 0; i < _I; i++) {
            _dephr[i]=0;
        }

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
                    _nicl(i,y(j, i + 2),l) += _xjkq(j, l);
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
                    for (int r = 0; r < _R; r++) {
                        _p(i, c, l, r) = sample_dir_pil(c);
                    }
                }
            }
        }

        //_dephr[0]=2;
        //_dephr[1]=2;
        //_dephr[2]=2;
        //_dephr[3]=1;
        //_dephr[4]=2;
        //_dephr[5]=2;
        //_dephr[6]=2;
        //_dephr[7]=3;

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
                    if (_zjq(j, r) == 1) {
                        _den6(l, r) += _xjkq(j, l);
                    }
                }
            }
        }


        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _f(j, l, r) = 1;
                }
            }
        }

        for (int i = 0; i < _I; i++) {

            if (item[i] == 0) {

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            _f(j, l,r) = _f(j, l,r) * _p(i, y(j, i + 2), l,r);
                        }
                    }
                }

            }

            if(item[i]==1){

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            _f(j, l,r) = (double) _f(j, l,r) * pow(2 * M_PI * _p(i, 1, l,r), -0.5) *
                                         exp(-(0.5 * pow(_p(i, 1, l,r), -1)) * pow((y(j, i + 2) - _p(i, 0, l,r)), 2));
                        }
                    }
                }

            }

            if(item[i]==2){

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        for (int r = 0; r < _R; r++) {
                            _f(j, l,r) = (double) _f(j, l,r) * exp(-_p(i, 0, l,r)) * pow(_p(i, 0, l,r), (y(j, i + 2))) /
                                         fattoriale(y(j, i + 2));
                        }
                    }
                }

            }

        }


        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _pf(j,l,r)= (double)_pi_l(l,r) * _f(j,l,r);
                }
            }
        }


        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _fykq(j,r)=0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _fykq(j,r) += _pf(j,l,r);
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 1;
            }
        }

        mm=0;
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int j = mm; j < _nq[q]+mm; j++) {
                    _ff2(q,r) =_ff2(q,r)* _fykq(_indvecqk[j],r);
                }
            }
            mm=mm+_nq[q];
        }


        double _logL=0;
        vector<double> _tt(_Q);
        for(int q=0;q<_Q;q++){
            _tt[q]=0;
        }

        for(int q=0;q<_Q;q++){
            for(int r=0;r<_R;r++){
                _tt[q]+=_p_r[r]*_ff2(q,r);
            }
        }

        for(int q=0;q<_Q;q++){
            _logL+=log(_tt[q]);
        }

//_logL_all(v,u1)=_logL;


        if(u1>it[mu]){
            mu=mu+1;
        }

        if(u1>=_burn && u1!=it[mu]){

            LogL+=_logL;

            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    nxz(l,r)+=_den6(l,r);
                }
            }

            for(int r=0;r<_R;r++){
                _p_r_hat[r] +=_p_r[r];
            }

            for (int r = 0; r < _R; r++) {
                _p_rc_hat[r] += _p_rc[r];
            }


            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _pi_l_hat(l,r) += _pi_l(l,r);
                }
            }


            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for(int l=0;l<_L;l++){
                        for (int r = 0; r < _R; r++) {
                            _p_hat(i, c, l, r) += _p(i, c, l, r);
                        }
                    }
                }
            }

            double explogL = 0;
            Tensor<double, 3> _Pxzy(_n, _L, _R);
            mat _Pzyest(_n, _R);

            for (int j = 0; j < _n; j++) {
                for (int q = 0; q < _Q; q++) {
                    if (y(j, 1) == (q+1)) {
                        for (int r = 0; r < _R; r++) {
                            _Pzyest(j, r) = _Prob2(q, r);
                        }
                    }
                }
            }

            for(int j=0;j<_n;j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _Pxzy(j,l,r)=_Pzyest(j,r)*_Prob32(j,l,r);
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    explogL+=_Prob2(q,r)*log(_p_rc[r]);
                }
            }


            for(int j=0;j<_n;j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if(_pi_l(l, r) * _f(j, l, r)!=0) {
                            explogL += _Pxzy(j, l, r) * log(_pi_l(l, r) * _f(j, l, r));
                            //explogL += _Pxzy(j, l, r) * (log(_pi_l(l, r)) + log(_f(j, l, r)));
                        }
                    }
                }
            }

            explogL2+=explogL;

            double Entropy = 0;
            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    if(_Prob2(q,r)!=0) {
                        Entropy += -_Prob2(q, r) * log(_Prob2(q, r));
                    }
                }
            }

            double Entropy1=0;
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if (_Prob32(j, l,r) != 0) {
                            Entropy1 += -_Pxzy(j, l, r) * log(_Prob32(j, l,r));
                        }
                    }
                }
            }

            Entropy=Entropy+Entropy1;
            Entropy2+=Entropy;

            double logLapprox;
            logLapprox=explogL+Entropy;
            logLapprox2+=logLapprox;

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob3_hat(j, l) += (double) _Prob3(j, l);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _Prob32_hat(j, l, r) += (double) _Prob32(j, l, r);
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _Prob2_hat(q, r) += (double) _Prob2(q, r);
                }
            }


        }


        if(u1==_iter2-1){

            LogL=LogL/(_iter2-_burn-nn);

            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    nxz(l,r)=(double)nxz(l,r)/(_iter2-_burn-nn);
                }
            }

            //_p_r_hat
            for(int r=0;r<_R;r++){
                _p_r_hat[r]=(double)_p_r_hat[r]/(_iter2-_burn-nn);
            }

            for (int r = 0; r < _R; r++) {
                _p_rc_hat[r] = (double) _p_rc_hat[r] / (_iter2-_burn-nn);
            }
            for (int r = 0; r < _R; r++) {
                _p_rcc_hat[r] = (double) _p_rc_hat[r];
            }


            //_pi_l_hat
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _pi_l_hat(l,r) = (double)_pi_l_hat(l,r)/(_iter2-_burn-nn);
                }
            }

            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _pi_lc_hat(l,r) = (double)_pi_l_hat(l,r);
                }
            }


            /*for (int i = 0; i < (_R - 1); i++) {
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
                    temp10[l] =(double) _pi_l_hat(l,i);
                    _pi_l_hat(l,i) = _pi_l_hat(l,min2);
                    _pi_l_hat(l,min2) = temp10[l];
                }

            }*/


            //_p_hat
            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for(int l=0;l<_L;l++){
                        for (int r = 0; r < _R; r++) {
                            _p_hat(i, c, l, r) = (double) _p_hat(i, c, l, r) / (_iter2 - _burn - nn);
                        }
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _Prob3_hat(j, l) = (double) _Prob3_hat(j, l) / (_iter2 - _burn - nn);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _Prob32_hat(j, l,r) = (double) _Prob32_hat(j, l,r) / (_iter2 - _burn - nn);
                    }
                }
            }

            for (int q = 0; q < _Q; q++) {
                for (int r = 0; r < _R; r++) {
                    _Prob2_hat(q, r) = (double) _Prob2_hat(q, r) / (_iter2 - _burn - nn);
                }
            }

            explogL2=(double)explogL2 / (_iter2 - _burn - nn);
            Entropy2=(double)Entropy2 / (_iter2 - _burn - nn);
            logLapprox2=(double)logLapprox2 / (_iter2 - _burn - nn);

        }

        for (int r = 0; r < _R; r++) {
            _p_r[r]=_p_rc[r];
        }

        //cout<<u1+2<<endl;
    }
    //four0.close();
    //foul0.close();
    //foul1.close();
    //foul2.close();

    //ofstream fout("C:/Users/nicol/Desktop/statistica/Progetto/filelogLmix.txt", ios::app);
    /*cout<<"LogL"<<endl;
    cout<<LogL<<endl;

    cout<<"explogL2"<<endl;
    cout<<explogL2<<endl;

    cout<<"Entropy2"<<endl;
    cout<<Entropy2<<endl;*/


    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _f(j, l, r) = 1;
            }
        }
    }

    for (int i = 0; i < _I; i++) {

        if (item[i] == 0) {

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = _f(j, l,r) * _p_hat(i, y(j, i + 2), l,r);
                    }
                }
            }

        }

        if(item[i]==1){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * pow(2 * M_PI * _p_hat(i, 1, l,r), -0.5) *
                                     exp(-(0.5 * pow(_p_hat(i, 1, l,r), -1)) * pow((y(j, i + 2) - _p_hat(i, 0, l,r)), 2));
                    }
                }
            }

        }

        if(item[i]==2){

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l,r) = (double) _f(j, l,r) * exp(-_p_hat(i, 0, l,r)) * pow(_p_hat(i, 0, l,r), (y(j, i + 2))) /
                                     fattoriale(y(j, i + 2));
                    }
                }
            }

        }

    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pf(j,l,r)= (double)_pi_lc_hat(l,r) * _f(j,l,r);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            _fykq(j,r)=0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _fykq(j,r) += _pf(j,l,r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _ff2(q, r) = 1;
        }
    }

    mm=0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int j = mm; j < _nq[q]+mm; j++) {
                _ff2(q,r) =_ff2(q,r)* _fykq(_indvecqk[j],r);
            }
        }
        mm=mm+_nq[q];
    }


    //double _logL=0;
    //vector<double> _tt(_Q);
    for(int q=0;q<_Q;q++){
        _tt[q]=0;
    }

    for(int q=0;q<_Q;q++){
        for(int r=0;r<_R;r++){
            _tt[q]+=_p_rc_hat[r]*_ff2(q,r);
        }
    }

    _logL=0;
    for(int q=0;q<_Q;q++){
        _logL+=log(_tt[q]);
    }

    explogL = 0;

    for (int j = 0; j < _n; j++) {
        for (int q = 0; q < _Q; q++) {
            if (y(j, 1) == (q+1)) {
                for (int r = 0; r < _R; r++) {
                    _Pzyest(j, r) = _Prob2_hat(q, r);
                }
            }
        }
    }

    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _Pxzy(j,l,r)=_Pzyest(j,r)*_Prob32_hat(j,l,r);
            }
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            explogL+=_Prob2_hat(q,r)*log(_p_rc_hat[r]);
        }
    }

    for(int j=0;j<_n;j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                if(_pi_l_hat(l, r) * _f(j, l, r)!=0) {
                    explogL += _Pxzy(j, l, r) * log(_pi_lc_hat(l, r) * _f(j, l, r));
                }
            }
        }
    }

    Entropy = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            if(_Prob2_hat(q,r)!=0) {
                Entropy += -_Prob2_hat(q, r) * log(_Prob2_hat(q, r));
            }
        }
    }

    Entropy1=0;
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                if (_Prob32_hat(j, l,r) != 0) {
                    Entropy1 += -_Pxzy(j, l, r) * log(_Prob32_hat(j, l,r));
                }
            }
        }
    }

    Entropy=Entropy+Entropy1;

    //double logLapprox;
    logLapprox=explogL+Entropy;

    /*cout<<"LogL"<<endl;
    cout<<_logL<<endl;

    cout<<"explogL"<<endl;
    cout<<explogL<<endl;

    cout<<"Entropy"<<endl;
    cout<<Entropy<<endl;*/

    cout<<"logLapprox Gibbs (a)"<<endl;
    cout<<logLapprox<<endl;

    cout<<"logLapprox Gibbs (b)"<<endl;
    cout<<logLapprox2<<endl;


    auto end = chrono::steady_clock::now();
    int elapsed_time = chrono::duration_cast<chrono::seconds>(end - start).count();
    std::cout << "Time: " << elapsed_time << " sec" << std::endl;

}
