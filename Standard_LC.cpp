// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include <boost/random.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <armadillo>
#include <unsupported/Eigen/CXX11/Tensor>
#define _USE_MATH_DEFINES
#include <chrono>

using namespace Eigen;
using namespace std;
using namespace arma;


//funzione per l'ordinamento dei dati rispetto a k e q.
bool comparator(std::vector<int> &a, std::vector<int> &b){
    if(a[1] < b[1]) return true;
    else if(a[1] == b[1])
    {
        if(a[2] < b[2]) return true;
        else return false;
    }
    else return false;
}

//lettura file
void parser(std::string filename, std::vector< std::vector<int> > &student_vector)
{
    std::fstream f(filename);

    if(!f.is_open())
    {
        std::cerr << "ERROR in opening file" << std::endl;
        exit(-1);
    }

    std::string line;

    while(std::getline(f, line))
    {
        std::vector<int> s(7);
        //sscanf(line.data(), "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d", &s[0], &s[1], &s[2], &s[3], &s[4], &s[5], &s[6], &s[7], &s[8], &s[9], &s[10]);
        sscanf(line.data(), "%d, %d, %d, %d, %d, %d, %d", &s[0], &s[1], &s[2], &s[3], &s[4], &s[5], &s[6]);
        student_vector.push_back(s);
    }

    f.close();
}

//funzione per convertire i dati letti come vettori in una matrice.
mat vector_to_matrix(const std::vector< std::vector<int> > &student_vector)
{
    mat y(student_vector.size(),student_vector[0].size());

    for(int i = 0; i < student_vector.size(); i++)
    {
        for(int j = 0; j < student_vector[i].size(); j++)
        {
            y(i,j) = student_vector[i][j];
        }
    }

    return y;
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

// [[Rcpp::export]]
int main() {

    //remove("C:/Users/nicol/Desktop/statistica/Progetto/fileBIC.txt");
    auto start = chrono::steady_clock::now();
    std::vector<std::vector<int> > student_vector;

    parser("C:/Users/nicol/Desktop/statistica/Progetto/dati_logLapprox_standard.txt", student_vector);
    std::sort(student_vector.begin(), student_vector.end(), comparator);

    int _n = student_vector.size();
    int _G;
    for (int j = 0; j < _n; j++) {
        _G = student_vector[j].size();
    }

    mat y;
    y = vector_to_matrix(student_vector);

    int _L, _iter2, _iter1, _vl;
    //_L=4;
    _iter2 =2000;
    _vl = 2;
    _iter1=2000;
    int burn = 500;

    int _iter0 = 1000;
    double eps0 = 0.000000000000001;
    double eps2 = 0.000000001;
    //vector<int> Seed(_starting);

    int _starting0 = 16;
    int _starting = 1;
    vector<double> logLik(_starting0);

    //double L, L1;
    vector<double> _flog(_n);

    int _I = _G - 1;
    vector<int> vC(_I);
    int ma;
    for(int i=0;i<_I;i++){
        ma = y(0,i+1);
        for (int j = 0; j < _n; j++) {
            if (y(j,i+1) > ma) {
                ma = y(j,i+1);
            }
        }
        vC[i]=ma+1;
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

    cout<< _C <<endl;

//different starting points
    for (int v = 1; v < _vl; v++) {
        _L = v + 1;
        _L = 4;

        //mat _f1(_n, _L);
        Tensor<double, 4> _P(_starting0, _I, _C, _L);
        mat _PL(_starting0, _L);
        arma::cube _Tau(_starting0, _n, _L);
        arma::cube _F(_starting0, _n, _L);
        mat _xj(_n, _L);

        //arma_rng::set_seed_random();
        //arma_rng::set_seed(122984);
        //arma::cube _p1(_I, _C, _L, fill::randu);
        arma::cube _p1(_I, _C, _L);
        //_p.setRandom();

        vector<double> _pl(_L);

       vector<double> _de(_n);
        vector<double> _de2(_L);
        mat _tau(_n, _L);


        for (int st = 0; st < _starting0; st++) {

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
                        _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 1), l);
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

                for (int i = 0; i < _I; i++) {
                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i, y(j, i + 1), l) += _tau(j, l);
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
                            if (_p1(i, c, l) < 0.0001) {
                                _p1(i, c, l) = 0.001;
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
                            _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 1), l);
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
            }


            logLik[st] = L;
            //cout<<logLik[st]<<endl;
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

            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        _P(st, i, c, l) = _p1(i, c, l);
                    }
                }
            }

            for (int l = 0; l < _L; l++) {
                _PL(st, l) = _pl[l];
            }


        }


        double _maxlogl;
        int indicestarting;
        _maxlogl = logLik[0];
        for (int st = 0; st < _starting0; st++) {
            if (logLik[st] > _maxlogl) {
                _maxlogl = logLik[st];
                indicestarting = st;
            }
        }


        double Entropy;
        Entropy = 0;
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                Entropy += (double) -_Tau(indicestarting, j, l) * log(_Tau(indicestarting, j, l));
            }
        }

        double completelogL;
        completelogL = 0;
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                completelogL += (double) _Tau(indicestarting, j, l) *
                                (log(_PL(indicestarting, l)) + log(_F(indicestarting, j, l)));
                //completelogL+=(double) _Tau(indicestarting,j,l)*log(_PL(indicestarting,l)*_F(indicestarting,j,l));
            }
        }
        //cout << "completelogL" << endl;
        //cout << completelogL << endl;

        //double CL;

        //int npar1;
        //double BIC,AIC,AIC3,CAIC,SABIC,CLC,ICL_BIC,AWE,ICL_AIC,ICL_AIC3,ICL_SABIC,ICL_CAIC;
        double L = logLik[indicestarting];
        //CL = L - Entropy;

        //cout << "logL" << endl;
        cout<<"logL"<<endl;
        cout << L << endl;

        //cout << "E" << endl;
        //cout << Entropy << endl;
        //cout << "Stocastico" << endl;


//----------------------stochastic EM------------------------
        vector<double> Entropy2(_starting);
        vector<double> completelogL2(_starting);
        vector<double> logL2(_starting);
        vector<double> logL3(_starting);
        for (int st1 = 0; st1 < _starting; st1++) {

            _pl[0]=0.4;
            _pl[1]=0.3;
            _pl[2]=0.2;
            _pl[3]=0.1;
            _p1(0,0,0)=0.8;
            _p1(0,1,0)=0.2;
            _p1(0,2,0)=0;
            _p1(0,3,0)=0;
            _p1(0,0,1)=0.9;
            _p1(0,1,1)=0.1;
            _p1(0,2,1)=0;
            _p1(0,3,1)=0;
            _p1(0,0,2)=0.3;
            _p1(0,1,2)=0.7;
            _p1(0,2,2)=0;
            _p1(0,3,2)=0;
            _p1(0,0,3)=0.1;
            _p1(0,1,3)=0.9;
            _p1(0,2,3)=0;
            _p1(0,3,3)=0;
            _p1(1,0,0)=0.7;
            _p1(1,1,0)=0.3;
            _p1(1,2,0)=0;
            _p1(1,3,0)=0;
            _p1(1,0,1)=0.8;
            _p1(1,1,1)=0.2;
            _p1(1,2,1)=0;
            _p1(1,3,1)=0;
            _p1(1,0,2)=0.2;
            _p1(1,1,2)=0.8;
            _p1(1,2,2)=0;
            _p1(1,3,2)=0;
            _p1(1,0,3)=0.3;
            _p1(1,1,3)=0.7;
            _p1(1,2,3)=0;
            _p1(1,3,3)=0;
            _p1(2,0,0)=0.7;
            _p1(2,1,0)=0.1;
            _p1(2,2,0)=0.2;
            _p1(2,3,0)=0;
            _p1(2,0,1)=0.8;
            _p1(2,1,1)=0.05;
            _p1(2,2,1)=0.15;
            _p1(2,3,1)=0;
            _p1(2,0,2)=0.1;
            _p1(2,1,2)=0.7;
            _p1(2,2,2)=0.2;
            _p1(2,3,2)=0;
            _p1(2,0,3)=0.15;
            _p1(2,1,3)=0.65;
            _p1(2,2,3)=0.2;
            _p1(2,3,3)=0;
            _p1(3,0,0)=0.8;
            _p1(3,1,0)=0.05;
            _p1(3,2,0)=0.15;
            _p1(3,3,0)=0;
            _p1(3,0,1)=0.1;
            _p1(3,1,1)=0.7;
            _p1(3,2,1)=0.2;
            _p1(3,3,1)=0;
            _p1(3,0,2)=0.7;
            _p1(3,1,2)=0.1;
            _p1(3,2,2)=0.2;
            _p1(3,3,2)=0;
            _p1(3,0,3)=0.2;
            _p1(3,1,3)=0.7;
            _p1(3,2,3)=0.1;
            _p1(3,3,3)=0;
            _p1(4,0,0)=0.65;
            _p1(4,1,0)=0.1;
            _p1(4,2,0)=0.1;
            _p1(4,3,0)=0.15;
            _p1(4,0,1)=0.1;
            _p1(4,1,1)=0.65;
            _p1(4,2,1)=0.1;
            _p1(4,3,1)=0.15;
            _p1(4,0,2)=0.6;
            _p1(4,1,2)=0.1;
            _p1(4,2,2)=0.2;
            _p1(4,3,2)=0.1;
            _p1(4,0,3)=0.2;
            _p1(4,1,3)=0.6;
            _p1(4,2,3)=0.1;
            _p1(4,3,3)=0.1;
            _p1(5,0,0)=0.75;
            _p1(5,1,0)=0.05;
            _p1(5,2,0)=0.05;
            _p1(5,3,0)=0.15;
            _p1(5,0,1)=0.2;
            _p1(5,1,1)=0.6;
            _p1(5,2,1)=0.05;
            _p1(5,3,1)=0.15;
            _p1(5,0,2)=0.7;
            _p1(5,1,2)=0.1;
            _p1(5,2,2)=0.1;
            _p1(5,3,2)=0.1;
            _p1(5,0,3)=0.1;
            _p1(5,1,3)=0.7;
            _p1(5,2,3)=0.15;
            _p1(5,3,3)=0.05;

            mat _f1(_n, _L);
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _f1(j, l) = 1;
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int i = 0; i < _I; i++) {
                        _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 1), l);
                    }
                }
            }

            vector<double> Entropy1(_iter1);
            vector<double> completelogL1(_iter1);

            mat plhat(_starting, _L);
            Tensor<double, 4> phat(_starting, _I, _C, _L);
            for (int st = 0; st < _starting; st++) {
                for (int l = 0; l < _L; l++) {
                    plhat(st, l) = 0;
                }
            }
            for (int st = 0; st < _starting; st++) {
                for (int i = 0; i < _I; i++) {
                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            phat(st, i, c, l) = 0;
                        }
                    }
                }
            }

            for (int u = 0; u < _iter1; u++) {

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

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _xj(j, l) = 0;
                    }
                }

                for (int j = 0; j < _n; j++) {
                    std::random_device rdtest1;
                    std::mt19937 generator1(rdtest1());
                    boost::random::discrete_distribution<int> distribution1(_tau.row(j).begin(), _tau.row(j).end());
                    int sample = distribution1(generator1);
                    _xj(j, sample) = 1;
                }

                for (int l = 0; l < _L; l++) {
                    _pl[l] = 0;
                }

                for (int l = 0; l < _L; l++) {
                    for (int j = 0; j < _n; j++) {
                        _pl[l] += _xj(j, l);
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
                        _de2[l] += _xj(j, l);
                    }
                }

                for (int i = 0; i < _I; i++) {
                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            _p1(i, c, l) = 0;
                        }
                    }
                }

                for (int i = 0; i < _I; i++) {
                    for (int l = 0; l < _L; l++) {
                        for (int j = 0; j < _n; j++) {
                            _p1(i, y(j, i + 1), l) += _xj(j, l);
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
                            if (_p1(i, c, l) < 0.0001) {
                                _p1(i, c, l) = 0.001;
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
                            _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 1), l);
                        }
                    }
                }

                if (u > burn) {

                    Entropy1[u] = 0;
                    for (int j = 0; j < _n; j++) {
                        for (int l = 0; l < _L; l++) {
                            Entropy1[u] += (double) -_tau(j, l) * log(_tau(j, l));
                        }
                    }

                    completelogL1[u] = 0;
                    for (int j = 0; j < _n; j++) {
                        for (int l = 0; l < _L; l++) {
                            completelogL1[u] += (double) _tau(j, l) * (log(_pl[l]) + log(_f1(j, l)));
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        plhat(st1, l) += _pl[l];
                    }

                    for (int i = 0; i < _I; i++) {
                        for (int c = 0; c < vC[i]; c++) {
                            for (int l = 0; l < _L; l++) {
                                phat(st1, i, c, l) += _p1(i, c, l);
                            }
                        }
                    }

                }

            }


            Entropy2[st1] = 0;
            completelogL2[st1] = 0;
            for (int u = 0; u < _iter1; u++) {
                Entropy2[st1] += Entropy1[u];
                completelogL2[st1] += completelogL1[u];
            }
            Entropy2[st1] = Entropy2[st1] / (_iter1 - burn);
            completelogL2[st1] = completelogL2[st1] / (_iter1 - burn);

            logL2[st1] = Entropy2[st1] + completelogL2[st1];

            for (int l = 0; l < _L; l++) {
                plhat(st1, l) = plhat(st1, l) / (_iter1 - burn);
            }

            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < vC[i]; c++) {
                    for (int l = 0; l < _L; l++) {
                        phat(st1, i, c, l) = phat(st1, i, c, l) / (_iter1 - burn);
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
                        _f1(j, l) =(double) _f1(j, l) * phat(st1, i, y(j, i + 1), l);
                    }
                }
            }
            mat num0(_n, _L);
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = (double) log(_f1(j, l)) + log(plhat(st1, l));
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
                    _Tau(st1, j, l) = (double) num0(j, l) / _den0[j];
                }
            }

            double Entropy3 = 0;
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    Entropy3 += (double) -_Tau(st1, j, l) * log(_Tau(st1, j, l));
                }
            }

            double completelogL3 = 0;
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    completelogL3 += (double) _Tau(st1, j, l) * (log(plhat(st1, l)) + log(_f1(j, l)));
                }
            }

            logL3[st1] = Entropy3 + completelogL3;

            cout<<"logLapprox stochastic EM (a)"<<endl;
            cout<<logL3[st1]<<endl;
            cout<<"logLapprox stochastic EM (b)"<<endl;
            cout<<logL2[st1]<<endl;

        }


 //----------------------------bayesian version--------------------------------

        //vector<double> Entropy2(_starting);
        //vector<double> completelogL2(_starting);
        //vector<double> logL2(_starting);
        //vector<double> logL3(_starting);

        for (int st1 = 0; st1 < _starting; st1++) {

            _pl[0]=0.4;
            _pl[1]=0.3;
            _pl[2]=0.2;
            _pl[3]=0.1;
            _p1(0,0,0)=0.8;
            _p1(0,1,0)=0.2;
            _p1(0,2,0)=0;
            _p1(0,3,0)=0;
            _p1(0,0,1)=0.9;
            _p1(0,1,1)=0.1;
            _p1(0,2,1)=0;
            _p1(0,3,1)=0;
            _p1(0,0,2)=0.3;
            _p1(0,1,2)=0.7;
            _p1(0,2,2)=0;
            _p1(0,3,2)=0;
            _p1(0,0,3)=0.1;
            _p1(0,1,3)=0.9;
            _p1(0,2,3)=0;
            _p1(0,3,3)=0;
            _p1(1,0,0)=0.7;
            _p1(1,1,0)=0.3;
            _p1(1,2,0)=0;
            _p1(1,3,0)=0;
            _p1(1,0,1)=0.8;
            _p1(1,1,1)=0.2;
            _p1(1,2,1)=0;
            _p1(1,3,1)=0;
            _p1(1,0,2)=0.2;
            _p1(1,1,2)=0.8;
            _p1(1,2,2)=0;
            _p1(1,3,2)=0;
            _p1(1,0,3)=0.3;
            _p1(1,1,3)=0.7;
            _p1(1,2,3)=0;
            _p1(1,3,3)=0;
            _p1(2,0,0)=0.7;
            _p1(2,1,0)=0.1;
            _p1(2,2,0)=0.2;
            _p1(2,3,0)=0;
            _p1(2,0,1)=0.8;
            _p1(2,1,1)=0.05;
            _p1(2,2,1)=0.15;
            _p1(2,3,1)=0;
            _p1(2,0,2)=0.1;
            _p1(2,1,2)=0.7;
            _p1(2,2,2)=0.2;
            _p1(2,3,2)=0;
            _p1(2,0,3)=0.15;
            _p1(2,1,3)=0.65;
            _p1(2,2,3)=0.2;
            _p1(2,3,3)=0;
            _p1(3,0,0)=0.8;
            _p1(3,1,0)=0.05;
            _p1(3,2,0)=0.15;
            _p1(3,3,0)=0;
            _p1(3,0,1)=0.1;
            _p1(3,1,1)=0.7;
            _p1(3,2,1)=0.2;
            _p1(3,3,1)=0;
            _p1(3,0,2)=0.7;
            _p1(3,1,2)=0.1;
            _p1(3,2,2)=0.2;
            _p1(3,3,2)=0;
            _p1(3,0,3)=0.2;
            _p1(3,1,3)=0.7;
            _p1(3,2,3)=0.1;
            _p1(3,3,3)=0;
            _p1(4,0,0)=0.65;
            _p1(4,1,0)=0.1;
            _p1(4,2,0)=0.1;
            _p1(4,3,0)=0.15;
            _p1(4,0,1)=0.1;
            _p1(4,1,1)=0.65;
            _p1(4,2,1)=0.1;
            _p1(4,3,1)=0.15;
            _p1(4,0,2)=0.6;
            _p1(4,1,2)=0.1;
            _p1(4,2,2)=0.2;
            _p1(4,3,2)=0.1;
            _p1(4,0,3)=0.2;
            _p1(4,1,3)=0.6;
            _p1(4,2,3)=0.1;
            _p1(4,3,3)=0.1;
            _p1(5,0,0)=0.75;
            _p1(5,1,0)=0.05;
            _p1(5,2,0)=0.05;
            _p1(5,3,0)=0.15;
            _p1(5,0,1)=0.2;
            _p1(5,1,1)=0.6;
            _p1(5,2,1)=0.05;
            _p1(5,3,1)=0.15;
            _p1(5,0,2)=0.7;
            _p1(5,1,2)=0.1;
            _p1(5,2,2)=0.1;
            _p1(5,3,2)=0.1;
            _p1(5,0,3)=0.1;
            _p1(5,1,3)=0.7;
            _p1(5,2,3)=0.15;
            _p1(5,3,3)=0.05;

            mat _f1(_n, _L);
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _f1(j, l) = 1;
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int i = 0; i < _I; i++) {
                        _f1(j, l) =(double) _f1(j, l) * _p1(i, y(j, i + 1), l);
                    }
                }
            }

            vector<double> Entropy1(_iter2);
            vector<double> completelogL1(_iter2);
            //int burn=100;
            mat plhat(_starting,_L);
            Tensor<double,4> phat(_starting, _I, _C, _L);
            for(int st=0;st<_starting;st++){
                for (int l = 0; l < _L; l++){
                    plhat(st,l)=0;
                }
            }
            for(int st=0;st<_starting;st++) {
                for (int i = 0; i < _I; i++) {
                    for (int c = 0; c < _C; c++) {
                        for (int l = 0; l < _L; l++) {
                            phat(st, i, c, l) = 0;
                        }
                    }
                }
            }

            arma::vec dir_param_pl(_L);
            for (int l = 0; l < _L; l++) {
                dir_param_pl(l) = _pl[l];
            }
            arma::cube dir_param_p(_I, _C, _L);
            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        dir_param_p(i, c, l) = _p1(i,c,l);
                    }
                }
            }

            for (int u = 0; u < _iter2; u++) {

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

                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _xj(j, l) = 0;
                    }
                }

                for (int j = 0; j < _n; j++) {
                    std::random_device rdtest1;
                    std::mt19937 generator1(rdtest1());
                    boost::random::discrete_distribution<int> distribution1(_tau.row(j).begin(), _tau.row(j).end());
                    int sample = distribution1(generator1);
                    _xj(j, sample) = 1;
                }


                vector<int> _nl(_L);
                for (int l = 0; l < _L; l++) {
                    _nl[l]=0;
                }
                for (int l = 0; l < _L; l++) {
                    for (int j = 0; j < _n; j++) {
                        _nl[l] += _xj(j,l);
                    }
                }


                for (int l = 0; l < _L; l++) {
                    dir_param_pl(l) += _nl[l];
                }
                arma::vec sample_dir_pl = rdirichlet(dir_param_pl);
                for (int l = 0; l < _L; l++) {
                    _pl[l] = sample_dir_pl(l);
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
                            _nicl(i,y(j, i + 1),l) += _xj(j, l);
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
                            _p1(i, c, l) = sample_dir_pil(c);
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
                            if (_p1(i, c, l) < 0.0001) {
                                _p1(i, c, l) = 0.001;
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
                            _f1(j, l) =(double) _f1(j, l) * _p1(i, y(j, i + 1), l);
                        }
                    }
                }

                if (u > burn) {

                    Entropy1[u] = 0;
                    for (int j = 0; j < _n; j++) {
                        for (int l = 0; l < _L; l++) {
                            Entropy1[u] += (double) -_tau(j, l) * log(_tau(j, l));
                        }
                    }

                    completelogL1[u] = 0;
                    for (int j = 0; j < _n; j++) {
                        for (int l = 0; l < _L; l++) {
                            completelogL1[u] += (double) _tau(j, l) * (log(_pl[l]) + log(_f1(j, l)));
                        }
                    }

                    for (int l = 0; l < _L; l++) {
                        plhat(st1, l) += _pl[l];
                    }

                    for (int i = 0; i < _I; i++) {
                        for (int c = 0; c < vC[i]; c++) {
                            for (int l = 0; l < _L; l++) {
                                phat(st1,i,c,l)+=_p1(i,c,l);
                            }
                        }
                    }

                }

            }


            Entropy2[st1]=0;
            completelogL2[st1]=0;
            for(int u=0;u<_iter2;u++){
                Entropy2[st1]+=Entropy1[u];
                completelogL2[st1]+=completelogL1[u];
            }
            Entropy2[st1]=Entropy2[st1]/(_iter2-burn);
            completelogL2[st1]=completelogL2[st1]/(_iter2-burn);

            logL2[st1]=Entropy2[st1]+completelogL2[st1];

            for (int l = 0; l < _L; l++) {
                plhat(st1, l) = plhat(st1, l)/(_iter2-burn);
            }

            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < vC[i]; c++) {
                    for (int l = 0; l < _L; l++) {
                        phat(st1,i,c,l)=phat(st1,i,c,l)/(_iter2-burn);
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
                        _f1(j, l) =(double) _f1(j, l) * phat(st1,i, y(j, i + 1), l);
                    }
                }
            }
            mat num0(_n, _L);
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = (double) log(_f1(j, l)) + log(plhat(st1,l));
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
                    _Tau(st1,j, l) = (double) num0(j, l) / _den0[j];
                }
            }

            double Entropy3 = 0;
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    Entropy3 += (double) -_Tau(st1,j, l) * log(_Tau(st1,j, l));
                }
            }

            double completelogL3 = 0;
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    completelogL3 += (double) _Tau(st1,j, l) * (log(plhat(st1,l)) + log(_f1(j, l)));
                }
            }

            logL3[st1]=Entropy3+completelogL3;

            cout<<"logLapprox Gibbs (a)"<<endl;
            cout<<logL3[st1]<<endl;
            cout<<"logLapprox Gibbs (b)"<<endl;
            cout<<logL2[st1]<<endl;

        }

}

}
