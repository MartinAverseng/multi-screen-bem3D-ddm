#include <iostream>
#include <iomanip>
#include "quadrule.hpp"
#include "hypersingular.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"
# include <omp.h>
using namespace std;

using namespace matlab::data;
using namespace matlab::mex;



void ismember(const std::vector<int>& A,
        const std::vector<int>& B,
        std::vector<bool>& b,
        std::vector<int>& I){
    for (int p = 0; p<3; ++p){
        b[p] = false;
        I[p] = -1;
        for (int q = 0; q < 3; ++q){
            if (A[p] == B[q]){
                b[p] = true;
                I[p] = q;
            }
        }
    }
}

class MexFunction: public matlab::mex::Function {
    matlab::data::ArrayFactory factory;
public:
    void operator()(ArgumentList outputs,ArgumentList inputs){
        long unsigned int Nf = (int) inputs[0][0];
        long unsigned int Ne = (int) inputs[1][0];
        Array J = std::move(inputs[2]); //
        Array Mvtx = std::move(inputs[3]);
        Array Melt = std::move(inputs[4]);
        //std::vector<double> vals(Nf*Nf,0.);
        TypedArray<double> v = std::move(inputs[5]);//factory.createArray({Nf,Nf},vals.begin(),vals.end());
        TypedArray<double> Wv = std::move(inputs[6]);//factory.createArray({Nf,Nf},vals.begin(),vals.end());
        // Algo remplissage de conteneur. Marche sur vector, list, map...
        // Et aussi n'importe quel type qui ait des itérateurs.
        // Voir normes C++ et version gcc. Flag -std=c++17
        
        /*for (int gf1 = 0; gf1 < Nf; ++gf1){
         * for (int gf2 = 0; gf2 < Nf; ++gf2){
         * A[gf1][gf2] = 0.;
         * }
         * }*/
        
        omp_set_num_threads(8);
#pragma omp parallel
        {
            std::vector<int> num(6);
            std::vector<double> vtx(18); // = double* qui sait comment gérer sa mém
            std::vector<double> res(9);
            std::vector<int> T1(3),T2(3),I(3);
            std::vector<bool> b(3);
            num[0] = 0; num[1] = 1; num[2] = 2;
            std::vector<double> Wv_local(Nf);
            int genFl,genFk;
            for (int i = 0; i<Nf; i++)
            {
                Wv_local[i] = 0;
            }
#pragma omp for
            for(int el1=0; el1<Ne; ++el1){
                std::cout << el1 << std::endl;
                for(int p=0; p<3; ++p){
                    T1[p] = (int)Melt[el1][p];
                    for(int q=0; q<3; ++q){
                        vtx[3*p+q] = Mvtx[T1[p]][q];
                    }
                }
                
                for(int el2=0; el2<Ne; ++el2){
                    for(int p=0; p<3; ++p){
                        T2[p] = (int)Melt[el2][p];
                        for(int q=0; q<3; ++q){
                            vtx[3*(3+p)+q] = Mvtx[T2[p]][q];
                        }
                    }
                    
                    ismember(T2,T1,b,I);
                    num[3]=3;num[4]=4;num[5]=5;
                    for(int p=0; p<3; ++p){
                        if(b[p]){num[3+p]=I[p];}
                    }
                    HsOp(vtx.data(),num.data(),res.data());
                    
                    for (int k = 0; k < 3; ++k){
                        for(int l= 0; l < 3; ++l){
                            genFk = (int)J[el1][k];
                            genFl = (int)J[el2][l];
                            
                            Wv_local[genFl]=Wv_local[genFl]+res[3*l + k]*v[genFk];
                        }
                    }
                    
                    
                }
                
                
            }
#pragma omp critical
            {
                for(int i=0; i<Nf; ++i){
                    Wv[i]=Wv[i] + Wv_local[i];
                }
            }
        }
        
        
        outputs[0] = Wv;
        
    }
};
