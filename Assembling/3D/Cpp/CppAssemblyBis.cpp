#include <iostream>
#include <iomanip>
#include "hypersingular.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <omp.h>

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
        Ne = (Ne + 1)/2 - 1;
        Array J = std::move(inputs[2]); 
        Array Mvtx = std::move(inputs[3]);
        Array Melt = std::move(inputs[4]);
        TypedArray<double> A = std::move(inputs[5]);
        
        
        omp_set_num_threads(8);
        
#pragma omp parallel
        {
            std::vector<int> num(6);
            std::vector<double> vtx(18); 
            std::vector<double> res(9);
            std::vector<int> T1(3),T2(3),T1bis(3),T2bis(3),I(3),sigma1(3),sigma2(3);
            std::vector<bool> b(3),b1(3),b2(3);
            std::vector<double> Aloc(Nf*Nf);
            for (int i = 0; i< Nf*Nf; ++i){
                Aloc[i] = 0;
            }
            int genFl,genFk;
            
            num[0] = 0; num[1] = 1; num[2] = 2;
#pragma omp for
            for(int el1=0; el1<Ne; ++el1){
                std::cout << el1 << std::endl;
                for(int p=0; p<3; ++p){
                    T1[p] = (int)Melt[el1][p];
                    T1bis[p] = (int)Melt[el1+Ne][p];
                    for(int q=0; q<3; ++q){
                        vtx[3*p+q] = Mvtx[T1[p]][q];
                    }
                }
                
                for(int el2=0; el2<Ne; ++el2){
                    
                    for(int p=0; p<3; ++p){
                        T2[p] = (int)Melt[el2][p];
                        T2bis[p] = (int)Melt[el2+Ne][p];
                        for(int q=0; q<3; ++q){
                            vtx[3*(3+p)+q] = Mvtx[T2[p]][q];
                        }
                    }
                    
                    ismember(T2,T1,b,I);
                    ismember(T1bis,T1,b1,sigma1);
                    std::cout << "sigma1 1 = " << sigma1[0];
                    std::cout << "sigma1 2 = " << sigma1[1];
                    std::cout << "sigma1 3 = " << sigma1[2] << std::endl;
                    ismember(T2bis,T2,b2,sigma2);
                    std::cout << "sigma2 1 = " << sigma2[0];
                    std::cout << "sigma2 2 = " << sigma2[1];
                    std::cout << "sigma2 3 = " << sigma2[2] << std::endl;
                    
                    num[3]=3;num[4]=4;num[5]=5;
                    for(int p=0; p<3; ++p){
                        if(b[p]){num[3+p]=I[p];}
                    }
                    HsOp(vtx.data(),num.data(),res.data());
                    
                    
                    for (int k = 0; k < 3; ++k){
                        for(int l= 0; l < 3; ++l){
                            genFk = (int)J[el1][k];
                            genFl = (int)J[el2][l];
                            Aloc[genFk+Nf*genFl] = Aloc[genFk + Nf*genFl] + res[3*l + k];
                            genFk = (int)J[el1+Ne][k];
                            genFl = (int)J[el2+Ne][l];
                            Aloc[genFk+Nf*genFl] = Aloc[genFk + Nf*genFl] + res[3*sigma2[l] + sigma1[k]];
                        }
                    }
                    
                    
                    
                    
                }
            }
#pragma omp critical
            {
                for (int i = 0; i<Nf; ++i){
                    for (int j=0; j<Nf; ++j){
                        A[i][j] = A[i][j] + Aloc[i*Nf + j];
                    }
                    
                }
            }
            
            
            
            
            outputs[0] = A;
        }
    };
};
