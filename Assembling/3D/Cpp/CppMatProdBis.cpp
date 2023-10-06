#include <iostream>
#include <iomanip>
#include "quadrule.hpp"
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
        int Ne = (int) inputs[1][0];
        int Nvtx = (int) inputs[2][0];
        const Array J = inputs[3];
        const Array MMvtx = inputs[4];
        const Array Melt = inputs[5];
        const TypedArray<double> vv = inputs[6];
        std::vector<double> Mvtx(3*Nvtx);
        std::vector<int> Jint(3*Ne);
        std::vector<int> MeltInt(3*Ne);
        std::vector<double> v(Nf);
        std::vector<double> Wvloc(Nf);
        for (int el = 0; el < Ne; ++el){
            for (int p =0; p < 3; ++p){
                Jint[el*3+p]=J[el][p];
                MeltInt[el*3+p]=Melt[el][p];
            }
        }
        for (int vt = 0; vt < Nvtx; ++vt){
            for (int p =0; p < 3; ++p){
                Mvtx[vt*3+p]=MMvtx[vt][p];
            }
        }
        //std::vector<double> vals(Nf*Nf,0.);
//         TypedArray<double> A = std::move(inputs[5]);
        
        for (int i = 0; i< Nf; ++i){
            Wvloc[i] = 0;
            v[i] = vv[i];
        }
        
        omp_set_num_threads(8);
#pragma omp parallel
        {
            std::vector<int> num(6);
            std::vector<double> vtx(18);
            std::vector<double> res(9);
            std::vector<int> T1(3),T2(3),I(3);
            std::vector<bool> b(3);
            int genFl,genFk;
            int index1,index2,index3,index4,index5,index6,index7;
            
            num[0] = 0; num[1] = 1; num[2] = 2;
//             std::cout << "MV prod" << std::endl;
#pragma omp for
            for(int el1=0; el1<Ne; ++el1){
                std::cout << el1 << std::endl;
                index1 = 3*el1; // Will be reused below
                for(int p=0; p<3; ++p){
                    T1[p] = MeltInt[index1+p]; // Don't change index1 here
                    for(int q=0; q<3; ++q){
                        vtx[3*p+q] = Mvtx[3*T1[p]+q];
                    }
                }
                
                for(int el2=0; el2<Ne; ++el2){
                    index2 = 3*el2;
                    for(int p=0; p<3; ++p){
                        T2[p] = MeltInt[index2];
                        index2++;
                        index3 = 3*(3+p);
                        for(int q=0; q<3; ++q){
                            vtx[index3] = Mvtx[3*T2[p]+q];
                            index3++;
                        }
                    }
                    
                    ismember(T2,T1,b,I);
                    index4=0;
                    for(int p=3; p<6; ++p){
                        if(b[index4]){
                            num[p]=I[index4];
                        }
                        else{
                            num[p]=p;
                        }
                        index4++;
                    }
                    HsOp(vtx.data(),num.data(),res.data());
                    index5=3*el2;
                    index6 = index1;
                    for (int k = 0; k < 3; ++k){
                        genFk = Jint[index6]; //
                        index6++;
                        for(int l= 0; l < 3; ++l){
                            genFl = Jint[index5];
                            index5++;
#pragma omp atomic
                            Wvloc[genFk] += res[3*l + k]*v[genFl];
                        }
                        index5-=3;//index5 = 3*el2+l
                    }
                }
            }
            
        }
        outputs[0] = factory.createArray({Nf,1},Wvloc.begin(),Wvloc.end());
    }
};
