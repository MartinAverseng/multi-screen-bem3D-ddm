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
        long unsigned int Nsub = (int) inputs[0][0];
        int NsubElt = (int) inputs[1][0];
        int Nf = (int) inputs[2][0];
        int Ne = (int) inputs[3][0];
        int Nvtx = (int) inputs[4][0];
        const Array J = inputs[5];
        const Array MMvtx = inputs[6];
        const Array Melt = inputs[7];
        const Array KKlist = inputs[8];
        const Array CCppDic = inputs[9];
        
        std::vector<double> Mvtx(3*Nvtx);
        std::vector<int> Jint(3*Ne);
        std::vector<int> MeltInt(3*Ne);
        std::vector<int> Klist(NsubElt);
        std::vector<int> CppDic(Nf);
        
        
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
        for (int elsub=0; elsub<NsubElt; elsub++){
            Klist[elsub] = KKlist[elsub];
        }
        for (int i=0;i<Nf;++i){
            CppDic[i] = CCppDic[i];
        }
        std::vector<double>Aloc(Nsub*Nsub);
        for (int i = 0; i < Nsub*Nsub; ++i){
            Aloc[i] = 0;
        }
        omp_set_num_threads(8);
        
#pragma omp parallel
        {
            
            std::vector<int> num(6);
            std::vector<double> vtx(18);
            std::vector<double> res(9);
            std::vector<int> T1(3),T2(3),I(3);
            std::vector<bool> b(3);
            int genFl,genFk,el1,el2,genFkSub,genFlSub;
            int index1,index2,index3,index4,index5,index6,index7;
            
            num[0] = 0; num[1] = 1; num[2] = 2;
#pragma omp for
            for(int el1sub=0; el1sub<NsubElt; ++el1sub){
                el1=Klist[el1sub];
                index1 = 3*el1; // Will be reused below
                for(int p=0; p<3; ++p){
                    T1[p] = MeltInt[index1+p]; // Don't change index1 here
                    for(int q=0; q<3; ++q){
                        vtx[3*p+q] = Mvtx[3*T1[p]+q];
                    }
                }
                
                for(int el2sub=0; el2sub<NsubElt; ++el2sub){
                    el2=Klist[el2sub];
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
                        genFk = Jint[index6];
                        genFkSub = CppDic[genFk];
                        index6++;
                        for(int l= 0; l < 3; ++l){
                            genFl = Jint[index5];
                            genFlSub = CppDic[genFl];
                            index5++;
                            if (genFkSub>=0 && genFlSub >=0){
                                index7=genFkSub+Nsub*genFlSub;
# pragma omp atomic
                                Aloc[index7] += res[3*l + k];
                            }
                        }
                        index5-=3;//index5 = 3*el2+l
                    }
                }
            }
            
        }
        outputs[0] = factory.createArray({Nsub,Nsub},Aloc.begin(),Aloc.end());
    }
};

