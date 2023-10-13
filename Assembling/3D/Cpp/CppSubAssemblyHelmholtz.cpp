#include <iostream>
#include <iomanip>
#include "quadrule.hpp"
#include "hypersingular_helmholtz.hpp"
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
        double wavenum = inputs[10][0];
        
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
        std::vector<double>Aloc_real(Nsub*Nsub);
        std::vector<double>Aloc_imag(Nsub*Nsub);
        for (int i = 0; i < Nsub*Nsub; ++i){
            Aloc_real[i] = 0;
            Aloc_imag[i] = 0;
        }
        omp_set_num_threads(8);
        
#pragma omp parallel
        {
            
            std::vector<int> num(6);
            std::vector<double> vtx(18);
            std::vector<double> res_real(9);
            std::vector<double> res_imag(9);
            std::vector<int> T1(3),T2(3),I(3);
            std::vector<bool> b(3);
            int genFl,genFk,el1,el2,genFkSub,genFlSub;
            int index1,index2,index3,index4,index5,index6,index7;
            
            num[0] = 0; num[1] = 1; num[2] = 2;
#pragma omp for
            for(int el1sub=0; el1sub<NsubElt; ++el1sub){
                el1=Klist[el1sub];
                //std::cout << el1 << std::endl;
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
                    HelmholtzHsOp(vtx.data(),num.data(),res_real.data(),res_imag.data(),wavenum);
                    index5=3*el2;
                    index6 = index1;
                    for (int k = 0; k < 3; ++k){
                        genFk = Jint[index6];
                        genFkSub = CppDic[genFk];
//
                        index6++;
                        for(int l= 0; l < 3; ++l){
                            genFl = Jint[index5];
                            genFlSub = CppDic[genFl];
                            index5++;
                            if (genFkSub>=0 && genFlSub >=0){
                                index7=genFkSub+Nsub*genFlSub;
# pragma omp atomic
                                Aloc_real[index7] += res_real[3*l + k];
                                Aloc_imag[index7] += res_imag[3*l + k];
                            }
                        }
                        index5-=3;//index5 = 3*el2+l
                    }
                }
            }
            
        }
        outputs[0] = factory.createArray({Nsub,Nsub},Aloc_real.begin(),Aloc_real.end());
        outputs[1] = factory.createArray({Nsub,Nsub},Aloc_imag.begin(),Aloc_imag.end());
    }
};


// #pragma omp parallel
//         {
//             std::vector<int> num(6);
//             std::vector<double> vtx(18);
//             std::vector<double> res(9);
//             std::vector<int> T1(3),T2(3),I(3);
//             std::vector<bool> b(3);
//             std::vector<double> Asubloc(Nsub*Nsub);
//             for (int i = 0; i< Nsub*Nsub; ++i){
//                 Asubloc[i] = 0;
//             }
//             int place;
//             double Rik, Rjl;
//             int genFl,genFk,genFlSub,genFkSub;
//             int el1, el2;
//
//             num[0] = 0; num[1] = 1; num[2] = 2;
// #pragma omp for
//             for(int el1ind=0; el1ind<Ne; ++el1ind){
//                 el1 = Klist[el1ind];
//                 std::cout << el1 << std::endl;
//                 for(int p=0; p<3; ++p){
//                     T1[p] = (int)Melt[el1][p];
//                     for(int q=0; q<3; ++q){
//                         vtx[3*p+q] = Mvtx[T1[p]][q];
//                     }
//                 }
//
//                 for(int el2ind=0; el2ind<Ne; ++el2ind){
//                     el2 = Klist[el2ind];
//                     for(int p=0; p<3; ++p){
//                         T2[p] = (int)Melt[el2][p];
//                         for(int q=0; q<3; ++q){
//                             vtx[3*(3+p)+q] = Mvtx[T2[p]][q];
//                         }
//                     }
//                     ismember(T2,T1,b,I);
//                     num[3]=3;num[4]=4;num[5]=5;
//                     for(int p=0; p<3; ++p){
//                         if(b[p]){num[3+p]=I[p];}
//                     }
//                     HsOp(vtx.data(),num.data(),res.data());
//
//
//                     for (int k = 0; k < 3; ++k){
//                         for(int l= 0; l < 3; ++l){
//                             genFk = (int)J[el1][k];
//                             genFl = (int)J[el2][l];
//                             genFkSub = CppDic[genFk];
//                             genFlSub = CppDic[genFl];
//                             if (genFkSub>=0 && genFlSub >=0){
//                                 place = genFkSub + Nsub*genFlSub;
//                                 Asubloc[place] = Asubloc[place] + res[3*l + k];
//                             }
//                         }
//                     }
//                 }
//             }
// #pragma omp critical
//             {
//                 for (int i = 0; i<Nsub; ++i){
//                     for (int j=0; j<Nsub; ++j){
//                         Asub[i][j] = Asub[i][j] + Asubloc[i*Nsub + j];
//                     }
//
//                 }
//             }
//
//
//
//
//             outputs[0] = Asub;
//         }
//     };
// };
