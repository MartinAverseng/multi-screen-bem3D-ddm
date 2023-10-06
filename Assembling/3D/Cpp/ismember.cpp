
#include <iostream>
#include <iomanip>
#include <vector>

int main(){
    std::vector<int> A(3),B(3),I(3);
    A[0] = 10; A[1] = 24; A[2] = 76;
    B[0] = 56; B[1] = 76; B[2] = 7;
    std::vector<bool> b(3);
    for (int p = 0; p<3; ++p){
        b[p] = false;
        I[p] = 0;
        for (int q = 0; q < 3; ++q){
            if (B[p] == A[q]){
                b[p] = true;
                I[p] = q;
            }
        }
    }
    std::cout << "I : " << I[0] << " " << I[1] << " " << I[2] << std::endl;
};
