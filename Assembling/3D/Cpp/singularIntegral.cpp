// Calls method HsOp in hypersinglar.hpp to compute the local hypersingular 
// matrix over two triangles.

#include <iostream>
#include <iomanip>
#include "hypersingular.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace std;
using namespace matlab::data;
using namespace matlab::mex;
class MexFunction : public matlab::mex::Function {
    matlab::data::ArrayFactory factory;
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        TypedArray<double> vtx = inputs[0];
        double vtxCpp[18];
        for (int i = 0; i < 6; i++){
            for (int j = 0; j<3;j++){
                vtxCpp[3*i + j] = vtx[i][j];
            }
        }

        TypedArray<double> tri = inputs[1];
        int triCpp[6];
        for (int i = 0; i< 2; i++){
            for (int j = 0; j < 3; j++){
                triCpp[3*i+j] = (int)tri[i][j];
            }
        }

        double res[9];

        HsOp(vtxCpp,triCpp,res);

        //transpose(res,3);
        std::vector<double> vec(res,res+9);
        ArrayFactory factory;
        TypedArray<double>  A = factory.createArray( {3,3},vec.begin(),vec.end());
        outputs[0] = A;
    }


};
