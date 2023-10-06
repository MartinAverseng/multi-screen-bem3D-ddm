#include <omp.h>
#include <iostream>
int main(){
double sum;
#pragma omp parallel reduction(+:sum)
for(int p = 0; p<10000; ++p){
sum += p;
}
std::cout << sum << std::endl;
std::cout << "Hello Mac world" << std::endl;
return 0;
}
