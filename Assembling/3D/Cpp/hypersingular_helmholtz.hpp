#ifndef HYPERSINGULAR_HPP
#define HYPERSINGULAR_HPP

#include <cmath>

double  PI      = 3.141592653589793;
double* rule[4] = {rule0+1,rule1+1,rule2+1,rule3+1};
int     nq[4]   = {(int)*rule0, (int)*rule1, (int)*rule2, (int)*rule3};

class R3{

private:

  double x[3];

public:

  R3(const double& x0=0.,const double& x1=0.,const double& x2=0.){
    x[0]=x0;x[1]=x1;x[2]=x2;
  }

  R3(const R3& p){x[0]=p.x[0];x[1]=p.x[1];x[2]=p.x[2];
  }

  R3& operator=(const R3& p){if(this!=&p){
    x[0]=p.x[0];x[1]=p.x[1];x[2]=p.x[2];}
    return *this;
  }

  double& operator[](const int& j){
    return x[j];
  }
  const double& operator[](const int& j) const {
    return x[j];
  }

  R3& operator+=(const R3& u){
    x[0]+=u[0];x[1]+=u[1];x[2]+=u[2]; return *this;
  }
  friend R3 operator+(R3 u, const R3& v){
    return u+=v;
  }

  R3& operator-=(const R3& u){
    x[0]-=u[0];x[1]-=u[1];x[2]-=u[2]; return *this;
  }
  friend R3 operator-(R3 u, const R3& v){
    return u-=v;
  }

  R3& operator*=(const double& a){
    x[0]*=a; x[1]*=a; x[2]*=a; return *this;
  }
  friend R3 operator*(const double& a, R3 u){
    return u*=a;
  }

  friend double operator,(const R3& u, const R3& v){
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
  }

  friend double norm(const R3& u){
    return sqrt( (u,u) );
  }

  friend void normalize(R3& u){
    double a=norm(u); u*=1./a;
  }

};

R3 vprod(const R3& u, const R3& v)
{
  return R3(u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]);
};


void HelmholtzHsOp(double* vtx, int* tri, double* res_real,double* res_imag,double wavenum){
  // We compute curls and (Nx,Ny) before permutation,
  // to avoid dealing with signs.

  // Then we permute and compute the -k^2 int_{T1xT2} phi_j(x) psi_k(y)Gk(x,y) dx dy
  double k2 = wavenum*wavenum;

  int* Jx = tri+3;
  int* Jy = tri;
  R3 x[3], y[3];
  for(int j=0; j<3; ++j){
    for(int k=0; k<3; ++k){
      x[j][k]=vtx[3*Jx[j]+k];
      y[j][k]=vtx[3*Jy[j]+k];
    }
  }
  //

  R3 Nx,Ny;
  Nx = vprod(x[1]-x[0],x[2]-x[1]);
  Ny = vprod(y[1]-y[0],y[2]-y[1]);
  double dx = norm(Nx); normalize(Nx);
  double dy = norm(Ny); normalize(Ny);
  double NxNy = (Nx,Ny);

  R3 nx[3], ny[3];
  for(int j=0; j<3; ++j){
    nx[(j+2)%3]=vprod(x[(j+1)%3]-x[j],Nx);
    ny[(j+2)%3]=vprod(y[(j+1)%3]-y[j],Ny);
  }

  R3 curlx[3], curly[3];
  for(int j=0; j<3; ++j){
    curlx[j] = (1./(x[(j+1)%3]-x[j],nx[j]))*vprod(nx[j],Nx);
    curly[j] = (1./(y[(j+1)%3]-y[j],ny[j]))*vprod(ny[j],Ny);
  }

  // Now compute permutation

  int Iy[] = {tri[0],tri[1],tri[2]};
  int Ix[] = {tri[3],tri[4],tri[5]};
  int sigma1[] = {0,1,2};
  int sigma2[] = {0,1,2};
  int sigma1inv[] = {0,1,2};
  int sigma2inv[] = {0,1,2};
  int r = 0;
  for (int j = 0; j< 3; ++j){
    for (int k = 0; k < 3; ++k){
      if (Ix[j]==Iy[k]){
        sigma1[sigma1inv[j]] = sigma1[r];
        sigma1inv[sigma1[r]] = sigma1inv[j];
        sigma1[r] = j;
        sigma1inv[j] = r;

        sigma2[sigma2inv[k]] = sigma2[r];
        sigma2inv[sigma2[r]] = sigma2inv[k];
        sigma2[r] = k;
        sigma2inv[k] = r;
        r++;
      }
    }
  }
  // Now apply the permutation
  for (int j = 0; j < 3; ++j){
    Ix[j] = tri[3+sigma1[j]];
    Iy[j] = tri[sigma2[j]];
  }


  R3 xsigma[3], ysigma[3];
  R3 nxsigma[3], nysigma[3];
  for(int j=0; j<3; ++j){
    for(int k=0; k<3; ++k){
      xsigma[j][k]=vtx[3*Ix[j]+k];
      ysigma[j][k]=vtx[3*Iy[j]+k];
    }
  }
  R3 Nxsigma,Nysigma;
  Nxsigma = vprod(xsigma[1]-xsigma[0],xsigma[2]-xsigma[1]);
  Nysigma = vprod(ysigma[1]-ysigma[0],ysigma[2]-ysigma[1]);
  for(int j=0; j<3; ++j){
    nxsigma[(j+2)%3]=vprod(xsigma[(j+1)%3]-xsigma[j],Nxsigma);
    nysigma[(j+2)%3]=vprod(ysigma[(j+1)%3]-ysigma[j],Nysigma);
  }


  double* qr = rule[r];
  double inter_real=0; double inter_imag=0; R3 Xq,Yq;
  double phijq, psikq, dj, dk, rXqYq;
  double inc_real = 0;
  double inc_imag = 0;
  double res1_real[] = {0,0,0,0,0,0,0,0,0};
  double res1_imag[] = {0,0,0,0,0,0,0,0,0};
  for(int q=0; q<nq[r]; q+=5){

    const double& u0 = qr[q  ];
    const double& u1 = qr[q+1];
    const double& v0 = qr[q+2];
    const double& v1 = qr[q+3];
    const double& wq = qr[q+4];

    Xq = xsigma[0]+u0*(xsigma[1]-xsigma[0])+u1*(xsigma[2]-xsigma[0]);
    Yq = ysigma[0]+v0*(ysigma[1]-ysigma[0])+v1*(ysigma[2]-ysigma[0]);

    rXqYq = norm(Xq-Yq);
    inc_real = wq/rXqYq*cos(wavenum*rXqYq);
    inc_imag = wq/rXqYq*sin(wavenum*rXqYq);
    inter_real += inc_real;
    inter_imag += inc_imag;


    for (int j = 0; j<3; ++j){
      dj = (xsigma[(j+1)%3] - xsigma[j],nxsigma[j]);
      phijq = (xsigma[(j+1)%3] - Xq,nxsigma[j])/dj;
      for (int k = 0; k<3; ++k){
        dk = (ysigma[(k+1)%3] - ysigma[k],nysigma[k]);
        psikq = (ysigma[(k+1)%3] - Yq,nysigma[k])/dk;
        res1_real[3*j + k] += inc_real*phijq*psikq;
        res1_imag[3*j + k] += inc_imag*phijq*psikq;
      }
    }
  }
  double prefactor1 = dx*dy/(4.*PI);
  double prefactor2 = NxNy*k2*prefactor1;
  inter_real = inter_real*prefactor1;
  inter_imag = inter_imag*prefactor1;


  for(int j=0; j<3; ++j){
    for(int k=0; k<3; ++k){
      res_real[3*j+k] = inter_real*(curlx[j],curly[k]) - res1_real[3*sigma1inv[j]+sigma2inv[k]]*prefactor2;
      res_imag[3*j+k] = inter_imag*(curlx[j],curly[k]) - res1_imag[3*sigma1inv[j]+sigma2inv[k]]*prefactor2;
    }
  }

};


#endif
