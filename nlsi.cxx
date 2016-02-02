#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void diff(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N)
{

  cmplx* d=new cmplx[N];

  cmplx u = cmplx(0, dt/(dx*dx)); // u=-alpha
  for(int i=0;i<N;i++) d[i] = 1.0 - 2.0*u;
  
  for(int i=1;i<N;i++) {
    d[i] -= u / d[i-1] * u;
    f0[i] -= u / d[i-1] * f0[i-1];
  }
  
  f1[N-1] = f0[N-1] / d[N-1];
  for(int i=N-2;i>=0;i--) f1[i] = ( f0[i] - u * f1[i+1] ) / d[i];

  delete[] d;
}
void abs2(cmplx* const f1, cmplx* const f0,
          const double dt, const int N)
{
  for (int i=0; i<N; i++) {
    cmplx arg = cmplx(0, -norm(f0[i])*dt);
    f1[i] = f0[i] * exp(arg);
  }
}
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) {
		  
		  diff(psi1, psi0, dt, dx, Nx);
		  abs2(psi1, psi1, dt, Nx);
		  h = psi0;
		  psi0 = psi1;
		  psi1 = h;
		  
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}

	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
