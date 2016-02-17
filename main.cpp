#include <cmath>
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <iostream>
#include <fstream>
#include <omp.h>

using namespace std;
typedef complex <double> cdouble;

/* Define physical parameters of simulation (point (a, b) in phasespace) */
double a = 1.2, b = 0.1;
double g = 1.;
double beta = 0.00005;
double R = 1.;
double gam_r = 1.;
double g_r = 2. * g;
double gam_c = b * gam_r;
double P_th = gam_c * gam_r / R;
double P = a * P_th;
double W = 1.0;

/* CGLE-specific values */
double phi_0 = sqrt((P-P_th)/gam_c);
double n0 = phi_0 * phi_0;
double ns = n0 / (P/gam_c - 1);
double D = gam_c;
double eta = (1 + ns / n0) / gam_c;
double cgle_var = 4 * D * (eta * eta) / n0;

/* Define numerical parameters of simulation */
double L = 600;
double T = 20000;
//double dx = 0.266667;
//double dt = 0.005;
 double dx = 0.4;
 double dt = 0.01;
int Nx = (int) (L / dx);
int Nt = (int) (T / dt);
double x0 = -L;
int freq = 0.02 * T;
int Ns = 2;

/* Some usefull stuff */
const complex<double> I(0,1);
const double TWO_PI = M_PI * 2.;
const double dxx = dx * dx;
const double thld = 0.8;

/* Noise statistics */
static cdouble n_mean = 0;
static cdouble n_var = 0;

/* Wigner noise variance */
double wigner_variance(){
  return 2 * beta * dt / dx * sqrt(D);
}

/* Gaussian noise generation - Natalia's work */
double gaussian_noise(const double &var)
{
	static bool haveSpare = false;
	static double rand1, rand2;

	if(haveSpare)
	{
		haveSpare = false;
		return sqrt(var * rand1) * sin(rand2);
	}

	haveSpare = true;

	rand1 = rand() / ((double) RAND_MAX);
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;
	return sqrt(var * rand1) * cos(rand2);
}

void f(vector<cdouble> &y, const vector<cdouble> &y0) {
  for (int ix = 1; ix < Nx-1; ix++){
    cdouble psi_sq = abs(y0[ix]) * abs(y0[ix]);
    cdouble derX = -0.5 * (y0[ix-1] - 2.0 * y0[ix] + y0[ix+1]) / dxx;
    cdouble cgle = cdouble(0,1) * 0.5 * (P / (1. + psi_sq / ns) - gam_c) * y0[ix];
    cdouble nlt = g * psi_sq * y0[ix];

    y[ix] = -cdouble(0, 1) * (derX + nlt + cgle);
  }

  // periodic boundary conditions

  y[0] = y[Nx-2];
  y[Nx-1] = y[1];

}

/* Compute norm */
double norm(const vector<cdouble> &y) {
    int N = y.size();

    double sum = 0;
    for (int i = 0; i < N; i++){
        sum += abs(y[i]) * abs(y[i]);
        }

    return sum/N;
}

/* Calculate correlation function */
vector<cdouble> corr(vector<cdouble>  &y) {
    int N = y.size();
    vector<cdouble> cfn(N, 0);

    // calculate norm squared
    double sum = 0;

    for (int i = 0; i < N; i++){
        sum += abs(y[i]) * abs(y[i]);
        }

    complex<double> mean = 0;

    // calculate correlation
    double sh = 0;

    for(int d = 0; d < N; d++)
    {

        for(int ix = 0; ix < N; ix++)
        {
            sh = ix + d;
            if (sh >= N)
                sh = sh - (N - 1);

            mean = mean + conj(y[ix]) * y[int(sh)];
        }

        double i = d - N/2;

        //if (i >= 0) cfn[i] = abs(mean/sum);
        //else cfn[d+N/2] = abs(mean/sum);

        if (i >= 0) cfn[i] = mean/sum;
        else cfn[d+N/2] = mean/sum;

	mean = 0;
    }

    return cfn;
}

void n_stats(vector <double> &y_re, vector <double> &y_im){

    //vector <double> stats(2, 0);

    cdouble abs = 0.;
    cdouble acc = 0.;
    cdouble acc2 = 0.;

    int N = y_re.size();


    for (int i = 0; i < N; i++) {
        abs = y_re[i] +  I * y_im[i];
        acc += abs;
        acc2 += acc * conj(acc);
    }

    n_mean = acc / double(N);
    n_var =  acc2 / double(N) - n_mean * n_mean;
}

/**************************************************************/
/************************   MAIN   ****************************/
/**************************************************************/

int main(){

  system("clear");
  srand(time(NULL));
  double start_time = omp_get_wtime();		// total simulation time (start)
  double time_s = 0;				// simulation time (iter_start)
  double time_e = 0;				// simulation time (iter_end)

/* Working variables */

  vector<cdouble> psi0(Nx, 0);
  vector<double> phase(Nx, 0);
  vector<cdouble> wig(Nx, 0);

  vector<double> fluct_total(Nx, 0);
  vector<double> fluct_noise(Nx, 0);
  vector<double> fluct_kinetic(Nx, 0);

/* Measurements */

  vector<double> no(Nt, 0);
  vector<cdouble> g1d(Nx, 0);

/*********************** Simulation loop ***********************/
for (int is = 1; is < Ns; is++)
{
    cout << "****************************************************" << endl;
    cout << "Simulation " << is << "/" << Ns - 1 << " - START" << endl;
    cout << "****************************************************" << endl;

    time_s = omp_get_wtime();

/* Define file handles for current simulation */
    ofstream ofile;
    ofstream oinit;
    ofstream ostats;
    ofstream ophase;
    ofstream ocorr;
    ofstream omodule;
    ofstream oinfo;
    ofstream ofluct;

    string fname = "S" + to_string(is);
    string fcommand = "mkdir " + fname;
    system(fcommand.c_str());

    ofile.open(fname + "/output.txt");
    ostats.open(fname + "/stats.txt");
    ophase.open(fname + "/phase.txt");
    omodule.open(fname + "/module.txt");
    oinit.open(fname + "/initial.txt");
    ocorr.open(fname + "/corr.txt");
    oinfo.open(fname + "/info.txt");
    ofluct.open(fname + "/fluct.txt");

    // copy gnuplot scripts

    string name = "cp scripts/corr.sh " + fname + "/corr.sh";
    system(name.c_str());

    name = "cp scripts/phase.sh " + fname + "/phase.sh";
    system(name.c_str());

    name = "cp scripts/module.sh " + fname + "/module.sh";
    system(name.c_str());

    name = "cp scripts/plot.sh " + fname + "/plot.sh";
    system(name.c_str());

    name = "cp scripts/stats.sh " + fname + "/stats.sh";
    system(name.c_str());

    name = "cp scripts/fluct.sh " + fname + "/fluct.sh";
    system(name.c_str());

    name = "cp scripts/animate.dem " + fname + "/animate.dem";
    system(name.c_str());

    name = "cp scripts/gnuplot.rot " + fname + "/gnuplot.rot";
    system(name.c_str());

/****************************************************************/
/* Set initial wave function and impose boundary conditions */
/****************************************************************/
double var = 0;

for (int ix = 0; ix < Nx; ix++){

  var = wigner_variance();

  psi0[ix] = phi_0 + 0.5 * gaussian_noise(var) + 0.5 * I * gaussian_noise(var);

}

  // imposing periodic boundary conditions

  psi0[0] = psi0[Nx-2];
  psi0[Nx-1] = psi0[1];

/****************************************************************/
/* Runge-Kutta 4th order algorithm */
/****************************************************************/

  // Allocate memory for coefficients
  vector<cdouble> psi(Nx, 0);
  vector<double> w_x(Nx, 0);
  vector<double> w_y(Nx, 0);
  vector<double> w_x_past(Nx, 0);
  vector<double> w_y_past(Nx, 0);

  vector<cdouble> k1(Nx, 0);
  vector<cdouble> k2(Nx, 0);
  vector<cdouble> k3(Nx, 0);
  vector<cdouble> k4(Nx, 0);

  vector<cdouble> psi_k1(Nx, 0);
  vector<cdouble> psi_k2(Nx, 0);
  vector<cdouble> psi_k3(Nx, 0);
  vector<cdouble> psi_k4(Nx, 0);

/****************************************************************/
/* Start dynamics */
/****************************************************************/

  for (int it = 0; it < Nt; it++) {
    double tt = it * dt;

    /*** RK4 step 1 ***/
    f(psi, psi0);

    for (int ix = 0; ix < Nx; ix++){
      k1[ix] = psi[ix];
      psi_k1[ix] = psi0[ix] + 0.5 * k1[ix] * dt;
    }

    // imposing boundary conditions

    k1[0] = k1[Nx-2];
    k1[Nx-1] = k1[1];

    /*** RK4 step 2 ***/
    f(psi, psi_k1);

    for (int ix = 0; ix < Nx; ix++) {
      k2[ix] = psi[ix];
      psi_k2[ix] = psi0[ix] + 0.5 * k2[ix] * dt;
    }

    // imposing boundary conditions

    k2[0] = k2[Nx-2];
    k2[Nx-1] = k2[1];

    /*** RK4 step 3 ***/
    f(psi, psi_k2);

    for (int ix = 0; ix < Nx; ix++) {
      k3[ix] = psi[ix];
      psi_k3[ix] = psi0[ix] + k3[ix] * dt;
    }

    // imposing boundary conditions

    k3[0] = k3[Nx-2];
    k3[Nx-1] = k3[1];

    /*** RK4 step 3 ***/
    f(psi, psi_k3);

    for (int ix = 0; ix < Nx; ix++) {
      k4[ix] = psi[ix];
    }

    // imposing boundary conditions

    k4[0] = k4[Nx-2];
    k4[Nx-1] = k4[1];

    // compute wigner noise
    for (int ix = 0; ix < Nx; ix++) {
      var = wigner_variance();

      w_x_past[ix] = w_x[ix];
      w_y_past[ix] = w_y[ix];

      w_x[ix] = gaussian_noise(var);
      w_y[ix] = gaussian_noise(var);

    }

    // imposing boundary conditions
    w_x[0] = w_x[Nx-2];
    w_x[Nx-1] = w_x[1];

    w_y[0] = w_y[Nx-2];
    w_y[Nx-1] = w_y[1];

    // update condensate state
    for (int ix = 0; ix < Nx; ix++) {
      psi[ix] = psi0[ix] + (k1[ix] + 2.0 * k2[ix] + 2.0 * k3[ix] + k4[ix]) * dt / 6.0
	 	+ (w_x[ix] + I * w_y[ix]);
      psi0[ix] = psi[ix];

    }

    // compute phase
    for (int ix = 0; ix < Nx; ix++) {
        phase[ix] = atan2(imag(psi[ix]), real(psi[ix]));
    }

    // imposing boundary conditions
    phase[0] = phase[Nx-2];
    phase[Nx-1] = phase[1];

    // calculate fluctuations
    if (it > thld * Nt) {

        for (int ix = 1; ix < Nx-1; ix++) {

            double f_tot = ((abs(psi[ix]) * abs(psi[ix])-n0)/n0);
            double f_kin =  - eta * (phase[ix+1] - 2 * phase[ix] + phase[ix-1]) / dxx;
            double f_noi = 2 * eta * sqrt(D/n0) * (w_x[ix] - w_x_past[ix]) / dt;

            fluct_total[ix] = fluct_total[ix] + f_tot * f_tot;
            fluct_kinetic[ix] = fluct_kinetic[ix] + f_kin * f_kin;
            fluct_noise[ix] = fluct_noise[ix] + f_noi * f_noi;
            }
    }
    // save time step to file
    if (it % freq == 0) {
      for (int ix = 0; ix < Nx; ix++) {

	double xx = x0 + ix * dx;

        ofile << xx << "\t" << abs(psi[ix]) * abs(psi[ix]) << "\t" << endl;
        ophase << tt << "\t" << xx << "\t" << phase[ix] << endl;
        omodule << tt << "\t" << xx << "\t" << abs(psi[ix]) * abs(psi[ix]) << endl;
      }

      ofile << endl << endl;
      ophase << endl;
      omodule << endl;

      n_stats(w_x, w_y);

      // save statistics to file
      no[it] = norm(psi);
      ostats << tt << "\t" << no[it] << "\t" << real(n_mean) << "\t" << imag(n_mean) << "\t" << real(n_var) << "\t"
             << imag(n_var) << "\t" << cgle_var << "\t" << endl;
    }

    no[it] = norm(psi);

    if (it % (Nt/100) == 0)
      cout << setprecision(2) << "[S = " << to_string(is) << "/" << Ns-1 << "|t = " << (double) it/Nt << "]: "
           << "\t\t N = " << no[it] << endl;
  }

/****************************************************************/
/* End dynamics */
/****************************************************************/
  // close handles
  ofile.close();
  ostats.close();
  ophase.close();
  omodule.close();

  /* Calculate correlation function at time t = T and save result to file */
  g1d = corr(psi);

  for (int id = 0; id < Nx-1; id++) {
    double dd = (-Nx/2 + id) * dx;
    ocorr << dd << " " << real(g1d[id]) << " " << imag(g1d[id]) << endl;
  }

  /* Save fluctuation information */
  ofluct << "# x\t total\t kinetic_term\t noise_term\t noise_term_theory" << endl;
  for (int ix = 0; ix < Nx; ix++) {
    double xx = x0 + ix * dx;

    ofluct << xx << "\t" << fluct_total[ix] / ((1-thld) * Nt) << "\t" << beta * fluct_kinetic[ix] / ((1-thld) * Nt) << "\t"
           << beta * fluct_noise[ix] / ((1-thld) * Nt) << "\t" << cgle_var * (dt/dx) << endl;
  }

  time_e = omp_get_wtime();
  double time = time_e - time_s;
  cout << endl << "Simulation ended (t_sim = " << (int) (time / 60)
       << " min " << (int) (time) % 60 << " sec)" << endl;

  ocorr.close();
  ofluct.close();

  /* Save simulation information to file*/
    oinfo << "# Physical parameters: " << endl
	  << "g =\t" << g << endl
	  << "R =\t" << R << endl
	  << "gam_r =\t" << gam_r << endl
	  << "g_r =\t" << g_r << endl
	  << "gam_c =\t" << gam_c << endl
      << "beta =\t" << beta << endl
	  << "P_th =\t" << P_th << endl
	  << "P =\t" << P << endl
	  << "W =\t" << W << endl
	  << "phi_0 = \t" << phi_0 << endl
	  << "n_0 = \t" << n0 << endl
	  << "n_s = \t" << ns << endl
	  << "D = \t" << D << endl
	  << "eta = \t" << eta << endl
	  << "cgle_var = \t" << cgle_var << endl
	  << "L =\t" << L << endl
	  << "T =\t" << T << endl
	  << "freq =\t" << freq << endl
	  << "dx =\t" << dx << endl
	  << "dt =\t" << dt << endl
	  << "tsim =\t" << (int) (time / 60) << " min "
	  << (int) (time) % 60 << " sec";
}
/**************************************************************************/

  double end_time = omp_get_wtime();
  double elapsed_time = end_time - start_time;
  cout << endl << "Total time elapsed: " << (int) (elapsed_time / 60)
       << " min " << (int) (elapsed_time) % 60 << " sec " << endl;

  return 0;
/**************************************************************************/
}
