#include <math.h>
#include <stdio.h>
#include <float.h>
#include "mtrand.h"
#define exponentiale M_E
#define pi M_PI
double max(double a, double b){
return a > b ? a : b;}
double min(double a, double b){
return a < b ? a : b;}
double root(double n,double x);

double root_0(double n,double x);

double root_1(double n,double x);

double cot(double x);

double cot_0(double x);

double arccot(double x);

double arccot_0(double x);

double coth(double x);

double coth_0(double x);

double csc(double x);

double csc_0(double x);

double arccsc(double x);

double arccsc_0(double x);

double csch(double x);

double csch_0(double x);

double sec(double x);

double sec_0(double x);

double arcsec(double x);

double arcsec_0(double x);

double sech(double x);

double sech_0(double x);

void res_function_(double *time_ptr, double *dynamicVars, double *yprime, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar);

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res);

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual);

void integrate_stochatic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory);

void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar);

double root(double n,double x){
return pow(x, ((1.0)/(n)));
}

double root_0(double n,double x){
return (-(((log(x))*(pow(x, ((1.0)/(n)))))*((1.0)/(pow(n, 2.0)))));
}

double root_1(double n,double x){
return ((pow(x, (((1.0)/(n)) - (1.0))))/(n));
}

double cot(double x){
return ((1.0)/(tan(x)));
}

double cot_0(double x){
return (-((1.0)/((pow(cos(x), 2.0))*(pow(tan(x), 2.0)))));
}

double arccot(double x){
return atan(((1.0)/(x)));
}

double arccot_0(double x){
return (-(((1.0)/(pow(x, 2.0)))/((pow(((1.0)/(x)), 2.0)) + (1.0))));
}

double coth(double x){
return ((1.0)/(tanh(x)));
}

double coth_0(double x){
return (-((1.0)/((pow(cosh(x), 2.0))*(pow(tanh(x), 2.0)))));
}

double csc(double x){
return ((1.0)/(sin(x)));
}

double csc_0(double x){
return (-((cos(x))/(pow(sin(x), 2.0))));
}

double arccsc(double x){
return asin(((1.0)/(x)));
}

double arccsc_0(double x){
return (-(((1.0)/(pow(x, 2.0)))/(sqrt(((1.0) - (pow(((1.0)/(x)), 2.0)))))));
}

double csch(double x){
return ((1.0)/(sinh(x)));
}

double csch_0(double x){
return (-((cosh(x))/(pow(sinh(x), 2.0))));
}

double sec(double x){
return ((1.0)/(cos(x)));
}

double sec_0(double x){
return ((sin(x))/(pow(cos(x), 2.0)));
}

double arcsec(double x){
return acos(((1.0)/(x)));
}

double arcsec_0(double x){
return (-(((1.0)/(pow(x, 2.0)))/(sqrt(((1.0) - (pow(((1.0)/(x)), 2.0)))))));
}

double sech(double x){
return ((1.0)/(cosh(x)));
}

double sech_0(double x){
return (-((sinh(x))/(pow(cosh(x), 2.0))));
}

void res_function_(double *time_ptr, double *dynamicVars, double *yprime, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];
double dummy_par = constants[4];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

residual[0] = (-((k1)*(X0))) - yprime[0];
residual[1] = ((k2)*(S2)) - yprime[1];
residual[2] = (((k1)*(X0)) - ((k2)*(S2))) - yprime[2];
residual[3] = -yprime[3];
}

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];
double dummy_par = constants[4];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

alg_derivs_res[0] = (((0.0) + ((-(1.0))*(yp[(, <, _, a, s, t, ., I, n, d, e, x,  , o, b, j, e, c, t,  , a, t,  , 0, x, 7, f, 3, b, 0, 7, d, 2, 4, 1, 6, 0, >, )]))) + (((Keq) + (1.0))*(yp[(, <, _, a, s, t, ., I, n, d, e, x,  , o, b, j, e, c, t,  , a, t,  , 0, x, 7, f, 3, b, 0, 7, d, 2, 4, 5, 2, 0, >, )])));
}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;

double X0, X1, T, S1;
double cell, Keq, k1, k2, dummy_par;
double S2;
dynamicVars[3] = alg_vals[0];

cell = constants[0];
Keq = constants[1];
k1 = constants[2];
k2 = constants[3];
dummy_par = constants[4];

X0 = dynamicVars[0];
X1 = dynamicVars[1];
T = dynamicVars[2];
S1 = dynamicVars[3];

S2 = ((Keq)*(S1));

residual[0] = (((S2) + (S1)) - (T));
}

void integrate_stochastic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory) {
return;}

void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];
double dummy_par = constants[4];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

}
