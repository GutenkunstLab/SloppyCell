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

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *delta, double *pd, double *cj_ptr, double *h_ptr, double *wt, double *constants, int *intpar);

void dres_da_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dA1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dB1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dc1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dr1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_ds1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dr_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD0_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dA2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dB2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dc2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dr2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_ds3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dA3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dB3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dc3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dr3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_ds5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dT3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dT4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dT2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dT1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dparameter_0000073_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dparameter_0000072_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_ds4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_ds6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_ds2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dD10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dL10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar);

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar);

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar);

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar);

void integrate_stochatic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory);

void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar);

double root(double n,double x){
return pow(x, 1.0/n);
}

double root_0(double n,double x){
return -(log(x)*pow(x, 1.0/n)*1.0/pow(n, 2.0));
}

double root_1(double n,double x){
return pow(x, 1.0/n - 1.0)/n;
}

double cot(double x){
return 1.0/tan(x);
}

double cot_0(double x){
return -(1.0/(pow(cos(x), 2.0)*pow(tan(x), 2.0)));
}

double arccot(double x){
return atan(1.0/x);
}

double arccot_0(double x){
return -(1.0/pow(x, 2.0)/(pow(1.0/x, 2.0) + 1.0));
}

double coth(double x){
return 1.0/tanh(x);
}

double coth_0(double x){
return -(1.0/(pow(cosh(x), 2.0)*pow(tanh(x), 2.0)));
}

double csc(double x){
return 1.0/sin(x);
}

double csc_0(double x){
return -(cos(x)/pow(sin(x), 2.0));
}

double arccsc(double x){
return asin(1.0/x);
}

double arccsc_0(double x){
return -(1.0/pow(x, 2.0)/sqrt(1.0 - pow(1.0/x, 2.0)));
}

double csch(double x){
return 1.0/sinh(x);
}

double csch_0(double x){
return -(cosh(x)/pow(sinh(x), 2.0));
}

double sec(double x){
return 1.0/cos(x);
}

double sec_0(double x){
return sin(x)/pow(cos(x), 2.0);
}

double arcsec(double x){
return acos(1.0/x);
}

double arcsec_0(double x){
return -(1.0/pow(x, 2.0)/sqrt(1.0 - pow(1.0/x, 2.0)));
}

double sech(double x){
return 1.0/cosh(x);
}

double sech_0(double x){
return -(sinh(x)/pow(cosh(x), 2.0));
}

void res_function_(double *time_ptr, double *dynamicVars, double *yprime, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


residual[0] = compartment_0000003*CCn*T4/(k4 + CCn) + compartment_0000003*(Clkc*v3*species_0000012 - parameter_0000073*CCc) - compartment_0000002*CCc*T3/(k3 + CCc) - Drosophilia*CCc*D0 - Drosophilia*CCc*D9/(CCc + L9) - yprime[0];
residual[1] = compartment_0000002*CCc*T3/(k3 + CCc) - compartment_0000003*CCn*T4/(k4 + CCn) - Drosophilia*CCn*D0 - Drosophilia*CCn*D10/(CCn + L10) - yprime[1];
residual[2] = compartment_0000003*Clkm*s6 - compartment_0000003*(Clkc*v3*species_0000012 - parameter_0000073*CCc) - Drosophilia*Clkc*D0 - Drosophilia*Clkc*D8/(Clkc + L8) - yprime[2];
residual[3] = compartment_0000003*(c3 + (B3 + pow(PTn/A3, a))*s5/(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0)) - Drosophilia*Clkm*D0 - Drosophilia*Clkm*D7/(Clkm + L7) - yprime[3];
residual[4] = compartment_0000003*s2*Perm - compartment_0000003*(Perc*Timc*v1 - parameter_0000072*PTc) - Drosophilia*D0*Perc - Drosophilia*D2*species_0000013*Perc/(L2 + Perc) - yprime[4];
residual[5] = compartment_0000003*(c1 + (B1 + pow(CCn/A1, a))*s1/(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0)) - compartment_0000003*D0*Perm - Drosophilia*D1*Perm/(L1 + Perm) - yprime[5];
residual[6] = compartment_0000003*PTn*T2/(k2 + PTn) + compartment_0000003*(Perc*Timc*v1 - parameter_0000072*PTc) - compartment_0000002*PTc*T1/(k1 + PTc) - Drosophilia*D0*PTc - Drosophilia*D5*PTc/(L5 + PTc) - yprime[6];
residual[7] = compartment_0000002*PTc*T1/(k1 + PTc) - compartment_0000003*PTn*T2/(k2 + PTn) - Drosophilia*D0*PTn - Drosophilia*D6*PTn/(L6 + PTn) - yprime[7];
residual[8] = compartment_0000003*s4*Timm - compartment_0000003*(Perc*Timc*v1 - parameter_0000072*PTc) - Drosophilia*D0*Timc - Drosophilia*D4*Timc/(L4 + Timc) - yprime[8];
residual[9] = compartment_0000003*(c2 + (B2 + pow(CCn/A2, a))*s3/(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0)) - Drosophilia*D0*Timm - Drosophilia*D3*Timm/(L3 + Timm) - yprime[9];
}

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;

double CCc, CCn, Clkc, Clkm, Perc, Perm, PTc, PTn, Timc, Timm;
double Drosophilia, compartment_0000002, compartment_0000003, EmptySet, species_0000012, species_0000013, a, A1, B1, c1, r1, s1, r, D0, A2, B2, c2, r2, s3, A3, B3, c3, r3, s5, k3, T3, k4, T4, k2, T2, k1, T1, v3, parameter_0000073, v1, parameter_0000072, s4, s6, s2, D1, L1, D2, L2, D3, L3, D4, L4, D5, L5, D6, L6, D7, L7, D8, L8, D9, L9, D10, L10;

Drosophilia = constants[0];
compartment_0000002 = constants[1];
compartment_0000003 = constants[2];
EmptySet = constants[3];
species_0000012 = constants[4];
species_0000013 = constants[5];
a = constants[6];
A1 = constants[7];
B1 = constants[8];
c1 = constants[9];
r1 = constants[10];
s1 = constants[11];
r = constants[12];
D0 = constants[13];
A2 = constants[14];
B2 = constants[15];
c2 = constants[16];
r2 = constants[17];
s3 = constants[18];
A3 = constants[19];
B3 = constants[20];
c3 = constants[21];
r3 = constants[22];
s5 = constants[23];
k3 = constants[24];
T3 = constants[25];
k4 = constants[26];
T4 = constants[27];
k2 = constants[28];
T2 = constants[29];
k1 = constants[30];
T1 = constants[31];
v3 = constants[32];
parameter_0000073 = constants[33];
v1 = constants[34];
parameter_0000072 = constants[35];
s4 = constants[36];
s6 = constants[37];
s2 = constants[38];
D1 = constants[39];
L1 = constants[40];
D2 = constants[41];
L2 = constants[42];
D3 = constants[43];
L3 = constants[44];
D4 = constants[45];
L4 = constants[46];
D5 = constants[47];
L5 = constants[48];
D6 = constants[49];
L6 = constants[50];
D7 = constants[51];
L7 = constants[52];
D8 = constants[53];
L8 = constants[54];
D9 = constants[55];
L9 = constants[56];
D10 = constants[57];
L10 = constants[58];

CCc = dynamicVars[0];
CCn = dynamicVars[1];
Clkc = dynamicVars[2];
Clkm = dynamicVars[3];
Perc = dynamicVars[4];
Perm = dynamicVars[5];
PTc = dynamicVars[6];
PTn = dynamicVars[7];
Timc = dynamicVars[8];
Timm = dynamicVars[9];


}

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = compartment_0000002*CCc*T3/pow(k3 + CCc, 2.0) + Drosophilia*CCc*D9/pow(CCc + L9, 2.0) - compartment_0000003*parameter_0000073 - compartment_0000002*T3/(k3 + CCc) - Drosophilia*D0 - Drosophilia*D9/(CCc + L9);
pd[10] = compartment_0000003*T4/(k4 + CCn) - compartment_0000003*CCn*T4/pow(k4 + CCn, 2.0);
pd[20] = compartment_0000003*v3*species_0000012;
pd[1] = compartment_0000002*T3/(k3 + CCc) - compartment_0000002*CCc*T3/pow(k3 + CCc, 2.0);
pd[11] = compartment_0000003*CCn*T4/pow(k4 + CCn, 2.0) + Drosophilia*CCn*D10/pow(CCn + L10, 2.0) - compartment_0000003*T4/(k4 + CCn) - Drosophilia*D0 - Drosophilia*D10/(CCn + L10);
pd[2] = compartment_0000003*parameter_0000073;
pd[22] = Drosophilia*Clkc*D8/pow(Clkc + L8, 2.0) - compartment_0000003*v3*species_0000012 - Drosophilia*D0 - Drosophilia*D8/(Clkc + L8);
pd[32] = compartment_0000003*s6;
pd[13] = -(compartment_0000003*(B3 + pow(PTn/A3, a))*s5*r*pow(CCn/r3, r - 1.0)/(r3*pow(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0, 2.0)));
pd[33] = Drosophilia*Clkm*D7/pow(Clkm + L7, 2.0) - Drosophilia*D0 - Drosophilia*D7/(Clkm + L7);
pd[73] = compartment_0000003*(s5*a*pow(PTn/A3, a - 1.0)/(A3*(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0)) - (B3 + pow(PTn/A3, a))*s5*a*pow(PTn/A3, a - 1.0)/(A3*pow(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0, 2.0)));
pd[44] = Drosophilia*D2*species_0000013*Perc/pow(L2 + Perc, 2.0) - compartment_0000003*Timc*v1 - Drosophilia*D0 - Drosophilia*D2*species_0000013/(L2 + Perc);
pd[54] = compartment_0000003*s2;
pd[64] = compartment_0000003*parameter_0000072;
pd[84] = -(compartment_0000003*Perc*v1);
pd[15] = compartment_0000003*(s1*a*pow(CCn/A1, a - 1.0)/(A1*(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0)) - (B1 + pow(CCn/A1, a))*s1*a*pow(CCn/A1, a - 1.0)/(A1*pow(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0, 2.0)));
pd[55] = Drosophilia*D1*Perm/pow(L1 + Perm, 2.0) - compartment_0000003*D0 - Drosophilia*D1/(L1 + Perm);
pd[75] = -(compartment_0000003*(B1 + pow(CCn/A1, a))*s1*r*pow(PTn/r1, r - 1.0)/(r1*pow(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0, 2.0)));
pd[46] = compartment_0000003*Timc*v1;
pd[66] = compartment_0000002*PTc*T1/pow(k1 + PTc, 2.0) + Drosophilia*D5*PTc/pow(L5 + PTc, 2.0) - compartment_0000003*parameter_0000072 - compartment_0000002*T1/(k1 + PTc) - Drosophilia*D0 - Drosophilia*D5/(L5 + PTc);
pd[76] = compartment_0000003*T2/(k2 + PTn) - compartment_0000003*PTn*T2/pow(k2 + PTn, 2.0);
pd[86] = compartment_0000003*Perc*v1;
pd[67] = compartment_0000002*T1/(k1 + PTc) - compartment_0000002*PTc*T1/pow(k1 + PTc, 2.0);
pd[77] = compartment_0000003*PTn*T2/pow(k2 + PTn, 2.0) + Drosophilia*D6*PTn/pow(L6 + PTn, 2.0) - compartment_0000003*T2/(k2 + PTn) - Drosophilia*D0 - Drosophilia*D6/(L6 + PTn);
pd[48] = -(compartment_0000003*Timc*v1);
pd[68] = compartment_0000003*parameter_0000072;
pd[88] = Drosophilia*D4*Timc/pow(L4 + Timc, 2.0) - compartment_0000003*Perc*v1 - Drosophilia*D0 - Drosophilia*D4/(L4 + Timc);
pd[98] = compartment_0000003*s4;
pd[19] = compartment_0000003*(s3*a*pow(CCn/A2, a - 1.0)/(A2*(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0)) - (B2 + pow(CCn/A2, a))*s3*a*pow(CCn/A2, a - 1.0)/(A2*pow(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0, 2.0)));
pd[79] = -(compartment_0000003*(B2 + pow(CCn/A2, a))*s3*r*pow(PTn/r2, r - 1.0)/(r2*pow(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0, 2.0)));
pd[99] = Drosophilia*D3*Timm/pow(L3 + Timm, 2.0) - Drosophilia*D0 - Drosophilia*D3/(L3 + Timm);
}

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = -1;
pd[11] = -1;
pd[22] = -1;
pd[33] = -1;
pd[44] = -1;
pd[55] = -1;
pd[66] = -1;
pd[77] = -1;
pd[88] = -1;
pd[99] = -1;
}

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *delta, double *pd, double *cj_ptr, double *h_ptr, double *wt, double *constants, int *intpar){
double cj = *cj_ptr;

double local_dres_dcdot[10*10] = {0};
int ii;

dres_dc_function_(time_ptr, dynamicVars, yprime, constants, pd);

dres_dcdot_function_(time_ptr, dynamicVars, yprime, constants, local_dres_dcdot);

for(ii=0; ii < 100; ii++){
  pd[ii] += cj*local_dres_dcdot[ii];}
}

void dres_da_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = compartment_0000003*(s5*log(PTn/A3)*pow(PTn/A3, a)/(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0) - (B3 + pow(PTn/A3, a))*s5*log(PTn/A3)*pow(PTn/A3, a)/pow(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0, 2.0));
pd[5] = compartment_0000003*(s1*log(CCn/A1)*pow(CCn/A1, a)/(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0) - (B1 + pow(CCn/A1, a))*s1*log(CCn/A1)*pow(CCn/A1, a)/pow(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0, 2.0));
pd[9] = compartment_0000003*(s3*log(CCn/A2)*pow(CCn/A2, a)/(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0) - (B2 + pow(CCn/A2, a))*s3*log(CCn/A2)*pow(CCn/A2, a)/pow(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0, 2.0));
}

void dres_dA1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[5] = compartment_0000003*((B1 + pow(CCn/A1, a))*s1*a*pow(CCn/A1, a - 1.0)*CCn/(pow(A1, 2.0)*pow(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0, 2.0)) - s1*a*pow(CCn/A1, a - 1.0)*CCn/(pow(A1, 2.0)*(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0)));
}

void dres_dB1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[5] = compartment_0000003*(s1/(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0) - (B1 + pow(CCn/A1, a))*s1/pow(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0, 2.0));
}

void dres_dc1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[5] = compartment_0000003;
}

void dres_dr1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[5] = compartment_0000003*(B1 + pow(CCn/A1, a))*s1*r*pow(PTn/r1, r - 1.0)*PTn/(pow(r1, 2.0)*pow(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0, 2.0));
}

void dres_ds1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[5] = compartment_0000003*(B1 + pow(CCn/A1, a))/(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0);
}

void dres_dr_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = -(compartment_0000003*(B3 + pow(PTn/A3, a))*s5*log(CCn/r3)*pow(CCn/r3, r)/pow(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0, 2.0));
pd[5] = -(compartment_0000003*(B1 + pow(CCn/A1, a))*s1*log(PTn/r1)*pow(PTn/r1, r)/pow(B1 + pow(CCn/A1, a) + pow(PTn/r1, r) + 1.0, 2.0));
pd[9] = -(compartment_0000003*(B2 + pow(CCn/A2, a))*s3*log(PTn/r2)*pow(PTn/r2, r)/pow(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0, 2.0));
}

void dres_dD0_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = -(Drosophilia*CCc);
pd[1] = -(Drosophilia*CCn);
pd[2] = -(Drosophilia*Clkc);
pd[3] = -(Drosophilia*Clkm);
pd[4] = -(Drosophilia*Perc);
pd[5] = -(compartment_0000003*Perm);
pd[6] = -(Drosophilia*PTc);
pd[7] = -(Drosophilia*PTn);
pd[8] = -(Drosophilia*Timc);
pd[9] = -(Drosophilia*Timm);
}

void dres_dA2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[9] = compartment_0000003*((B2 + pow(CCn/A2, a))*s3*a*pow(CCn/A2, a - 1.0)*CCn/(pow(A2, 2.0)*pow(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0, 2.0)) - s3*a*pow(CCn/A2, a - 1.0)*CCn/(pow(A2, 2.0)*(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0)));
}

void dres_dB2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[9] = compartment_0000003*(s3/(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0) - (B2 + pow(CCn/A2, a))*s3/pow(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0, 2.0));
}

void dres_dc2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[9] = compartment_0000003;
}

void dres_dr2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[9] = compartment_0000003*(B2 + pow(CCn/A2, a))*s3*r*pow(PTn/r2, r - 1.0)*PTn/(pow(r2, 2.0)*pow(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0, 2.0));
}

void dres_ds3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[9] = compartment_0000003*(B2 + pow(CCn/A2, a))/(B2 + pow(CCn/A2, a) + pow(PTn/r2, r) + 1.0);
}

void dres_dA3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = compartment_0000003*((B3 + pow(PTn/A3, a))*s5*a*pow(PTn/A3, a - 1.0)*PTn/(pow(A3, 2.0)*pow(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0, 2.0)) - s5*a*pow(PTn/A3, a - 1.0)*PTn/(pow(A3, 2.0)*(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0)));
}

void dres_dB3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = compartment_0000003*(s5/(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0) - (B3 + pow(PTn/A3, a))*s5/pow(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0, 2.0));
}

void dres_dc3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = compartment_0000003;
}

void dres_dr3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = compartment_0000003*(B3 + pow(PTn/A3, a))*s5*r*pow(CCn/r3, r - 1.0)*CCn/(pow(r3, 2.0)*pow(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0, 2.0));
}

void dres_ds5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = compartment_0000003*(B3 + pow(PTn/A3, a))/(B3 + pow(PTn/A3, a) + pow(CCn/r3, r) + 1.0);
}

void dres_dk3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = compartment_0000002*CCc*T3/pow(k3 + CCc, 2.0);
pd[1] = -(compartment_0000002*CCc*T3/pow(k3 + CCc, 2.0));
}

void dres_dT3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = -(compartment_0000002*CCc/(k3 + CCc));
pd[1] = compartment_0000002*CCc/(k3 + CCc);
}

void dres_dk4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = -(compartment_0000003*CCn*T4/pow(k4 + CCn, 2.0));
pd[1] = compartment_0000003*CCn*T4/pow(k4 + CCn, 2.0);
}

void dres_dT4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = compartment_0000003*CCn/(k4 + CCn);
pd[1] = -(compartment_0000003*CCn/(k4 + CCn));
}

void dres_dk2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[6] = -(compartment_0000003*PTn*T2/pow(k2 + PTn, 2.0));
pd[7] = compartment_0000003*PTn*T2/pow(k2 + PTn, 2.0);
}

void dres_dT2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[6] = compartment_0000003*PTn/(k2 + PTn);
pd[7] = -(compartment_0000003*PTn/(k2 + PTn));
}

void dres_dk1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[6] = compartment_0000002*PTc*T1/pow(k1 + PTc, 2.0);
pd[7] = -(compartment_0000002*PTc*T1/pow(k1 + PTc, 2.0));
}

void dres_dT1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[6] = -(compartment_0000002*PTc/(k1 + PTc));
pd[7] = compartment_0000002*PTc/(k1 + PTc);
}

void dres_dv3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = compartment_0000003*Clkc*species_0000012;
pd[2] = -(compartment_0000003*Clkc*species_0000012);
}

void dres_dparameter_0000073_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = -(compartment_0000003*CCc);
pd[2] = compartment_0000003*CCc;
}

void dres_dv1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[4] = -(compartment_0000003*Perc*Timc);
pd[6] = compartment_0000003*Perc*Timc;
pd[8] = -(compartment_0000003*Perc*Timc);
}

void dres_dparameter_0000072_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[4] = compartment_0000003*PTc;
pd[6] = -(compartment_0000003*PTc);
pd[8] = compartment_0000003*PTc;
}

void dres_ds4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[8] = compartment_0000003*Timm;
}

void dres_ds6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[2] = compartment_0000003*Clkm;
}

void dres_ds2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[4] = compartment_0000003*Perm;
}

void dres_dD1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[5] = -(Drosophilia*Perm/(L1 + Perm));
}

void dres_dL1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[5] = Drosophilia*D1*Perm/pow(L1 + Perm, 2.0);
}

void dres_dD2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[4] = -(Drosophilia*species_0000013*Perc/(L2 + Perc));
}

void dres_dL2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[4] = Drosophilia*D2*species_0000013*Perc/pow(L2 + Perc, 2.0);
}

void dres_dD3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[9] = -(Drosophilia*Timm/(L3 + Timm));
}

void dres_dL3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[9] = Drosophilia*D3*Timm/pow(L3 + Timm, 2.0);
}

void dres_dD4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[8] = -(Drosophilia*Timc/(L4 + Timc));
}

void dres_dL4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[8] = Drosophilia*D4*Timc/pow(L4 + Timc, 2.0);
}

void dres_dD5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[6] = -(Drosophilia*PTc/(L5 + PTc));
}

void dres_dL5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[6] = Drosophilia*D5*PTc/pow(L5 + PTc, 2.0);
}

void dres_dD6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[7] = -(Drosophilia*PTn/(L6 + PTn));
}

void dres_dL6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[7] = Drosophilia*D6*PTn/pow(L6 + PTn, 2.0);
}

void dres_dD7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = -(Drosophilia*Clkm/(Clkm + L7));
}

void dres_dL7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[3] = Drosophilia*Clkm*D7/pow(Clkm + L7, 2.0);
}

void dres_dD8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[2] = -(Drosophilia*Clkc/(Clkc + L8));
}

void dres_dL8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[2] = Drosophilia*Clkc*D8/pow(Clkc + L8, 2.0);
}

void dres_dD9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = -(Drosophilia*CCc/(CCc + L9));
}

void dres_dL9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[0] = Drosophilia*CCc*D9/pow(CCc + L9, 2.0);
}

void dres_dD10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[1] = -(Drosophilia*CCn/(CCn + L10));
}

void dres_dL10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCn = dynamicVars[1];
double Clkc = dynamicVars[2];
double Clkm = dynamicVars[3];
double Perc = dynamicVars[4];
double Perm = dynamicVars[5];
double PTc = dynamicVars[6];
double PTn = dynamicVars[7];
double Timc = dynamicVars[8];
double Timm = dynamicVars[9];


pd[1] = Drosophilia*CCn*D10/pow(CCn + L10, 2.0);
}

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){

int p_index = (int)constants[59];
double constants_only[59];
int jj;
double *dc_dp;
double *dcdot_dp;
double *local_dres_dp;
int ii;
double local_dres_dc[100] = {0};
double local_dres_dcdot[100] = {0};
int row, col;

res_function_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);

for (jj = 0; jj < 59; jj++){
constants_only[jj] = constants[jj];}
dc_dp = &sens_y[10];
dcdot_dp = &sens_yp[10];
local_dres_dp = &sens_res[10];
for(ii = 0; ii < 10; ii++){
local_dres_dp[ii] = 0;}
switch(p_index)
{
case 0 : dres_da_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 1 : dres_dA1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 2 : dres_dB1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 3 : dres_dc1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 4 : dres_dr1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 5 : dres_ds1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 6 : dres_dr_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 7 : dres_dD0_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 8 : dres_dA2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 9 : dres_dB2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 10 : dres_dc2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 11 : dres_dr2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 12 : dres_ds3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 13 : dres_dA3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 14 : dres_dB3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 15 : dres_dc3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 16 : dres_dr3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 17 : dres_ds5_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 18 : dres_dk3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 19 : dres_dT3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 20 : dres_dk4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 21 : dres_dT4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 22 : dres_dk2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 23 : dres_dT2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 24 : dres_dk1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 25 : dres_dT1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 26 : dres_dv3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 27 : dres_dparameter_0000073_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 28 : dres_dv1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 29 : dres_dparameter_0000072_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 30 : dres_ds4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 31 : dres_ds6_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 32 : dres_ds2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 33 : dres_dD1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 34 : dres_dL1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 35 : dres_dD2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 36 : dres_dL2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 37 : dres_dD3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 38 : dres_dL3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 39 : dres_dD4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 40 : dres_dL4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 41 : dres_dD5_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 42 : dres_dL5_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 43 : dres_dD6_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 44 : dres_dL6_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 45 : dres_dD7_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 46 : dres_dL7_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 47 : dres_dD8_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 48 : dres_dL8_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 49 : dres_dD9_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 50 : dres_dL9_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 51 : dres_dD10_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 52 : dres_dL10_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
}
dres_dc_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dc);
for(row = 0; row < 10; row++){
for(col = 0; col < 10; col++){
sens_res[row+10] += local_dres_dc[row + col*10]*dc_dp[col];}}
dres_dcdot_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dcdot);
for(row = 0; row < 10; row++){
for(col = 0; col < 10; col++){
sens_res[row+10] += local_dres_dcdot[row + col*10]*dcdot_dp[col];}}
}

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double dynamicVars[10];
double yprime[10];
int ii;
for(ii = 0; ii < 10; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
res_function_(time_ptr, dynamicVars, yprime, cj_ptr, residual, ires_ptr, constants, ipar);
}

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double dynamicVars[10];
double yprime[10];
int ii;
for(ii = 0; ii < 10; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
root_func_(neq_ptr, time_ptr, dynamicVars, yprime, nrt_ptr, root_devs, constants, ipar);
}

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){
double sens_y[20];
double sens_yp[20];
int ii;
for(ii = 0; ii < 10; ii++){
sens_y[ii] = max(exp(sens_y_log[ii]), DBL_MIN);
sens_yp[ii] = sens_yp_log[ii] * sens_y[ii];}
for(ii = 10; ii < 20; ii++){
sens_y[ii] = sens_y_log[ii];
sens_yp[ii] = sens_yp_log[ii];}
sens_rhs_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);
}

void integrate_stochastic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory) {
  int i; /* Temp variables */

  unsigned long seed = *seed_ptr;
  short stch[32][10] = {{0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, -1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, -1},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
                    {-1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {1, -1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, -1, 0, 0},
                    {0, 0, 0, 0, 0, 0, -1, 1, 0, 0},
                    {1, 0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, -1, 0, 1, 0, -1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
                    {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
                    {0, 0, 0, 0, 0, -1, 0, 0, 0, 0},
                    {0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, -1},
                    {0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
                    {0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
                    {0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
                    {0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, -1, 0, 0, 0, 0, 0, 0, 0, 0}};
  short depd[32+1][32] = {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                    {1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
                    {1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
                    {1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
                    {1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                    {1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                    {1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                    {1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

  double time = *time_ptr;
  double sd = (*rmsd_ptr)*(*rmsd_ptr)*10.;
  double stop_time = *stop_time_ptr;
  double dt=0.0;

  double dv0[10];
  int rxnInd = 32;
  double propensity, selection, props[32], av[0];
    double _sd = 0.0;
  if (*reseed) {init_genrand(seed);}

  for (i=0;i<10;i++) {dv0[i]=dv[i];}

  while (time < stop_time) {

    if (depd[rxnInd][0]) {props[0]=cv[2]*(cv[9] + (cv[8] + pow(dv[1]/cv[7], cv[6]))*cv[11]/(1.0 + cv[8] + pow(dv[1]/cv[7], cv[6]) + pow(dv[7]/cv[10], cv[12])));}
    if (depd[rxnInd][1]) {props[1]=cv[2]*cv[13]*dv[5];}
    if (depd[rxnInd][2]) {props[2]=cv[2]*(cv[16] + (cv[15] + pow(dv[1]/cv[14], cv[6]))*cv[18]/(1.0 + cv[15] + pow(dv[1]/cv[14], cv[6]) + pow(dv[7]/cv[17], cv[12])));}
    if (depd[rxnInd][3]) {props[3]=cv[0]*cv[13]*dv[9];}
    if (depd[rxnInd][4]) {props[4]=cv[2]*(cv[21] + (cv[20] + pow(dv[7]/cv[19], cv[6]))*cv[23]/(1.0 + cv[20] + pow(dv[7]/cv[19], cv[6]) + pow(dv[1]/cv[22], cv[12])));}
    if (depd[rxnInd][5]) {props[5]=cv[0]*dv[3]*cv[13];}
    if (depd[rxnInd][6]) {props[6]=cv[1]*dv[0]*cv[25]/(cv[24] + dv[0]);}
    if (depd[rxnInd][7]) {props[7]=cv[2]*dv[1]*cv[27]/(cv[26] + dv[1]);}
    if (depd[rxnInd][8]) {props[8]=cv[2]*dv[7]*cv[29]/(cv[28] + dv[7]);}
    if (depd[rxnInd][9]) {props[9]=cv[1]*dv[6]*cv[31]/(cv[30] + dv[6]);}
    if (depd[rxnInd][10]) {props[10]=cv[2]*(dv[2]*cv[32]*cv[4] - cv[33]*dv[0]);}
    if (depd[rxnInd][11]) {props[11]=cv[2]*(dv[4]*dv[8]*cv[34] - cv[35]*dv[6]);}
    if (depd[rxnInd][12]) {props[12]=cv[2]*cv[36]*dv[9];}
    if (depd[rxnInd][13]) {props[13]=cv[2]*dv[3]*cv[37];}
    if (depd[rxnInd][14]) {props[14]=cv[2]*cv[38]*dv[5];}
    if (depd[rxnInd][15]) {props[15]=cv[0]*cv[13]*dv[4];}
    if (depd[rxnInd][16]) {props[16]=cv[0]*cv[13]*dv[6];}
    if (depd[rxnInd][17]) {props[17]=cv[0]*cv[13]*dv[7];}
    if (depd[rxnInd][18]) {props[18]=cv[0]*dv[0]*cv[13];}
    if (depd[rxnInd][19]) {props[19]=cv[0]*dv[2]*cv[13];}
    if (depd[rxnInd][20]) {props[20]=cv[0]*dv[1]*cv[13];}
    if (depd[rxnInd][21]) {props[21]=cv[0]*cv[13]*dv[8];}
    if (depd[rxnInd][22]) {props[22]=cv[0]*cv[39]*dv[5]/(cv[40] + dv[5]);}
    if (depd[rxnInd][23]) {props[23]=cv[0]*cv[41]*cv[5]*dv[4]/(cv[42] + dv[4]);}
    if (depd[rxnInd][24]) {props[24]=cv[0]*cv[43]*dv[9]/(cv[44] + dv[9]);}
    if (depd[rxnInd][25]) {props[25]=cv[0]*cv[45]*dv[8]/(cv[46] + dv[8]);}
    if (depd[rxnInd][26]) {props[26]=cv[0]*cv[47]*dv[6]/(cv[48] + dv[6]);}
    if (depd[rxnInd][27]) {props[27]=cv[0]*cv[49]*dv[7]/(cv[50] + dv[7]);}
    if (depd[rxnInd][28]) {props[28]=cv[0]*dv[3]*cv[51]/(dv[3] + cv[52]);}
    if (depd[rxnInd][29]) {props[29]=cv[0]*dv[2]*cv[53]/(dv[2] + cv[54]);}
    if (depd[rxnInd][30]) {props[30]=cv[0]*dv[0]*cv[55]/(dv[0] + cv[56]);}
    if (depd[rxnInd][31]) {props[31]=cv[0]*dv[1]*cv[57]/(dv[1] + cv[58]);}

    propensity = 0.0;
    for (i=0;i<32;i++) {
      propensity += props[i];}
   if (propensity<=0.0) {
      dt = stop_time-time;
      time = stop_time;
      break;
   }

    dt = -log(1.0-genrand_real32())/propensity;
    time += dt;

    selection = propensity * genrand_real32();

    for (rxnInd=0; rxnInd<32; rxnInd++) {
      if (selection < props[rxnInd]) {break;}
      else {selection -= props[rxnInd];}}

    for (i=0;i<10;i++) {dv[i]+=stch[rxnInd][i];}

    for (i=0;i<10;i++) {
        _sd += (dv0[i]-dv[i])*(dv0[i]-dv[i]);
    }
    if (_sd > sd) {break;}
  }

  for (i=0;i<10;i++) {
    trajectory[i]=(double)dv[i];
    if (time > stop_time) {
      trajectory[i] += (double)stch[rxnInd][i]/dt*(stop_time - time);
    }
  }
  if (time>stop_time) {(*stop_time_ptr) = stop_time;}
  else {(*stop_time_ptr) = time;}

  (*time_ptr) = time;
}


void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double time = *time_ptr;

double Drosophilia = constants[0];
double compartment_0000002 = constants[1];
double compartment_0000003 = constants[2];
double EmptySet = constants[3];
double species_0000012 = constants[4];
double species_0000013 = constants[5];
double a = constants[6];
double A1 = constants[7];
double B1 = constants[8];
double c1 = constants[9];
double r1 = constants[10];
double s1 = constants[11];
double r = constants[12];
double D0 = constants[13];
double A2 = constants[14];
double B2 = constants[15];
double c2 = constants[16];
double r2 = constants[17];
double s3 = constants[18];
double A3 = constants[19];
double B3 = constants[20];
double c3 = constants[21];
double r3 = constants[22];
double s5 = constants[23];
double k3 = constants[24];
double T3 = constants[25];
double k4 = constants[26];
double T4 = constants[27];
double k2 = constants[28];
double T2 = constants[29];
double k1 = constants[30];
double T1 = constants[31];
double v3 = constants[32];
double parameter_0000073 = constants[33];
double v1 = constants[34];
double parameter_0000072 = constants[35];
double s4 = constants[36];
double s6 = constants[37];
double s2 = constants[38];
double D1 = constants[39];
double L1 = constants[40];
double D2 = constants[41];
double L2 = constants[42];
double D3 = constants[43];
double L3 = constants[44];
double D4 = constants[45];
double L4 = constants[46];
double D5 = constants[47];
double L5 = constants[48];
double D6 = constants[49];
double L6 = constants[50];
double D7 = constants[51];
double L7 = constants[52];
double D8 = constants[53];
double L8 = constants[54];
double D9 = constants[55];
double L9 = constants[56];
double D10 = constants[57];
double L10 = constants[58];

double CCc = dynamicVars[0];
double CCc_deriv_wrt_time = yprime[0];
double CCn = dynamicVars[1];
double CCn_deriv_wrt_time = yprime[1];
double Clkc = dynamicVars[2];
double Clkc_deriv_wrt_time = yprime[2];
double Clkm = dynamicVars[3];
double Clkm_deriv_wrt_time = yprime[3];
double Perc = dynamicVars[4];
double Perc_deriv_wrt_time = yprime[4];
double Perm = dynamicVars[5];
double Perm_deriv_wrt_time = yprime[5];
double PTc = dynamicVars[6];
double PTc_deriv_wrt_time = yprime[6];
double PTn = dynamicVars[7];
double PTn_deriv_wrt_time = yprime[7];
double Timc = dynamicVars[8];
double Timc_deriv_wrt_time = yprime[8];
double Timm = dynamicVars[9];
double Timm_deriv_wrt_time = yprime[9];


}
