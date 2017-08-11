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

void dres_dV1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKi_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dKK10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

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

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


residual[0] = V2*MKKK_P/(KK2 + MKKK_P) - V1*MKKK/((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK)) - yprime[0];
residual[1] = V1*MKKK/((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK)) - V2*MKKK_P/(KK2 + MKKK_P) - yprime[1];
residual[2] = V6*MKK_P/(KK6 + MKK_P) - k3*MKKK_P*MKK/(KK3 + MKK) - yprime[2];
residual[3] = k3*MKKK_P*MKK/(KK3 + MKK) + V5*MKK_PP/(KK5 + MKK_PP) - k4*MKKK_P*MKK_P/(KK4 + MKK_P) - V6*MKK_P/(KK6 + MKK_P) - yprime[3];
residual[4] = k4*MKKK_P*MKK_P/(KK4 + MKK_P) - V5*MKK_PP/(KK5 + MKK_PP) - yprime[4];
residual[5] = V10*MAPK_P/(KK10 + MAPK_P) - k7*MKK_PP*MAPK/(KK7 + MAPK) - yprime[5];
residual[6] = k7*MKK_PP*MAPK/(KK7 + MAPK) + V9*MAPK_PP/(KK9 + MAPK_PP) - k8*MKK_PP*MAPK_P/(KK8 + MAPK_P) - V10*MAPK_P/(KK10 + MAPK_P) - yprime[6];
residual[7] = k8*MKK_PP*MAPK_P/(KK8 + MAPK_P) - V9*MAPK_PP/(KK9 + MAPK_PP) - yprime[7];
}

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;

double MKKK, MKKK_P, MKK, MKK_P, MKK_PP, MAPK, MAPK_P, MAPK_PP;
double uVol, V1, Ki, n, K1, V2, KK2, k3, KK3, k4, KK4, V5, KK5, V6, KK6, k7, KK7, k8, KK8, V9, KK9, V10, KK10;

uVol = constants[0];
V1 = constants[1];
Ki = constants[2];
n = constants[3];
K1 = constants[4];
V2 = constants[5];
KK2 = constants[6];
k3 = constants[7];
KK3 = constants[8];
k4 = constants[9];
KK4 = constants[10];
V5 = constants[11];
KK5 = constants[12];
V6 = constants[13];
KK6 = constants[14];
k7 = constants[15];
KK7 = constants[16];
k8 = constants[17];
KK8 = constants[18];
V9 = constants[19];
KK9 = constants[20];
V10 = constants[21];
KK10 = constants[22];

MKKK = dynamicVars[0];
MKKK_P = dynamicVars[1];
MKK = dynamicVars[2];
MKK_P = dynamicVars[3];
MKK_PP = dynamicVars[4];
MAPK = dynamicVars[5];
MAPK_P = dynamicVars[6];
MAPK_PP = dynamicVars[7];


}

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[0] = V1*MKKK*(pow(MAPK_PP/Ki, n) + 1.0)/pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0) - V1/((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK));
pd[8] = V2/(KK2 + MKKK_P) - V2*MKKK_P/pow(KK2 + MKKK_P, 2.0);
pd[56] = V1*MKKK*(K1 + MKKK)*n*pow(MAPK_PP/Ki, n - 1.0)/(Ki*pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0));
pd[1] = V1/((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK)) - V1*MKKK*(pow(MAPK_PP/Ki, n) + 1.0)/pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0);
pd[9] = V2*MKKK_P/pow(KK2 + MKKK_P, 2.0) - V2/(KK2 + MKKK_P);
pd[57] = -(V1*MKKK*(K1 + MKKK)*n*pow(MAPK_PP/Ki, n - 1.0)/(Ki*pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0)));
pd[10] = -(k3*MKK/(KK3 + MKK));
pd[18] = k3*MKKK_P*MKK/pow(KK3 + MKK, 2.0) - k3*MKKK_P/(KK3 + MKK);
pd[26] = V6/(KK6 + MKK_P) - V6*MKK_P/pow(KK6 + MKK_P, 2.0);
pd[11] = k3*MKK/(KK3 + MKK) - k4*MKK_P/(KK4 + MKK_P);
pd[19] = k3*MKKK_P/(KK3 + MKK) - k3*MKKK_P*MKK/pow(KK3 + MKK, 2.0);
pd[27] = k4*MKKK_P*MKK_P/pow(KK4 + MKK_P, 2.0) + V6*MKK_P/pow(KK6 + MKK_P, 2.0) - k4*MKKK_P/(KK4 + MKK_P) - V6/(KK6 + MKK_P);
pd[35] = V5/(KK5 + MKK_PP) - V5*MKK_PP/pow(KK5 + MKK_PP, 2.0);
pd[12] = k4*MKK_P/(KK4 + MKK_P);
pd[28] = k4*MKKK_P/(KK4 + MKK_P) - k4*MKKK_P*MKK_P/pow(KK4 + MKK_P, 2.0);
pd[36] = V5*MKK_PP/pow(KK5 + MKK_PP, 2.0) - V5/(KK5 + MKK_PP);
pd[37] = -(k7*MAPK/(KK7 + MAPK));
pd[45] = k7*MKK_PP*MAPK/pow(KK7 + MAPK, 2.0) - k7*MKK_PP/(KK7 + MAPK);
pd[53] = V10/(KK10 + MAPK_P) - V10*MAPK_P/pow(KK10 + MAPK_P, 2.0);
pd[38] = k7*MAPK/(KK7 + MAPK) - k8*MAPK_P/(KK8 + MAPK_P);
pd[46] = k7*MKK_PP/(KK7 + MAPK) - k7*MKK_PP*MAPK/pow(KK7 + MAPK, 2.0);
pd[54] = k8*MKK_PP*MAPK_P/pow(KK8 + MAPK_P, 2.0) + V10*MAPK_P/pow(KK10 + MAPK_P, 2.0) - k8*MKK_PP/(KK8 + MAPK_P) - V10/(KK10 + MAPK_P);
pd[62] = V9/(KK9 + MAPK_PP) - V9*MAPK_PP/pow(KK9 + MAPK_PP, 2.0);
pd[39] = k8*MAPK_P/(KK8 + MAPK_P);
pd[55] = k8*MKK_PP/(KK8 + MAPK_P) - k8*MKK_PP*MAPK_P/pow(KK8 + MAPK_P, 2.0);
pd[63] = V9*MAPK_PP/pow(KK9 + MAPK_PP, 2.0) - V9/(KK9 + MAPK_PP);
}

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[0] = -1;
pd[9] = -1;
pd[18] = -1;
pd[27] = -1;
pd[36] = -1;
pd[45] = -1;
pd[54] = -1;
pd[63] = -1;
}

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *delta, double *pd, double *cj_ptr, double *h_ptr, double *wt, double *constants, int *intpar){
double cj = *cj_ptr;

double local_dres_dcdot[8*8] = {0};
int ii;

dres_dc_function_(time_ptr, dynamicVars, yprime, constants, pd);

dres_dcdot_function_(time_ptr, dynamicVars, yprime, constants, local_dres_dcdot);

for(ii=0; ii < 64; ii++){
  pd[ii] += cj*local_dres_dcdot[ii];}
}

void dres_dV1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[0] = -(MKKK/((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK)));
pd[1] = MKKK/((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK));
}

void dres_dKi_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[0] = -(V1*MKKK*(K1 + MKKK)*n*pow(MAPK_PP/Ki, n - 1.0)*MAPK_PP/(pow(Ki, 2.0)*pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0)));
pd[1] = V1*MKKK*(K1 + MKKK)*n*pow(MAPK_PP/Ki, n - 1.0)*MAPK_PP/(pow(Ki, 2.0)*pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0));
}

void dres_dK1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[0] = V1*MKKK*(pow(MAPK_PP/Ki, n) + 1.0)/pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0);
pd[1] = -(V1*MKKK*(pow(MAPK_PP/Ki, n) + 1.0)/pow((pow(MAPK_PP/Ki, n) + 1.0)*(K1 + MKKK), 2.0));
}

void dres_dV2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[0] = MKKK_P/(KK2 + MKKK_P);
pd[1] = -(MKKK_P/(KK2 + MKKK_P));
}

void dres_dKK2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[0] = -(V2*MKKK_P/pow(KK2 + MKKK_P, 2.0));
pd[1] = V2*MKKK_P/pow(KK2 + MKKK_P, 2.0);
}

void dres_dk3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[2] = -(MKKK_P*MKK/(KK3 + MKK));
pd[3] = MKKK_P*MKK/(KK3 + MKK);
}

void dres_dKK3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[2] = k3*MKKK_P*MKK/pow(KK3 + MKK, 2.0);
pd[3] = -(k3*MKKK_P*MKK/pow(KK3 + MKK, 2.0));
}

void dres_dk4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[3] = -(MKKK_P*MKK_P/(KK4 + MKK_P));
pd[4] = MKKK_P*MKK_P/(KK4 + MKK_P);
}

void dres_dKK4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[3] = k4*MKKK_P*MKK_P/pow(KK4 + MKK_P, 2.0);
pd[4] = -(k4*MKKK_P*MKK_P/pow(KK4 + MKK_P, 2.0));
}

void dres_dV5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[3] = MKK_PP/(KK5 + MKK_PP);
pd[4] = -(MKK_PP/(KK5 + MKK_PP));
}

void dres_dKK5_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[3] = -(V5*MKK_PP/pow(KK5 + MKK_PP, 2.0));
pd[4] = V5*MKK_PP/pow(KK5 + MKK_PP, 2.0);
}

void dres_dV6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[2] = MKK_P/(KK6 + MKK_P);
pd[3] = -(MKK_P/(KK6 + MKK_P));
}

void dres_dKK6_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[2] = -(V6*MKK_P/pow(KK6 + MKK_P, 2.0));
pd[3] = V6*MKK_P/pow(KK6 + MKK_P, 2.0);
}

void dres_dk7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[5] = -(MKK_PP*MAPK/(KK7 + MAPK));
pd[6] = MKK_PP*MAPK/(KK7 + MAPK);
}

void dres_dKK7_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[5] = k7*MKK_PP*MAPK/pow(KK7 + MAPK, 2.0);
pd[6] = -(k7*MKK_PP*MAPK/pow(KK7 + MAPK, 2.0));
}

void dres_dk8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[6] = -(MKK_PP*MAPK_P/(KK8 + MAPK_P));
pd[7] = MKK_PP*MAPK_P/(KK8 + MAPK_P);
}

void dres_dKK8_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[6] = k8*MKK_PP*MAPK_P/pow(KK8 + MAPK_P, 2.0);
pd[7] = -(k8*MKK_PP*MAPK_P/pow(KK8 + MAPK_P, 2.0));
}

void dres_dV9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[6] = MAPK_PP/(KK9 + MAPK_PP);
pd[7] = -(MAPK_PP/(KK9 + MAPK_PP));
}

void dres_dKK9_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[6] = -(V9*MAPK_PP/pow(KK9 + MAPK_PP, 2.0));
pd[7] = V9*MAPK_PP/pow(KK9 + MAPK_PP, 2.0);
}

void dres_dV10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[5] = MAPK_P/(KK10 + MAPK_P);
pd[6] = -(MAPK_P/(KK10 + MAPK_P));
}

void dres_dKK10_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_P = dynamicVars[1];
double MKK = dynamicVars[2];
double MKK_P = dynamicVars[3];
double MKK_PP = dynamicVars[4];
double MAPK = dynamicVars[5];
double MAPK_P = dynamicVars[6];
double MAPK_PP = dynamicVars[7];


pd[5] = -(V10*MAPK_P/pow(KK10 + MAPK_P, 2.0));
pd[6] = V10*MAPK_P/pow(KK10 + MAPK_P, 2.0);
}

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){

int p_index = (int)constants[23];
double constants_only[23];
int jj;
double *dc_dp;
double *dcdot_dp;
double *local_dres_dp;
int ii;
double local_dres_dc[64] = {0};
double local_dres_dcdot[64] = {0};
int row, col;

res_function_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);

for (jj = 0; jj < 23; jj++){
constants_only[jj] = constants[jj];}
dc_dp = &sens_y[8];
dcdot_dp = &sens_yp[8];
local_dres_dp = &sens_res[8];
for(ii = 0; ii < 8; ii++){
local_dres_dp[ii] = 0;}
switch(p_index)
{
case 0 : dres_dV1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 1 : dres_dKi_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 2 : dres_dK1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 3 : dres_dV2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 4 : dres_dKK2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 5 : dres_dk3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 6 : dres_dKK3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 7 : dres_dk4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 8 : dres_dKK4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 9 : dres_dV5_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 10 : dres_dKK5_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 11 : dres_dV6_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 12 : dres_dKK6_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 13 : dres_dk7_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 14 : dres_dKK7_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 15 : dres_dk8_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 16 : dres_dKK8_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 17 : dres_dV9_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 18 : dres_dKK9_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 19 : dres_dV10_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 20 : dres_dKK10_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
}
dres_dc_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dc);
for(row = 0; row < 8; row++){
for(col = 0; col < 8; col++){
sens_res[row+8] += local_dres_dc[row + col*8]*dc_dp[col];}}
dres_dcdot_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dcdot);
for(row = 0; row < 8; row++){
for(col = 0; col < 8; col++){
sens_res[row+8] += local_dres_dcdot[row + col*8]*dcdot_dp[col];}}
}

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double dynamicVars[8];
double yprime[8];
int ii;
for(ii = 0; ii < 8; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
res_function_(time_ptr, dynamicVars, yprime, cj_ptr, residual, ires_ptr, constants, ipar);
}

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double dynamicVars[8];
double yprime[8];
int ii;
for(ii = 0; ii < 8; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
root_func_(neq_ptr, time_ptr, dynamicVars, yprime, nrt_ptr, root_devs, constants, ipar);
}

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){
double sens_y[16];
double sens_yp[16];
int ii;
for(ii = 0; ii < 8; ii++){
sens_y[ii] = max(exp(sens_y_log[ii]), DBL_MIN);
sens_yp[ii] = sens_yp_log[ii] * sens_y[ii];}
for(ii = 8; ii < 16; ii++){
sens_y[ii] = sens_y_log[ii];
sens_yp[ii] = sens_yp_log[ii];}
sens_rhs_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);
}

void integrate_stochastic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory) {
  int i; /* Temp variables */

  unsigned long seed = *seed_ptr;
  short stch[10][8] = {{-1, 1, 0, 0, 0, 0, 0, 0},
                    {1, -1, 0, 0, 0, 0, 0, 0},
                    {0, 0, -1, 1, 0, 0, 0, 0},
                    {0, 0, 0, -1, 1, 0, 0, 0},
                    {0, 0, 0, 1, -1, 0, 0, 0},
                    {0, 0, 1, -1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, -1, 1, 0},
                    {0, 0, 0, 0, 0, 0, -1, 1},
                    {0, 0, 0, 0, 0, 0, 1, -1},
                    {0, 0, 0, 0, 0, 1, -1, 0}};
  short depd[10+1][10] = {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0},
                    {1, 1, 1, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, 1, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 1, 1, 1, 1, 1, 0, 0},
                    {0, 0, 0, 1, 1, 1, 1, 1, 0, 0},
                    {0, 0, 1, 1, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 1, 0, 1},
                    {1, 0, 0, 0, 0, 0, 0, 1, 1, 1},
                    {1, 0, 0, 0, 0, 0, 0, 1, 1, 1},
                    {0, 0, 0, 0, 0, 0, 1, 1, 0, 1},
                    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

  double time = *time_ptr;
  double sd = (*rmsd_ptr)*(*rmsd_ptr)*8.;
  double stop_time = *stop_time_ptr;
  double dt=0.0;

  double dv0[8];
  int rxnInd = 10;
  double propensity, selection, props[10], av[0];
    double _sd = 0.0;
  if (*reseed) {init_genrand(seed);}

  for (i=0;i<8;i++) {dv0[i]=dv[i];}

  while (time < stop_time) {

    if (depd[rxnInd][0]) {props[0]=cv[1]*dv[0]/((1.0 + pow(dv[7]/cv[2], cv[3]))*(cv[4] + dv[0]));}
    if (depd[rxnInd][1]) {props[1]=cv[5]*dv[1]/(cv[6] + dv[1]);}
    if (depd[rxnInd][2]) {props[2]=cv[7]*dv[1]*dv[2]/(cv[8] + dv[2]);}
    if (depd[rxnInd][3]) {props[3]=cv[9]*dv[1]*dv[3]/(cv[10] + dv[3]);}
    if (depd[rxnInd][4]) {props[4]=cv[11]*dv[4]/(cv[12] + dv[4]);}
    if (depd[rxnInd][5]) {props[5]=cv[13]*dv[3]/(cv[14] + dv[3]);}
    if (depd[rxnInd][6]) {props[6]=cv[15]*dv[4]*dv[5]/(cv[16] + dv[5]);}
    if (depd[rxnInd][7]) {props[7]=cv[17]*dv[4]*dv[6]/(cv[18] + dv[6]);}
    if (depd[rxnInd][8]) {props[8]=cv[19]*dv[7]/(cv[20] + dv[7]);}
    if (depd[rxnInd][9]) {props[9]=cv[21]*dv[6]/(cv[22] + dv[6]);}

    propensity = 0.0;
    for (i=0;i<10;i++) {
      propensity += props[i];}
   if (propensity<=0.0) {
      dt = stop_time-time;
      time = stop_time;
      break;
   }

    dt = -log(1.0-genrand_real32())/propensity;
    time += dt;

    selection = propensity * genrand_real32();

    for (rxnInd=0; rxnInd<10; rxnInd++) {
      if (selection < props[rxnInd]) {break;}
      else {selection -= props[rxnInd];}}

    for (i=0;i<8;i++) {dv[i]+=stch[rxnInd][i];}

    for (i=0;i<8;i++) {
        _sd += (dv0[i]-dv[i])*(dv0[i]-dv[i]);
    }
    if (_sd > sd) {break;}
  }

  for (i=0;i<8;i++) {
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

double uVol = constants[0];
double V1 = constants[1];
double Ki = constants[2];
double n = constants[3];
double K1 = constants[4];
double V2 = constants[5];
double KK2 = constants[6];
double k3 = constants[7];
double KK3 = constants[8];
double k4 = constants[9];
double KK4 = constants[10];
double V5 = constants[11];
double KK5 = constants[12];
double V6 = constants[13];
double KK6 = constants[14];
double k7 = constants[15];
double KK7 = constants[16];
double k8 = constants[17];
double KK8 = constants[18];
double V9 = constants[19];
double KK9 = constants[20];
double V10 = constants[21];
double KK10 = constants[22];

double MKKK = dynamicVars[0];
double MKKK_deriv_wrt_time = yprime[0];
double MKKK_P = dynamicVars[1];
double MKKK_P_deriv_wrt_time = yprime[1];
double MKK = dynamicVars[2];
double MKK_deriv_wrt_time = yprime[2];
double MKK_P = dynamicVars[3];
double MKK_P_deriv_wrt_time = yprime[3];
double MKK_PP = dynamicVars[4];
double MKK_PP_deriv_wrt_time = yprime[4];
double MAPK = dynamicVars[5];
double MAPK_deriv_wrt_time = yprime[5];
double MAPK_P = dynamicVars[6];
double MAPK_P_deriv_wrt_time = yprime[6];
double MAPK_PP = dynamicVars[7];
double MAPK_PP_deriv_wrt_time = yprime[7];


}
