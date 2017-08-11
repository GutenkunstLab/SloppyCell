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

void dres_dV_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK1_P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_1P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_d_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_dC_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_dN_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_IP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dn_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_IT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

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

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

residual[0] = Cell*V_2P*P1/(K_2P + P1) + Cell*k_sP*Mp - Cell*V_1P*P0/(K1_P + P0) - Cell*k_d*P0 - yprime[0];
residual[1] = Cell*V_2T*T1/(K_2T + T1) + Cell*k_sT*Mt - Cell*V_1T*T0/(K_1T + T0) - Cell*k_d*T0 - yprime[1];
residual[2] = Cell*V_1P*P0/(K1_P + P0) + Cell*V_4P*P2/(K_4P + P2) - Cell*V_2P*P1/(K_2P + P1) - Cell*V_3P*P1/(K_3P + P1) - Cell*k_d*P1 - yprime[2];
residual[3] = Cell*V_1T*T0/(K_1T + T0) + Cell*V_4T*T2/(K_4T + T2) - Cell*V_2T*T1/(K_2T + T1) - Cell*V_3T*T1/(K_3T + T1) - Cell*k_d*T1 - yprime[3];
residual[4] = Cell*V_3P*P1/(K_3P + P1) - Cell*V_4P*P2/(K_4P + P2) - (Cell*k_d*P2 + Cell*V_dP*P2/(K_dP + P2)) - (Cell*k3*P2*T2 - Cell*k4*CC) - yprime[4];
residual[5] = Cell*V_3T*T1/(K_3T + T1) - Cell*V_4T*T2/(K_4T + T2) - (Cell*k_d*T2 + Cell*V_dT*T2/(K_dT + T2)) - (Cell*k3*P2*T2 - Cell*k4*CC) - yprime[5];
residual[6] = Cell*k3*P2*T2 - Cell*k4*CC - (Cell*k1*CC - compartment_0000002*k2*Cn) - Cell*k_dC*CC - yprime[6];
residual[7] = Cell*k1*CC - compartment_0000002*k2*Cn - compartment_0000002*k_dN*Cn - yprime[7];
residual[8] = Cell*v_sP*pow(K_IP, n)/(pow(K_IP, n) + pow(Cn, n)) - (Cell*k_d*Mp + Cell*V_mP*Mp/(K_mP + Mp)) - yprime[8];
residual[9] = Cell*V_sT*pow(K_IT, n)/(pow(K_IT, n) + pow(Cn, n)) - (Cell*k_d*Mt + Cell*V_mT*Mt/(K_mT + Mt)) - yprime[9];
}

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;

double P0, T0, P1, T1, P2, T2, CC, Cn, Mp, Mt;
double Cell, compartment_0000002, V_mT, V_dT, K1_P, V_1P, K_1T, V_1T, K_2P, V_2P, K_2T, V_2T, K_3P, V_3P, K_3T, V_3T, K_4P, V_4P, K_4T, V_4T, k_d, V_dP, K_dP, K_dT, k3, k4, k1, k2, k_dC, k_dN, v_sP, K_IP, n, V_sT, K_IT, k_sP, k_sT, V_mP, K_mP, K_mT;
double Pt, Tt;

Cell = constants[0];
compartment_0000002 = constants[1];
V_mT = constants[2];
V_dT = constants[3];
K1_P = constants[4];
V_1P = constants[5];
K_1T = constants[6];
V_1T = constants[7];
K_2P = constants[8];
V_2P = constants[9];
K_2T = constants[10];
V_2T = constants[11];
K_3P = constants[12];
V_3P = constants[13];
K_3T = constants[14];
V_3T = constants[15];
K_4P = constants[16];
V_4P = constants[17];
K_4T = constants[18];
V_4T = constants[19];
k_d = constants[20];
V_dP = constants[21];
K_dP = constants[22];
K_dT = constants[23];
k3 = constants[24];
k4 = constants[25];
k1 = constants[26];
k2 = constants[27];
k_dC = constants[28];
k_dN = constants[29];
v_sP = constants[30];
K_IP = constants[31];
n = constants[32];
V_sT = constants[33];
K_IT = constants[34];
k_sP = constants[35];
k_sT = constants[36];
V_mP = constants[37];
K_mP = constants[38];
K_mT = constants[39];

P0 = dynamicVars[0];
T0 = dynamicVars[1];
P1 = dynamicVars[2];
T1 = dynamicVars[3];
P2 = dynamicVars[4];
T2 = dynamicVars[5];
CC = dynamicVars[6];
Cn = dynamicVars[7];
Mp = dynamicVars[8];
Mt = dynamicVars[9];

Pt = P0 + P1 + P2 + CC + Cn;
Tt = T0 + T1 + T2 + CC + Cn;

}

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[0] = Cell*V_1P*P0/pow(K1_P + P0, 2.0) - Cell*V_1P/(K1_P + P0) - Cell*k_d;
pd[20] = Cell*V_2P/(K_2P + P1) - Cell*V_2P*P1/pow(K_2P + P1, 2.0);
pd[80] = Cell*k_sP;
pd[11] = Cell*V_1T*T0/pow(K_1T + T0, 2.0) - Cell*V_1T/(K_1T + T0) - Cell*k_d;
pd[31] = Cell*V_2T/(K_2T + T1) - Cell*V_2T*T1/pow(K_2T + T1, 2.0);
pd[91] = Cell*k_sT;
pd[2] = Cell*V_1P/(K1_P + P0) - Cell*V_1P*P0/pow(K1_P + P0, 2.0);
pd[22] = Cell*V_2P*P1/pow(K_2P + P1, 2.0) + Cell*V_3P*P1/pow(K_3P + P1, 2.0) - Cell*V_2P/(K_2P + P1) - Cell*V_3P/(K_3P + P1) - Cell*k_d;
pd[42] = Cell*V_4P/(K_4P + P2) - Cell*V_4P*P2/pow(K_4P + P2, 2.0);
pd[13] = Cell*V_1T/(K_1T + T0) - Cell*V_1T*T0/pow(K_1T + T0, 2.0);
pd[33] = Cell*V_2T*T1/pow(K_2T + T1, 2.0) + Cell*V_3T*T1/pow(K_3T + T1, 2.0) - Cell*V_2T/(K_2T + T1) - Cell*V_3T/(K_3T + T1) - Cell*k_d;
pd[53] = Cell*V_4T/(K_4T + T2) - Cell*V_4T*T2/pow(K_4T + T2, 2.0);
pd[24] = Cell*V_3P/(K_3P + P1) - Cell*V_3P*P1/pow(K_3P + P1, 2.0);
pd[44] = Cell*V_4P*P2/pow(K_4P + P2, 2.0) + Cell*V_dP*P2/pow(K_dP + P2, 2.0) - Cell*V_4P/(K_4P + P2) - Cell*k_d - Cell*V_dP/(K_dP + P2) - Cell*k3*T2;
pd[54] = -(Cell*k3*P2);
pd[64] = Cell*k4;
pd[35] = Cell*V_3T/(K_3T + T1) - Cell*V_3T*T1/pow(K_3T + T1, 2.0);
pd[45] = -(Cell*k3*T2);
pd[55] = Cell*V_4T*T2/pow(K_4T + T2, 2.0) + Cell*V_dT*T2/pow(K_dT + T2, 2.0) - Cell*V_4T/(K_4T + T2) - Cell*k_d - Cell*V_dT/(K_dT + T2) - Cell*k3*P2;
pd[65] = Cell*k4;
pd[46] = Cell*k3*T2;
pd[56] = Cell*k3*P2;
pd[66] = -(Cell*k4) - Cell*k1 - Cell*k_dC;
pd[76] = compartment_0000002*k2;
pd[67] = Cell*k1;
pd[77] = -(compartment_0000002*k2) - compartment_0000002*k_dN;
pd[78] = -(Cell*v_sP*pow(K_IP, n)*n*pow(Cn, n - 1.0)/pow(pow(K_IP, n) + pow(Cn, n), 2.0));
pd[88] = Cell*V_mP*Mp/pow(K_mP + Mp, 2.0) - Cell*k_d - Cell*V_mP/(K_mP + Mp);
pd[79] = -(Cell*V_sT*pow(K_IT, n)*n*pow(Cn, n - 1.0)/pow(pow(K_IT, n) + pow(Cn, n), 2.0));
pd[99] = Cell*V_mT*Mt/pow(K_mT + Mt, 2.0) - Cell*k_d - Cell*V_mT/(K_mT + Mt);
}

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

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

void dres_dV_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[9] = -(Cell*Mt/(K_mT + Mt));
}

void dres_dV_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[5] = -(Cell*T2/(K_dT + T2));
}

void dres_dK1_P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[0] = Cell*V_1P*P0/pow(K1_P + P0, 2.0);
pd[2] = -(Cell*V_1P*P0/pow(K1_P + P0, 2.0));
}

void dres_dV_1P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[0] = -(Cell*P0/(K1_P + P0));
pd[2] = Cell*P0/(K1_P + P0);
}

void dres_dK_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[1] = Cell*V_1T*T0/pow(K_1T + T0, 2.0);
pd[3] = -(Cell*V_1T*T0/pow(K_1T + T0, 2.0));
}

void dres_dV_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[1] = -(Cell*T0/(K_1T + T0));
pd[3] = Cell*T0/(K_1T + T0);
}

void dres_dK_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[0] = -(Cell*V_2P*P1/pow(K_2P + P1, 2.0));
pd[2] = Cell*V_2P*P1/pow(K_2P + P1, 2.0);
}

void dres_dV_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[0] = Cell*P1/(K_2P + P1);
pd[2] = -(Cell*P1/(K_2P + P1));
}

void dres_dK_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[1] = -(Cell*V_2T*T1/pow(K_2T + T1, 2.0));
pd[3] = Cell*V_2T*T1/pow(K_2T + T1, 2.0);
}

void dres_dV_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[1] = Cell*T1/(K_2T + T1);
pd[3] = -(Cell*T1/(K_2T + T1));
}

void dres_dK_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[2] = Cell*V_3P*P1/pow(K_3P + P1, 2.0);
pd[4] = -(Cell*V_3P*P1/pow(K_3P + P1, 2.0));
}

void dres_dV_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[2] = -(Cell*P1/(K_3P + P1));
pd[4] = Cell*P1/(K_3P + P1);
}

void dres_dK_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[3] = Cell*V_3T*T1/pow(K_3T + T1, 2.0);
pd[5] = -(Cell*V_3T*T1/pow(K_3T + T1, 2.0));
}

void dres_dV_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[3] = -(Cell*T1/(K_3T + T1));
pd[5] = Cell*T1/(K_3T + T1);
}

void dres_dK_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[2] = -(Cell*V_4P*P2/pow(K_4P + P2, 2.0));
pd[4] = Cell*V_4P*P2/pow(K_4P + P2, 2.0);
}

void dres_dV_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[2] = Cell*P2/(K_4P + P2);
pd[4] = -(Cell*P2/(K_4P + P2));
}

void dres_dK_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[3] = -(Cell*V_4T*T2/pow(K_4T + T2, 2.0));
pd[5] = Cell*V_4T*T2/pow(K_4T + T2, 2.0);
}

void dres_dV_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[3] = Cell*T2/(K_4T + T2);
pd[5] = -(Cell*T2/(K_4T + T2));
}

void dres_dk_d_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[0] = -(Cell*P0);
pd[1] = -(Cell*T0);
pd[2] = -(Cell*P1);
pd[3] = -(Cell*T1);
pd[4] = -(Cell*P2);
pd[5] = -(Cell*T2);
pd[8] = -(Cell*Mp);
pd[9] = -(Cell*Mt);
}

void dres_dV_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[4] = -(Cell*P2/(K_dP + P2));
}

void dres_dK_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[4] = Cell*V_dP*P2/pow(K_dP + P2, 2.0);
}

void dres_dK_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[5] = Cell*V_dT*T2/pow(K_dT + T2, 2.0);
}

void dres_dk3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[4] = -(Cell*P2*T2);
pd[5] = -(Cell*P2*T2);
pd[6] = Cell*P2*T2;
}

void dres_dk4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[4] = Cell*CC;
pd[5] = Cell*CC;
pd[6] = -(Cell*CC);
}

void dres_dk1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[6] = -(Cell*CC);
pd[7] = Cell*CC;
}

void dres_dk2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[6] = compartment_0000002*Cn;
pd[7] = -(compartment_0000002*Cn);
}

void dres_dk_dC_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[6] = -(Cell*CC);
}

void dres_dk_dN_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[7] = -(compartment_0000002*Cn);
}

void dres_dv_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[8] = Cell*pow(K_IP, n)/(pow(K_IP, n) + pow(Cn, n));
}

void dres_dK_IP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[8] = Cell*v_sP*n*pow(K_IP, n - 1.0)/(pow(K_IP, n) + pow(Cn, n)) - Cell*v_sP*pow(K_IP, n)*n*pow(K_IP, n - 1.0)/pow(pow(K_IP, n) + pow(Cn, n), 2.0);
}

void dres_dn_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[8] = Cell*v_sP*log(K_IP)*pow(K_IP, n)/(pow(K_IP, n) + pow(Cn, n)) - Cell*v_sP*pow(K_IP, n)*(log(K_IP)*pow(K_IP, n) + log(Cn)*pow(Cn, n))/pow(pow(K_IP, n) + pow(Cn, n), 2.0);
pd[9] = Cell*V_sT*log(K_IT)*pow(K_IT, n)/(pow(K_IT, n) + pow(Cn, n)) - Cell*V_sT*pow(K_IT, n)*(log(K_IT)*pow(K_IT, n) + log(Cn)*pow(Cn, n))/pow(pow(K_IT, n) + pow(Cn, n), 2.0);
}

void dres_dV_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[9] = Cell*pow(K_IT, n)/(pow(K_IT, n) + pow(Cn, n));
}

void dres_dK_IT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[9] = Cell*V_sT*n*pow(K_IT, n - 1.0)/(pow(K_IT, n) + pow(Cn, n)) - Cell*V_sT*pow(K_IT, n)*n*pow(K_IT, n - 1.0)/pow(pow(K_IT, n) + pow(Cn, n), 2.0);
}

void dres_dk_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[0] = Cell*Mp;
}

void dres_dk_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[1] = Cell*Mt;
}

void dres_dV_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[8] = -(Cell*Mp/(K_mP + Mp));
}

void dres_dK_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[8] = Cell*V_mP*Mp/pow(K_mP + Mp, 2.0);
}

void dres_dK_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double T0 = dynamicVars[1];
double P1 = dynamicVars[2];
double T1 = dynamicVars[3];
double P2 = dynamicVars[4];
double T2 = dynamicVars[5];
double CC = dynamicVars[6];
double Cn = dynamicVars[7];
double Mp = dynamicVars[8];
double Mt = dynamicVars[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Tt = T0 + T1 + T2 + CC + Cn;

pd[9] = Cell*V_mT*Mt/pow(K_mT + Mt, 2.0);
}

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){

int p_index = (int)constants[40];
double constants_only[40];
int jj;
double *dc_dp;
double *dcdot_dp;
double *local_dres_dp;
int ii;
double local_dres_dc[100] = {0};
double local_dres_dcdot[100] = {0};
int row, col;

res_function_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);

for (jj = 0; jj < 40; jj++){
constants_only[jj] = constants[jj];}
dc_dp = &sens_y[10];
dcdot_dp = &sens_yp[10];
local_dres_dp = &sens_res[10];
for(ii = 0; ii < 10; ii++){
local_dres_dp[ii] = 0;}
switch(p_index)
{
case 0 : dres_dV_mT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 1 : dres_dV_dT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 2 : dres_dK1_P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 3 : dres_dV_1P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 4 : dres_dK_1T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 5 : dres_dV_1T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 6 : dres_dK_2P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 7 : dres_dV_2P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 8 : dres_dK_2T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 9 : dres_dV_2T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 10 : dres_dK_3P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 11 : dres_dV_3P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 12 : dres_dK_3T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 13 : dres_dV_3T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 14 : dres_dK_4P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 15 : dres_dV_4P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 16 : dres_dK_4T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 17 : dres_dV_4T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 18 : dres_dk_d_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 19 : dres_dV_dP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 20 : dres_dK_dP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 21 : dres_dK_dT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 22 : dres_dk3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 23 : dres_dk4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 24 : dres_dk1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 25 : dres_dk2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 26 : dres_dk_dC_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 27 : dres_dk_dN_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 28 : dres_dv_sP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 29 : dres_dK_IP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 30 : dres_dn_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 31 : dres_dV_sT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 32 : dres_dK_IT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 33 : dres_dk_sP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 34 : dres_dk_sT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 35 : dres_dV_mP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 36 : dres_dK_mP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 37 : dres_dK_mT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
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
  short stch[24][10] = {{-1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, -1, 0, 1, 0, 0, 0, 0, 0, 0},
                    {1, 0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, -1, 0, 0, 0, 0, 0, 0},
                    {0, 0, -1, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, -1, 0, 1, 0, 0, 0, 0},
                    {0, 0, 1, 0, -1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, -1, 0, 0, 0, 0},
                    {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, -1, 0, 0, 0, 0},
                    {0, 0, 0, 0, -1, -1, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, -1, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, -1}};
  short depd[24+1][24] = {{1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1},
                    {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1},
                    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

  double time = *time_ptr;
  double sd = (*rmsd_ptr)*(*rmsd_ptr)*10.;
  double stop_time = *stop_time_ptr;
  double dt=0.0;

  double dv0[10];
  int rxnInd = 24;
  double propensity, selection, props[24], av[2];
    double _sd = 0.0;
  if (*reseed) {init_genrand(seed);}

  for (i=0;i<10;i++) {dv0[i]=dv[i];}

  while (time < stop_time) {
    av[0]=dv[0] + dv[2] + dv[4] + dv[6] + dv[7];
    av[1]=dv[1] + dv[3] + dv[5] + dv[6] + dv[7];

    if (depd[rxnInd][0]) {props[0]=cv[0]*cv[5]*dv[0]/(cv[4] + dv[0]);}
    if (depd[rxnInd][1]) {props[1]=cv[0]*cv[7]*dv[1]/(cv[6] + dv[1]);}
    if (depd[rxnInd][2]) {props[2]=cv[0]*cv[9]*dv[2]/(cv[8] + dv[2]);}
    if (depd[rxnInd][3]) {props[3]=cv[0]*cv[11]*dv[3]/(cv[10] + dv[3]);}
    if (depd[rxnInd][4]) {props[4]=cv[0]*cv[13]*dv[2]/(cv[12] + dv[2]);}
    if (depd[rxnInd][5]) {props[5]=cv[0]*cv[15]*dv[3]/(cv[14] + dv[3]);}
    if (depd[rxnInd][6]) {props[6]=cv[0]*cv[17]*dv[4]/(cv[16] + dv[4]);}
    if (depd[rxnInd][7]) {props[7]=cv[0]*cv[19]*dv[5]/(cv[18] + dv[5]);}
    if (depd[rxnInd][8]) {props[8]=cv[0]*cv[20]*dv[0];}
    if (depd[rxnInd][9]) {props[9]=cv[0]*cv[20]*dv[1];}
    if (depd[rxnInd][10]) {props[10]=cv[0]*cv[20]*dv[2];}
    if (depd[rxnInd][11]) {props[11]=cv[0]*cv[20]*dv[3];}
    if (depd[rxnInd][12]) {props[12]=cv[0]*cv[20]*dv[4] + cv[0]*cv[21]*dv[4]/(cv[22] + dv[4]);}
    if (depd[rxnInd][13]) {props[13]=cv[0]*cv[20]*dv[5] + cv[0]*cv[3]*dv[5]/(cv[23] + dv[5]);}
    if (depd[rxnInd][14]) {props[14]=cv[0]*cv[24]*dv[4]*dv[5] - cv[0]*cv[25]*dv[6];}
    if (depd[rxnInd][15]) {props[15]=cv[0]*cv[26]*dv[6] - cv[1]*cv[27]*dv[7];}
    if (depd[rxnInd][16]) {props[16]=cv[0]*cv[28]*dv[6];}
    if (depd[rxnInd][17]) {props[17]=cv[1]*cv[29]*dv[7];}
    if (depd[rxnInd][18]) {props[18]=cv[0]*cv[30]*pow(cv[31], cv[32])/(pow(cv[31], cv[32]) + pow(dv[7], cv[32]));}
    if (depd[rxnInd][19]) {props[19]=cv[0]*cv[33]*pow(cv[34], cv[32])/(pow(cv[34], cv[32]) + pow(dv[7], cv[32]));}
    if (depd[rxnInd][20]) {props[20]=cv[0]*cv[35]*dv[8];}
    if (depd[rxnInd][21]) {props[21]=cv[0]*cv[36]*dv[9];}
    if (depd[rxnInd][22]) {props[22]=cv[0]*cv[20]*dv[8] + cv[0]*cv[37]*dv[8]/(cv[38] + dv[8]);}
    if (depd[rxnInd][23]) {props[23]=cv[0]*cv[20]*dv[9] + cv[0]*cv[2]*dv[9]/(cv[39] + dv[9]);}

    propensity = 0.0;
    for (i=0;i<24;i++) {
      propensity += props[i];}
   if (propensity<=0.0) {
      dt = stop_time-time;
      time = stop_time;
      break;
   }

    dt = -log(1.0-genrand_real32())/propensity;
    time += dt;

    selection = propensity * genrand_real32();

    for (rxnInd=0; rxnInd<24; rxnInd++) {
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

double Cell = constants[0];
double compartment_0000002 = constants[1];
double V_mT = constants[2];
double V_dT = constants[3];
double K1_P = constants[4];
double V_1P = constants[5];
double K_1T = constants[6];
double V_1T = constants[7];
double K_2P = constants[8];
double V_2P = constants[9];
double K_2T = constants[10];
double V_2T = constants[11];
double K_3P = constants[12];
double V_3P = constants[13];
double K_3T = constants[14];
double V_3T = constants[15];
double K_4P = constants[16];
double V_4P = constants[17];
double K_4T = constants[18];
double V_4T = constants[19];
double k_d = constants[20];
double V_dP = constants[21];
double K_dP = constants[22];
double K_dT = constants[23];
double k3 = constants[24];
double k4 = constants[25];
double k1 = constants[26];
double k2 = constants[27];
double k_dC = constants[28];
double k_dN = constants[29];
double v_sP = constants[30];
double K_IP = constants[31];
double n = constants[32];
double V_sT = constants[33];
double K_IT = constants[34];
double k_sP = constants[35];
double k_sT = constants[36];
double V_mP = constants[37];
double K_mP = constants[38];
double K_mT = constants[39];

double P0 = dynamicVars[0];
double P0_deriv_wrt_time = yprime[0];
double T0 = dynamicVars[1];
double T0_deriv_wrt_time = yprime[1];
double P1 = dynamicVars[2];
double P1_deriv_wrt_time = yprime[2];
double T1 = dynamicVars[3];
double T1_deriv_wrt_time = yprime[3];
double P2 = dynamicVars[4];
double P2_deriv_wrt_time = yprime[4];
double T2 = dynamicVars[5];
double T2_deriv_wrt_time = yprime[5];
double CC = dynamicVars[6];
double CC_deriv_wrt_time = yprime[6];
double Cn = dynamicVars[7];
double Cn_deriv_wrt_time = yprime[7];
double Mp = dynamicVars[8];
double Mp_deriv_wrt_time = yprime[8];
double Mt = dynamicVars[9];
double Mt_deriv_wrt_time = yprime[9];

double Pt = P0 + P1 + P2 + CC + Cn;
double Pt_deriv_wrt_time = P2_deriv_wrt_time + CC_deriv_wrt_time + P0_deriv_wrt_time + P1_deriv_wrt_time + Cn_deriv_wrt_time;
double Tt = T0 + T1 + T2 + CC + Cn;
double Tt_deriv_wrt_time = CC_deriv_wrt_time + T2_deriv_wrt_time + Cn_deriv_wrt_time + T0_deriv_wrt_time + T1_deriv_wrt_time;

}
