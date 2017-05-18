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

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *pd, double *cj_ptr, double *constants, int *intpar);

void dres_dk_2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_IT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_IP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dn_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_1P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_dN_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_1P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_dC_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk_d_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dK_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dV_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

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

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

residual[0] = v_sP*pow(K_IP, n)/(pow(K_IP, n) + pow(C_N, n)) - v_mP*M_P/(K_mP + M_P) - k_d*M_P - yprime[0];
residual[1] = k_sP*M_P + V_2P*P_1/(K_2P + P_1) - V_1P*P_0/(K_1P + P_0) - k_d*P_0 - yprime[1];
residual[2] = V_1P*P_0/(K_1P + P_0) + V_4P*P_2/(K_4P + P_2) - V_2P*P_1/(K_2P + P_1) - V_3P*P_1/(K_3P + P_1) - k_d*P_1 - yprime[2];
residual[3] = V_3P*P_1/(K_3P + P_1) + k_4*C - V_4P*P_2/(K_4P + P_2) - k_3*P_2*T_2 - v_dP*P_2/(K_dP + P_2) - k_d*P_2 - yprime[3];
residual[4] = v_sT*pow(K_IT, n)/(pow(K_IT, n) + pow(C_N, n)) - v_mT*M_T/(K_mT + M_T) - k_d*M_T - yprime[4];
residual[5] = k_sT*M_T + V_2T*T_1/(K_2T + T_1) - V_1T*T_0/(K_1T + T_0) - k_d*T_0 - yprime[5];
residual[6] = V_1T*T_0/(K_1T + T_0) + V_4T*T_2/(K_4T + T_2) - V_2T*T_1/(K_2T + T_1) - V_3T*T_1/(K_3T + T_1) - k_d*T_1 - yprime[6];
residual[7] = V_3T*T_1/(K_3T + T_1) + k_4*C - V_4T*T_2/(K_4T + T_2) - k_3*P_2*T_2 - v_dT*T_2/(K_dT + T_2) - k_d*T_2 - yprime[7];
residual[8] = k_3*P_2*T_2 + k_2*C_N - k_4*C - k_1*C - k_dC*C - yprime[8];
residual[9] = k_1*C - k_2*C_N - k_dN*C_N - yprime[9];
}

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;


double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

}

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[0] = v_mP*M_P/pow(K_mP + M_P, 2.0) - v_mP/(K_mP + M_P) - k_d;
pd[90] = -(v_sP*pow(K_IP, n)*n*pow(C_N, n - 1.0)/pow(pow(K_IP, n) + pow(C_N, n), 2.0));
pd[1] = k_sP;
pd[11] = V_1P*P_0/pow(K_1P + P_0, 2.0) - V_1P/(K_1P + P_0) - k_d;
pd[21] = V_2P/(K_2P + P_1) - V_2P*P_1/pow(K_2P + P_1, 2.0);
pd[12] = V_1P/(K_1P + P_0) - V_1P*P_0/pow(K_1P + P_0, 2.0);
pd[22] = V_2P*P_1/pow(K_2P + P_1, 2.0) + V_3P*P_1/pow(K_3P + P_1, 2.0) - V_2P/(K_2P + P_1) - V_3P/(K_3P + P_1) - k_d;
pd[32] = V_4P/(K_4P + P_2) - V_4P*P_2/pow(K_4P + P_2, 2.0);
pd[23] = V_3P/(K_3P + P_1) - V_3P*P_1/pow(K_3P + P_1, 2.0);
pd[33] = V_4P*P_2/pow(K_4P + P_2, 2.0) + v_dP*P_2/pow(K_dP + P_2, 2.0) - V_4P/(K_4P + P_2) - k_3*T_2 - v_dP/(K_dP + P_2) - k_d;
pd[73] = -(k_3*P_2);
pd[83] = k_4;
pd[44] = v_mT*M_T/pow(K_mT + M_T, 2.0) - v_mT/(K_mT + M_T) - k_d;
pd[94] = -(v_sT*pow(K_IT, n)*n*pow(C_N, n - 1.0)/pow(pow(K_IT, n) + pow(C_N, n), 2.0));
pd[45] = k_sT;
pd[55] = V_1T*T_0/pow(K_1T + T_0, 2.0) - V_1T/(K_1T + T_0) - k_d;
pd[65] = V_2T/(K_2T + T_1) - V_2T*T_1/pow(K_2T + T_1, 2.0);
pd[56] = V_1T/(K_1T + T_0) - V_1T*T_0/pow(K_1T + T_0, 2.0);
pd[66] = V_2T*T_1/pow(K_2T + T_1, 2.0) + V_3T*T_1/pow(K_3T + T_1, 2.0) - V_2T/(K_2T + T_1) - V_3T/(K_3T + T_1) - k_d;
pd[76] = V_4T/(K_4T + T_2) - V_4T*T_2/pow(K_4T + T_2, 2.0);
pd[37] = -(k_3*T_2);
pd[67] = V_3T/(K_3T + T_1) - V_3T*T_1/pow(K_3T + T_1, 2.0);
pd[77] = V_4T*T_2/pow(K_4T + T_2, 2.0) + v_dT*T_2/pow(K_dT + T_2, 2.0) - V_4T/(K_4T + T_2) - k_3*P_2 - v_dT/(K_dT + T_2) - k_d;
pd[87] = k_4;
pd[38] = k_3*T_2;
pd[78] = k_3*P_2;
pd[88] = -k_4 - k_1 - k_dC;
pd[98] = k_2;
pd[89] = k_1;
pd[99] = -k_2 - k_dN;
}

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

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

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *pd, double *cj_ptr, double *constants, int *intpar){
double cj = *cj_ptr;

dres_dc_function_(time_ptr, dynamicVars, yprime, constants, pd);

double local_dres_dcdot[10*10] = {0};
dres_dcdot_function_(time_ptr, dynamicVars, yprime, constants, local_dres_dcdot);

int ii;
for(ii=0; ii < 100; ii++){
  pd[ii] += cj*local_dres_dcdot[ii];}
}

void dres_dk_2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[8] = C_N;
pd[9] = -C_N;
}

void dres_dK_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[5] = -(V_2T*T_1/pow(K_2T + T_1, 2.0));
pd[6] = V_2T*T_1/pow(K_2T + T_1, 2.0);
}

void dres_dK_IT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[4] = v_sT*n*pow(K_IT, n - 1.0)/(pow(K_IT, n) + pow(C_N, n)) - v_sT*pow(K_IT, n)*n*pow(K_IT, n - 1.0)/pow(pow(K_IT, n) + pow(C_N, n), 2.0);
}

void dres_dK_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[7] = v_dT*T_2/pow(K_dT + T_2, 2.0);
}

void dres_dK_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[6] = V_3T*T_1/pow(K_3T + T_1, 2.0);
pd[7] = -(V_3T*T_1/pow(K_3T + T_1, 2.0));
}

void dres_dK_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[1] = -(V_2P*P_1/pow(K_2P + P_1, 2.0));
pd[2] = V_2P*P_1/pow(K_2P + P_1, 2.0);
}

void dres_dK_IP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[0] = v_sP*n*pow(K_IP, n - 1.0)/(pow(K_IP, n) + pow(C_N, n)) - v_sP*pow(K_IP, n)*n*pow(K_IP, n - 1.0)/pow(pow(K_IP, n) + pow(C_N, n), 2.0);
}

void dres_dK_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[3] = v_dP*P_2/pow(K_dP + P_2, 2.0);
}

void dres_dk_4_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[3] = C;
pd[7] = C;
pd[8] = -C;
}

void dres_dK_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[0] = v_mP*M_P/pow(K_mP + M_P, 2.0);
}

void dres_dK_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[4] = v_mT*M_T/pow(K_mT + M_T, 2.0);
}

void dres_dV_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[6] = T_2/(K_4T + T_2);
pd[7] = -(T_2/(K_4T + T_2));
}

void dres_dv_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[0] = pow(K_IP, n)/(pow(K_IP, n) + pow(C_N, n));
}

void dres_dv_mP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[0] = -(M_P/(K_mP + M_P));
}

void dres_dn_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[0] = v_sP*log(K_IP)*pow(K_IP, n)/(pow(K_IP, n) + pow(C_N, n)) - v_sP*pow(K_IP, n)*(log(K_IP)*pow(K_IP, n) + log(C_N)*pow(C_N, n))/pow(pow(K_IP, n) + pow(C_N, n), 2.0);
pd[4] = v_sT*log(K_IT)*pow(K_IT, n)/(pow(K_IT, n) + pow(C_N, n)) - v_sT*pow(K_IT, n)*(log(K_IT)*pow(K_IT, n) + log(C_N)*pow(C_N, n))/pow(pow(K_IT, n) + pow(C_N, n), 2.0);
}

void dres_dv_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[4] = pow(K_IT, n)/(pow(K_IT, n) + pow(C_N, n));
}

void dres_dv_mT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[4] = -(M_T/(K_mT + M_T));
}

void dres_dV_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[5] = -(T_0/(K_1T + T_0));
pd[6] = T_0/(K_1T + T_0);
}

void dres_dV_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[2] = P_2/(K_4P + P_2);
pd[3] = -(P_2/(K_4P + P_2));
}

void dres_dV_1P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[1] = -(P_0/(K_1P + P_0));
pd[2] = P_0/(K_1P + P_0);
}

void dres_dk_1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[8] = -C;
pd[9] = C;
}

void dres_dv_dT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[7] = -(T_2/(K_dT + T_2));
}

void dres_dk_3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[3] = -(P_2*T_2);
pd[7] = -(P_2*T_2);
pd[8] = P_2*T_2;
}

void dres_dk_sT_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[5] = M_T;
}

void dres_dv_dP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[3] = -(P_2/(K_dP + P_2));
}

void dres_dk_sP_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[1] = M_P;
}

void dres_dk_dN_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[9] = -C_N;
}

void dres_dK_1T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[5] = V_1T*T_0/pow(K_1T + T_0, 2.0);
pd[6] = -(V_1T*T_0/pow(K_1T + T_0, 2.0));
}

void dres_dK_1P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[1] = V_1P*P_0/pow(K_1P + P_0, 2.0);
pd[2] = -(V_1P*P_0/pow(K_1P + P_0, 2.0));
}

void dres_dK_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[2] = V_3P*P_1/pow(K_3P + P_1, 2.0);
pd[3] = -(V_3P*P_1/pow(K_3P + P_1, 2.0));
}

void dres_dV_3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[2] = -(P_1/(K_3P + P_1));
pd[3] = P_1/(K_3P + P_1);
}

void dres_dk_dC_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[8] = -C;
}

void dres_dK_4T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[6] = -(V_4T*T_2/pow(K_4T + T_2, 2.0));
pd[7] = V_4T*T_2/pow(K_4T + T_2, 2.0);
}

void dres_dV_2P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[1] = P_1/(K_2P + P_1);
pd[2] = -(P_1/(K_2P + P_1));
}

void dres_dV_3T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[6] = -(T_1/(K_3T + T_1));
pd[7] = T_1/(K_3T + T_1);
}

void dres_dk_d_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[0] = -M_P;
pd[1] = -P_0;
pd[2] = -P_1;
pd[3] = -P_2;
pd[4] = -M_T;
pd[5] = -T_0;
pd[6] = -T_1;
pd[7] = -T_2;
}

void dres_dK_4P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[2] = -(V_4P*P_2/pow(K_4P + P_2, 2.0));
pd[3] = V_4P*P_2/pow(K_4P + P_2, 2.0);
}

void dres_dV_2T_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double P_0 = dynamicVars[1];
double P_1 = dynamicVars[2];
double P_2 = dynamicVars[3];
double M_T = dynamicVars[4];
double T_0 = dynamicVars[5];
double T_1 = dynamicVars[6];
double T_2 = dynamicVars[7];
double C = dynamicVars[8];
double C_N = dynamicVars[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double T_t = T_0 + T_1 + T_2 + C + C_N;

pd[5] = T_1/(K_2T + T_1);
pd[6] = -(T_1/(K_2T + T_1));
}

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){

res_function_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);

int p_index = (int)constants[40];
double constants_only[40];
int jj;
for (jj = 0; jj < 40; jj++){
constants_only[jj] = constants[jj];}
double *dc_dp = &sens_y[10];
double *dcdot_dp = &sens_yp[10];
double *local_dres_dp = &sens_res[10];
int ii;
for(ii = 0; ii < 10; ii++){
local_dres_dp[ii] = 0;}
switch(p_index)
{
case 0 : dres_dk_2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 1 : dres_dK_2T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 2 : dres_dK_IT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 3 : dres_dK_dT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 4 : dres_dK_3T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 5 : dres_dK_2P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 6 : dres_dK_IP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 7 : dres_dK_dP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 8 : dres_dk_4_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 9 : dres_dK_mP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 10 : dres_dK_mT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 11 : dres_dV_4T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 12 : dres_dv_sP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 13 : dres_dv_mP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 14 : dres_dn_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 15 : dres_dv_sT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 16 : dres_dv_mT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 17 : dres_dV_1T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 18 : dres_dV_4P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 19 : dres_dV_1P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 20 : dres_dk_1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 21 : dres_dv_dT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 22 : dres_dk_3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 23 : dres_dk_sT_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 24 : dres_dv_dP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 25 : dres_dk_sP_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 26 : dres_dk_dN_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 27 : dres_dK_1T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 28 : dres_dK_1P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 29 : dres_dK_3P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 30 : dres_dV_3P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 31 : dres_dk_dC_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 32 : dres_dK_4T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 33 : dres_dV_2P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 34 : dres_dV_3T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 35 : dres_dk_d_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 36 : dres_dK_4P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 37 : dres_dV_2T_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
}
double local_dres_dc[100] = {0};
dres_dc_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dc);
int row, col;
for(row = 0; row < 10; row++){
for(col = 0; col < 10; col++){
sens_res[row+10] += local_dres_dc[row + col*10]*dc_dp[col];}}
double local_dres_dcdot[100] = {0};
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
  if (*reseed) {init_genrand(seed);}

  short stch[30][10] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                    {0, -1, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, -1, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, -1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, -1, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, -1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, -1, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, -1, 0, 0},
                    {0, 0, 0, -1, 0, 0, 0, -1, 1, 0},
                    {0, 0, 0, 1, 0, 0, 0, 1, -1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, -1, 1},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, -1},
                    {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
                    {0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
                    {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, -1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, -1}};
  short depd[30+1][30] = {{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0},
                    {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
                    {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
                    {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                    {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

  double time = *time_ptr;
  double sd = (*rmsd_ptr)*(*rmsd_ptr)*10.;
  double stop_time = *stop_time_ptr;
  double dt=0.0;

  double dv0[10];
  for (i=0;i<10;i++) {dv0[i]=dv[i];}

  int rxnInd = 30;
  double propensity, selection, props[30], av[2];
  while (time < stop_time) {
    av[0]=dv[1] + dv[2] + dv[3] + dv[8] + dv[9];
    av[1]=dv[5] + dv[6] + dv[7] + dv[8] + dv[9];

    if (depd[rxnInd][0]) {props[0]=cv[13]*pow(cv[7], cv[15])/(pow(cv[7], cv[15]) + pow(dv[9], cv[15]));}
    if (depd[rxnInd][1]) {props[1]=cv[16]*pow(cv[3], cv[15])/(pow(cv[3], cv[15]) + pow(dv[9], cv[15]));}
    if (depd[rxnInd][2]) {props[2]=cv[26]*dv[0];}
    if (depd[rxnInd][3]) {props[3]=cv[24]*dv[4];}
    if (depd[rxnInd][4]) {props[4]=cv[20]*dv[1]/(cv[29] + dv[1]);}
    if (depd[rxnInd][5]) {props[5]=cv[34]*dv[2]/(cv[6] + dv[2]);}
    if (depd[rxnInd][6]) {props[6]=cv[31]*dv[2]/(cv[30] + dv[2]);}
    if (depd[rxnInd][7]) {props[7]=cv[19]*dv[3]/(cv[37] + dv[3]);}
    if (depd[rxnInd][8]) {props[8]=cv[18]*dv[5]/(cv[28] + dv[5]);}
    if (depd[rxnInd][9]) {props[9]=cv[38]*dv[6]/(cv[2] + dv[6]);}
    if (depd[rxnInd][10]) {props[10]=cv[35]*dv[6]/(cv[5] + dv[6]);}
    if (depd[rxnInd][11]) {props[11]=cv[12]*dv[7]/(cv[33] + dv[7]);}
    if (depd[rxnInd][12]) {props[12]=cv[23]*dv[3]*dv[7];}
    if (depd[rxnInd][13]) {props[13]=cv[9]*dv[8];}
    if (depd[rxnInd][14]) {props[14]=cv[21]*dv[8];}
    if (depd[rxnInd][15]) {props[15]=cv[1]*dv[9];}
    if (depd[rxnInd][16]) {props[16]=cv[14]*dv[0]/(cv[10] + dv[0]);}
    if (depd[rxnInd][17]) {props[17]=cv[17]*dv[4]/(cv[11] + dv[4]);}
    if (depd[rxnInd][18]) {props[18]=cv[25]*dv[3]/(cv[8] + dv[3]);}
    if (depd[rxnInd][19]) {props[19]=cv[22]*dv[7]/(cv[4] + dv[7]);}
    if (depd[rxnInd][20]) {props[20]=cv[36]*dv[0];}
    if (depd[rxnInd][21]) {props[21]=cv[36]*dv[1];}
    if (depd[rxnInd][22]) {props[22]=cv[36]*dv[2];}
    if (depd[rxnInd][23]) {props[23]=cv[36]*dv[3];}
    if (depd[rxnInd][24]) {props[24]=cv[36]*dv[4];}
    if (depd[rxnInd][25]) {props[25]=cv[36]*dv[5];}
    if (depd[rxnInd][26]) {props[26]=cv[36]*dv[6];}
    if (depd[rxnInd][27]) {props[27]=cv[36]*dv[7];}
    if (depd[rxnInd][28]) {props[28]=cv[32]*dv[8];}
    if (depd[rxnInd][29]) {props[29]=cv[27]*dv[9];}

    propensity = 0.0;
    for (i=0;i<30;i++) {
      propensity += props[i];}
   if (propensity<=0.0) {
      dt = stop_time-time;
      time = stop_time;
      break;
   }

    dt = -log(1.0-genrand_real32())/propensity;
    time += dt;

    selection = propensity * genrand_real32();

    for (rxnInd=0; rxnInd<30; rxnInd++) {
      if (selection < props[rxnInd]) {break;}
      else {selection -= props[rxnInd];}}

    for (i=0;i<10;i++) {dv[i]+=stch[rxnInd][i];}

    double _sd = 0.0;
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

double cell = constants[0];
double k_2 = constants[1];
double K_2T = constants[2];
double K_IT = constants[3];
double K_dT = constants[4];
double K_3T = constants[5];
double K_2P = constants[6];
double K_IP = constants[7];
double K_dP = constants[8];
double k_4 = constants[9];
double K_mP = constants[10];
double K_mT = constants[11];
double V_4T = constants[12];
double v_sP = constants[13];
double v_mP = constants[14];
double n = constants[15];
double v_sT = constants[16];
double v_mT = constants[17];
double V_1T = constants[18];
double V_4P = constants[19];
double V_1P = constants[20];
double k_1 = constants[21];
double v_dT = constants[22];
double k_3 = constants[23];
double k_sT = constants[24];
double v_dP = constants[25];
double k_sP = constants[26];
double k_dN = constants[27];
double K_1T = constants[28];
double K_1P = constants[29];
double K_3P = constants[30];
double V_3P = constants[31];
double k_dC = constants[32];
double K_4T = constants[33];
double V_2P = constants[34];
double V_3T = constants[35];
double k_d = constants[36];
double K_4P = constants[37];
double V_2T = constants[38];
double _tau = constants[39];

double M_P = dynamicVars[0];
double M_P_deriv_wrt_time = yprime[0];
double P_0 = dynamicVars[1];
double P_0_deriv_wrt_time = yprime[1];
double P_1 = dynamicVars[2];
double P_1_deriv_wrt_time = yprime[2];
double P_2 = dynamicVars[3];
double P_2_deriv_wrt_time = yprime[3];
double M_T = dynamicVars[4];
double M_T_deriv_wrt_time = yprime[4];
double T_0 = dynamicVars[5];
double T_0_deriv_wrt_time = yprime[5];
double T_1 = dynamicVars[6];
double T_1_deriv_wrt_time = yprime[6];
double T_2 = dynamicVars[7];
double T_2_deriv_wrt_time = yprime[7];
double C = dynamicVars[8];
double C_deriv_wrt_time = yprime[8];
double C_N = dynamicVars[9];
double C_N_deriv_wrt_time = yprime[9];

double P_t = P_0 + P_1 + P_2 + C + C_N;
double P_t_deriv_wrt_time = C_deriv_wrt_time + P_1_deriv_wrt_time + P_2_deriv_wrt_time + P_0_deriv_wrt_time + C_N_deriv_wrt_time;
double T_t = T_0 + T_1 + T_2 + C + C_N;
double T_t_deriv_wrt_time = T_2_deriv_wrt_time + T_0_deriv_wrt_time + T_1_deriv_wrt_time + C_deriv_wrt_time + C_N_deriv_wrt_time;

}
