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

void dres_dKeq_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dk2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar);

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar);

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar);

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar);

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

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

alg_derivs_res[0] = (((0.0) + ((-(1.0))*(yp[(, <, _, a, s, t, ., I, n, d, e, x,  , o, b, j, e, c, t,  , a, t,  , 0, x, 7, f, 8, 8, 5, 0, 1, 6, c, 7, 3, 0, >, )]))) + (((Keq) + (1.0))*(yp[(, <, _, a, s, t, ., I, n, d, e, x,  , o, b, j, e, c, t,  , a, t,  , 0, x, 7, f, 8, 8, 5, 0, 1, 6, c, 9, 4, 0, >, )])));
}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;

double X0, X1, T, S1;
double cell, Keq, k1, k2;
double S2;
dynamicVars[3] = alg_vals[0];

cell = constants[0];
Keq = constants[1];
k1 = constants[2];
k2 = constants[3];

X0 = dynamicVars[0];
X1 = dynamicVars[1];
T = dynamicVars[2];
S1 = dynamicVars[3];

S2 = ((Keq)*(S1));

residual[0] = (((S2) + (S1)) - (T));
}

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

pd[0] = (-(k1));
pd[13] = ((k2)*(Keq));
pd[2] = k1;
pd[14] = (-((k2)*(Keq)));
}

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

pd[0] = -1;
pd[5] = -1;
pd[10] = -1;
pd[15] = -1;
}

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *delta, double *pd, double *cj_ptr, double *h_ptr, double *wt, double *constants, int *intpar){
double cj = *cj_ptr;

double local_dres_dcdot[4*4] = {0};
int ii;

dres_dc_function_(time_ptr, dynamicVars, yprime, constants, pd);

dres_dcdot_function_(time_ptr, dynamicVars, yprime, constants, local_dres_dcdot);

for(ii=0; ii < 16; ii++){
  pd[ii] += cj*local_dres_dcdot[ii];}
}

void dres_dKeq_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

pd[1] = ((k2)*(S1));
pd[2] = (-((k2)*(S1)));
}

void dres_dk1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

pd[0] = (-(X0));
pd[2] = X0;
}

void dres_dk2_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];

double X0 = dynamicVars[0];
double X1 = dynamicVars[1];
double T = dynamicVars[2];
double S1 = dynamicVars[3];

double S2 = ((Keq)*(S1));

pd[1] = S2;
pd[2] = (-(S2));
}

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){

int p_index = (int)constants[4];
double constants_only[4];
int jj;
double *dc_dp;
double *dcdot_dp;
double *local_dres_dp;
int ii;
double local_dres_dc[16] = {0};
double local_dres_dcdot[16] = {0};
int row, col;

res_function_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);

for (jj = 0; jj < 4; jj++){
constants_only[jj] = constants[jj];}
dc_dp = &sens_y[4];
dcdot_dp = &sens_yp[4];
local_dres_dp = &sens_res[4];
for(ii = 0; ii < 4; ii++){
local_dres_dp[ii] = 0;}
switch(p_index)
{
case 0 : dres_dKeq_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 1 : dres_dk1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 2 : dres_dk2_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
}
dres_dc_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dc);
for(row = 0; row < 4; row++){
for(col = 0; col < 4; col++){
sens_res[row+4] += local_dres_dc[row + col*4]*dc_dp[col];}}
dres_dcdot_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dcdot);
for(row = 0; row < 4; row++){
for(col = 0; col < 4; col++){
sens_res[row+4] += local_dres_dcdot[row + col*4]*dcdot_dp[col];}}
}

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double dynamicVars[4];
double yprime[4];
int ii;
for(ii = 0; ii < 4; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
res_function_(time_ptr, dynamicVars, yprime, cj_ptr, residual, ires_ptr, constants, ipar);
}

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double dynamicVars[4];
double yprime[4];
int ii;
for(ii = 0; ii < 4; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
root_func_(neq_ptr, time_ptr, dynamicVars, yprime, nrt_ptr, root_devs, constants, ipar);
}

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){
double sens_y[8];
double sens_yp[8];
int ii;
for(ii = 0; ii < 4; ii++){
sens_y[ii] = max(exp(sens_y_log[ii]), DBL_MIN);
sens_yp[ii] = sens_yp_log[ii] * sens_y[ii];}
for(ii = 4; ii < 8; ii++){
sens_y[ii] = sens_y_log[ii];
sens_yp[ii] = sens_yp_log[ii];}
sens_rhs_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);
}

void integrate_stochastic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory) {
return;}

void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double time = *time_ptr;

double cell = constants[0];
double Keq = constants[1];
double k1 = constants[2];
double k2 = constants[3];

double X0 = dynamicVars[0];
double X0_deriv_wrt_time = yprime[0];
double X1 = dynamicVars[1];
double X1_deriv_wrt_time = yprime[1];
double T = dynamicVars[2];
double T_deriv_wrt_time = yprime[2];
double S1 = dynamicVars[3];
double S1_deriv_wrt_time = yprime[3];

double S2 = ((Keq)*(S1));
double S2_deriv_wrt_time = ((Keq)*(S1_deriv_wrt_time));

root_devs[0] = ((time > 18.0) && (X1 > 0.4)) - 0.5;
root_devs[1] = ((time > 10.0) && (X1 > 0.4)) - 0.5;
root_devs[2] = (X1 <= 0.5) - 0.5;
root_devs[3] = (time < 1.0) - 0.5;
root_devs[4] = (X1 > 0.4) - 0.5;

root_devs[5] = (time > 18.0) - 0.5;

root_devs[6] = (X1 > 0.4) - 0.5;

root_devs[7] = (time > 10.0) - 0.5;

}
