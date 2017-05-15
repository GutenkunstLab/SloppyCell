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

void dres_dr1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dr3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dtao_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dr4_0_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dv1_0_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

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

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

residual[0] = Nucleus*r4*v4*2.0 - Cytoplasm*r1*v1*D*60.0 - yprime[0];
residual[1] = Cytoplasm*r1*v1*D*60.0 - Cytoplasm*pow(v2, 2.0)*2.0 - yprime[1];
residual[2] = Cytoplasm*pow(v2, 2.0) - Cytoplasm*r3*v3 - yprime[2];
residual[3] = Cytoplasm*r3*v3 - Nucleus*r4*v4 - yprime[3];
residual[4] = -yprime[4];
residual[5] = -yprime[5];
residual[6] = -yprime[6];
}

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;


double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

}

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

pd[0] = -(Cytoplasm*r1*D*60.0);
pd[21] = Nucleus*r4*2.0;
pd[28] = Nucleus*v4*2.0;
pd[35] = -(Cytoplasm*r1*v1*60.0*time);
pd[42] = -(Cytoplasm*r1*v1*60.0);
pd[1] = Cytoplasm*r1*D*60.0;
pd[8] = -(Cytoplasm*v2*4.0);
pd[36] = Cytoplasm*r1*v1*time*60.0;
pd[43] = Cytoplasm*r1*v1*60.0;
pd[9] = Cytoplasm*v2*2.0;
pd[16] = -(Cytoplasm*r3);
pd[17] = Cytoplasm*r3;
pd[24] = -(Nucleus*r4);
pd[31] = -(Nucleus*v4);
}

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

pd[0] = -1;
pd[8] = -1;
pd[16] = -1;
pd[24] = -1;
pd[32] = -1;
pd[40] = -1;
pd[48] = -1;
}

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *pd, double *cj_ptr, double *constants, int *intpar){
double cj = *cj_ptr;

dres_dc_function_(time_ptr, dynamicVars, yprime, constants, pd);

double local_dres_dcdot[7*7] = {0};
dres_dcdot_function_(time_ptr, dynamicVars, yprime, constants, local_dres_dcdot);

int ii;
for(ii=0; ii < 49; ii++){
  pd[ii] += cj*local_dres_dcdot[ii];}
}

void dres_dr1_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

pd[0] = -(Cytoplasm*v1*D*60.0);
pd[1] = Cytoplasm*v1*D*60.0;
}

void dres_dr3_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

pd[2] = -(Cytoplasm*v3);
pd[3] = Cytoplasm*v3;
}

void dres_dtao_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

}

void dres_dr4_0_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

}

void dres_dv1_0_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v2 = dynamicVars[1];
double v3 = dynamicVars[2];
double v4 = dynamicVars[3];
double r4 = dynamicVars[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_offset = dynamicVars[6];

double data1 = v2 + 2.0*v3;
double data2 = v1 + v2 + 2.0*v3;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double D = smooth_D;
double frac_v1 = v1/v1_0;
double frac_v2 = v2/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v4 = 2.0*v4/v1_0;

}

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){

res_function_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);

int p_index = (int)constants[7];
double constants_only[7];
int jj;
for (jj = 0; jj < 7; jj++){
constants_only[jj] = constants[jj];}
double *dc_dp = &sens_y[7];
double *dcdot_dp = &sens_yp[7];
double *local_dres_dp = &sens_res[7];
int ii;
for(ii = 0; ii < 7; ii++){
local_dres_dp[ii] = 0;}
switch(p_index)
{
case 0 : dres_dr1_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 1 : dres_dr3_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 2 : dres_dtao_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 3 : dres_dr4_0_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
case 4 : dres_dv1_0_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
}
double local_dres_dc[49] = {0};
dres_dc_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dc);
int row, col;
for(row = 0; row < 7; row++){
for(col = 0; col < 7; col++){
sens_res[row+7] += local_dres_dc[row + col*7]*dc_dp[col];}}
double local_dres_dcdot[49] = {0};
dres_dcdot_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dcdot);
for(row = 0; row < 7; row++){
for(col = 0; col < 7; col++){
sens_res[row+7] += local_dres_dcdot[row + col*7]*dcdot_dp[col];}}
}

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double dynamicVars[7];
double yprime[7];
int ii;
for(ii = 0; ii < 7; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
res_function_(time_ptr, dynamicVars, yprime, cj_ptr, residual, ires_ptr, constants, ipar);
}

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double dynamicVars[7];
double yprime[7];
int ii;
for(ii = 0; ii < 7; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
root_func_(neq_ptr, time_ptr, dynamicVars, yprime, nrt_ptr, root_devs, constants, ipar);
}

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){
double sens_y[14];
double sens_yp[14];
int ii;
for(ii = 0; ii < 7; ii++){
sens_y[ii] = max(exp(sens_y_log[ii]), DBL_MIN);
sens_yp[ii] = sens_yp_log[ii] * sens_y[ii];}
for(ii = 7; ii < 14; ii++){
sens_y[ii] = sens_y_log[ii];
sens_yp[ii] = sens_yp_log[ii];}
sens_rhs_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);
}

void integrate_stochastic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory) {
return;}

void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double time = *time_ptr;

double Cytoplasm = constants[0];
double Nucleus = constants[1];
double r1 = constants[2];
double r3 = constants[3];
double tao = constants[4];
double r4_0 = constants[5];
double v1_0 = constants[6];

double v1 = dynamicVars[0];
double v1_deriv_wrt_time = yprime[0];
double v2 = dynamicVars[1];
double v2_deriv_wrt_time = yprime[1];
double v3 = dynamicVars[2];
double v3_deriv_wrt_time = yprime[2];
double v4 = dynamicVars[3];
double v4_deriv_wrt_time = yprime[3];
double r4 = dynamicVars[4];
double r4_deriv_wrt_time = yprime[4];
double smooth_D_slope = dynamicVars[5];
double smooth_D_slope_deriv_wrt_time = yprime[5];
double smooth_D_offset = dynamicVars[6];
double smooth_D_offset_deriv_wrt_time = yprime[6];

double data1 = v2 + 2.0*v3;
double data1_deriv_wrt_time = v2_deriv_wrt_time + v3_deriv_wrt_time*2.0;
double data2 = v1 + v2 + 2.0*v3;
double data2_deriv_wrt_time = v1_deriv_wrt_time + v2_deriv_wrt_time + v3_deriv_wrt_time*2.0;
double smooth_D = smooth_D_slope*time + smooth_D_offset;
double smooth_D_deriv_wrt_time = smooth_D_offset_deriv_wrt_time + time*smooth_D_slope_deriv_wrt_time + smooth_D_slope;
double D = smooth_D;
double D_deriv_wrt_time = 0.0;
double frac_v1 = v1/v1_0;
double frac_v1_deriv_wrt_time = v1_deriv_wrt_time/v1_0;
double frac_v2 = v2/v1_0;
double frac_v2_deriv_wrt_time = v2_deriv_wrt_time/v1_0;
double frac_v3 = 2.0*v3/v1_0;
double frac_v3_deriv_wrt_time = v3_deriv_wrt_time*2.0/v1_0;
double frac_v4 = 2.0*v4/v1_0;
double frac_v4_deriv_wrt_time = v4_deriv_wrt_time*2.0/v1_0;

root_devs[0] = (time > 0.0) - 0.5;
root_devs[1] = (time > 2.0) - 0.5;
root_devs[2] = (time > 4.0) - 0.5;
root_devs[3] = (time > 6.0) - 0.5;
root_devs[4] = (time > 8.0) - 0.5;
root_devs[5] = (time > 10.0) - 0.5;
root_devs[6] = (time > 12.0) - 0.5;
root_devs[7] = (time > 14.0) - 0.5;
root_devs[8] = (time > 16.0) - 0.5;
root_devs[9] = (time > 18.0) - 0.5;
root_devs[10] = (time > 20.0) - 0.5;
root_devs[11] = (time > 25.0) - 0.5;
root_devs[12] = (time > 30.0) - 0.5;
root_devs[13] = (time > 40.0) - 0.5;
root_devs[14] = (time > 50.0) - 0.5;
root_devs[15] = (time > tao) - 0.5;
}
