#include <math.h>

double log_gaussian_prior_integrand(double u, double ak, double bk, double mulB, double siglB, double T, double B_best, double lB_best);

double log_gaussian_prior_integrand(double u, double ak, double bk, double mulB, double siglB, double T, double B_best, double lB_best){
    double B = exp(u) * B_best;
    double lB = u + lB_best;
    double ret = exp(-ak/(2*T) * pow(B-B_best,2) - pow(lB-mulB,2)/(2 * pow(siglB,2)));
    return ret;
}
