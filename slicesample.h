#ifndef _SLICESAMPLE
#define _SLICESAMPLE

double slicesample(double init, double (*logpdf)(double), double width,
                   int burnin);

#endif
