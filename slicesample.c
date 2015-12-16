/*******************************************************************************
  An implementation of the slice sampling in standard C.
  Only the case of sampling 1-D data is considered here. However, it should
  be convenient to extend the algorithm for other types of data.

  Each time slicesample() is called, a random number is returned after
  specified iterations of burn-in by starting from a passed initial point.

  Youwei Lu
  tom.l.08.yui [at] gmail [dot] com
******************************************************************************/

#include <math.h>
#include <stdlib.h>
#include "slicesample.h"

double slicesample(double init, double (*logpdf)(double), double width0,
                   int burnin) {
  double llb, rrb, z, lx, width;
  int ii;

  for(ii=0; ii<burnin; ii++) {
    lx = (*logpdf)(init)+log(drand48());
    llb = init-drand48()*width0;
    rrb = init+drand48()*width0;

    width = width0;
    while((*logpdf)(llb) > lx) {
      llb -= drand48()*width;
      width *= 1.5;
    }

    width = width0;
    while((*logpdf)(rrb) > lx) {
      rrb += drand48()*width;
      width *= 1.5;
    }

    do {
      z = llb + drand48()*(rrb-llb);
      if((*logpdf)(z) >= lx) {
        break;
      }

      if(z < init) {
        llb = z;
      }
      else {
        rrb = z;
      }
    }
    while(1);

    init = z;
  }

  return init;
}

#ifdef SLICESAMPLE_TEST

/*******************************************************************************
   Draw a random sample from an unormalized Gauss with mean and standard
   deviation specified by global variables mm and sigma.
   After showing each drawn point, output the sample mean and standard deviation.
 ******************************************************************************/
#include <stdio.h>

double mm, sigma;

double log_unnormalized_normal(double xx) {
  return ( -0.5*(xx-mm)*(xx-mm)/(sigma*sigma) );
}

int main() {
  mm = 10;
  sigma = 100;

  const int sample_size = 100000;
  double sample[sample_size], mean=0, var=0;
  int ii;

  for(ii=0; ii<sample_size; ii++) {
    sample[ii] = slicesample(0, &log_unnormalized_normal, 5, 100);
    printf("Point %d: %.4f\n", ii, sample[ii]);
    mean += sample[ii];
  }

  mean /= sample_size;

  for(ii=0; ii<sample_size; ii++) {
    var += (sample[ii]-mean)*(sample[ii]-mean);
  }
  var /= sample_size-1;

  printf("Mean=%.4f. Std=%.4f.\n", mean, sqrt(var));

  return 0;
}

#endif
