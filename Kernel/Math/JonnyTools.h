#ifndef JONNYTOOLS_H
#define JONNYTOOLS_H
#include <stdlib.h>

/* For use with qsort */
int qsort_cmp(const void *x, const void *y)
{
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}

/* Returns the 100*perc-th percentile of array.
 * array should be sorted, e.g. 
 *   qsort(xs, n, sizeof(double), qsort_cmp);
 * and perc should be a number between 0 and 1. 
 */
double getPercentile(const double* array, size_t nel, double perc)
{
    size_t lower_ind = floor(perc * (nel-1));
    size_t upper_ind = ceil(perc * (nel-1));
    double lambda = (perc * (nel-1)) - lower_ind;
    double lower_x = array[lower_ind];
    double upper_x = array[upper_ind];
    double percentile = (1-lambda)*lower_x + lambda*upper_x;
    // fprintf(stderr, "nel %d perc %f lower_ind %d upper_ind %d lambda %f lower_x %f upper_x %f percentile %f\n", 
    //         nel, perc, lower_ind, upper_ind, lambda, lower_x, upper_x, percentile);
    return percentile;
}

#endif
