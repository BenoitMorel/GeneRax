#include <cstdio> 
#include <cstdlib> 
#include <cassert>
#include <cmath>
#include "AUTest.hpp"
#include <cstring>
#include <string>
#include <limits>
#include <iostream>

#include <gsl/mygsl.h>

// this was originally for vectorization
size_t get_safe_upper_limit(size_t s)
{
  return s;
}

#define BIGX            20.0                                 /* max value to represent exp (x) */
#define LOG_SQRT_PI     0.5723649429247000870717135          /* log (sqrt (pi)) */
#define I_SQRT_PI       0.5641895835477562869480795          /* 1 / sqrt (pi) */
#define Z_MAX           6.0                                  /* maximum meaningful z value */
#define ex(x)           (((x) < -BIGX) ? 0.0 : exp (x))

/************** Normalz: probability of normal z value *********************/

/*
 * ALGORITHM: Adapted from a polynomial approximation in:
 *                         Ibbetson D, Algorithm 209
 *                                                 Collected Algorithms of the CACM 1963 p. 616
 *                                                                 Note:
 *                                                                                         This routine has six digit accuracy, so it is only useful for absolute
 *                                                                                                                 z values < 6.  For z values >= to 6.0, Normalz() returns 0.0.
 *                                                                                                                  */

double Normalz(double z) /*VAR returns cumulative probability from -oo to z VAR normal z value */ {
  double y, x, w;

  if (z == 0.0)
    x = 0.0;
  else {
    y = 0.5 * fabs(z);
    if (y >= (Z_MAX * 0.5))
      x = 1.0;
    else if (y < 1.0) {
      w = y*y;
      x = ((((((((0.000124818987 * w
                        - 0.001075204047) * w + 0.005198775019) * w
                    - 0.019198292004) * w + 0.059054035642) * w
                - 0.151968751364) * w + 0.319152932694) * w
            - 0.531923007300) * w + 0.797884560593) * y * 2.0;
    } else {
      y -= 2.0;
      x = (((((((((((((-0.000045255659 * y
                                  + 0.000152529290) * y - 0.000019538132) * y
                              - 0.000676904986) * y + 0.001390604284) * y
                          - 0.000794620820) * y - 0.002034254874) * y
                      + 0.006549791214) * y - 0.010557625006) * y
                  + 0.011630447319) * y - 0.009279453341) * y
              + 0.005353579108) * y - 0.002141268741) * y
          + 0.000535310849) * y + 0.999936657524;
    }
  }
  return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}


double computePValueChiSquare(double x, int df) /* x: obtained chi-square value,  df: degrees of freedom */ {
  double a, y, s;
  double e, c, z;
  int even; /* true if df is an even number */

  if (x <= 0.0 || df < 1)
    return (1.0);

  y = 1;

  a = 0.5 * x;
  even = (2 * (df / 2)) == df;
  if (df > 1)
    y = ex(-a);
  s = (even ? y : (2.0 * Normalz(-sqrt(x))));
  if (df > 2) {
    x = 0.5 * (df - 1.0);
    z = (even ? 1.0 : 0.5);
    if (a > BIGX) {
      e = (even ? 0.0 : LOG_SQRT_PI);
      c = log(a);
      while (z <= x) {
        e = log(z) + e;
        s += ex(c * z - a - e);
        z += 1.0;
      }
      return (s);
    } else {
      e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
      c = 0.0;
      while (z <= x) {
        e = e * (a / z);
        c = c + e;
        z += 1.0;
      }
      return (c * y + s);
    }
  } else
    return (s);
}


int cntdist2(double *vec, int bb, double t)
{
  int i,i0,i1;

  i0=0; i1=bb-1;
  if(t < vec[0]) return 0;
  else if(vec[bb-1] <= t) return bb;

  while(i1-i0>1) {
    i=(i0+i1)/2;
    if(vec[i] <= t) i0=i;
    else i1=i;
  }

  return i1;
}

double log3(double x)
{
  double y,z1,z2,z3,z4,z5;
  if(fabs(x)>1.0e-3) {
    y=-log(1.0-x);
  } else {
    z1=x; z2=z1*x; z3=z2*x; z4=z3*x; z5=z4*x;
    y=((((z5/5.0)+z4/4.0)+z3/3.0)+z2/2.0)+z1;
  }
  return y;
}

int mleloopmax=30;
double mleeps=1e-10;
int mlecoef(double *cnts, double *rr, double bb, int kk,
	    double *coef0, /* set initinal value (size=2) */
	    double *lrt, int *df, /* LRT statistic */
        double *se
	    )
{
  int i,m,loop;
  double coef[2], update[2];
  double d1f, d2f, d11f, d12f, d22f; /* derivatives */
  double v11, v12, v22; /* inverse of -d??f */
  double a,e;
  double s[kk], r[kk],c[kk], b[kk],z[kk],p[kk],d[kk],g[kk],h[kk];

  m=0;
  for(i=0;i<kk;i++)
    {
      r[m]=rr[i]; s[m]=sqrt(rr[i]); c[m]=cnts[i]*bb; b[m]=bb;
      m++;
    }
  if(m<2) return 1;

  coef[0]=coef0[0]; /* signed distance */
  coef[1]=coef0[1]; /* curvature */

  for(loop=0;loop<mleloopmax;loop++) {
    d1f=d2f=d11f=d12f=d22f=0.0;
    for(i=0;i<m;i++) {
      z[i]=coef[0]*s[i]+coef[1]/s[i];
      p[i]=gsl_cdf_ugaussian_P(-z[i]);
      d[i]=gsl_ran_ugaussian_pdf(z[i]);
      if(p[i]>0.0 && p[i]<1.0) {
	g[i]=d[i]*( d[i]*(-c[i]+2.0*c[i]*p[i]-b[i]*p[i]*p[i])/
		    (p[i]*p[i]*(1.0-p[i])*(1.0-p[i]))
		    + z[i]*(c[i]-b[i]*p[i])/(p[i]*(1.0-p[i])) );
	h[i]=d[i]*(c[i]-b[i]*p[i])/(p[i]*(1.0-p[i]));
      } else { g[i]=h[i]=0.0; }
      d1f+= -h[i]*s[i]; d2f+= -h[i]/s[i];
      d11f+= g[i]*r[i]; d12f+= g[i]; d22f+= g[i]/r[i];
    }

    a=d11f*d22f-d12f*d12f;
    if(a==0.0) {
      return 2;
    }
    v11=-d22f/a; v12=d12f/a; v22=-d11f/a;

    /* Newton-Raphson update */
    update[0]=v11*d1f+v12*d2f; update[1]=v12*d1f+v22*d2f;
    coef[0]+=update[0]; coef[1]+=update[1];

    /* check convergence */
    e=-d11f*update[0]*update[0]-2.0*d12f*update[0]*update[1]
      -d22f*update[1]*update[1];

    if(e<mleeps) break;
  }

  /* calc log-likelihood */
  *lrt=0.0; *df=0;
  for(i=0;i<m;i++) {
    if(p[i]>0.0 && p[i]<1.0) {
      *df+=1;
      if(c[i]>0.0) a=c[i]*log(c[i]/b[i]/p[i]); else a=0.0;
      if(c[i]<b[i]) a+=(b[i]-c[i])*(log3(p[i])-log3(c[i]/b[i]));
      *lrt += a;
    }
  }
  *lrt *= 2.0; *df -= 2;

  /* write back the results */
  coef0[0]=coef[0]; coef0[1]=coef[1];
  *se = v11 + v22 - 2*v12;
//  vmat[0][0]=v11;vmat[0][1]=vmat[1][0]=v12;vmat[1][1]=v22; 
  if(loop==mleloopmax || *df< 0) i=1; else i=0;
  return i;
}

double cntdist3(double *vec, int bb, double t)
{
  double p,n;
  int i;
  i=cntdist2(vec,bb,t)-1; /* to find vec[i] <= t < vec[i+1] */
  n=(double)bb;
  if(i<0) {
    if(vec[1]>vec[0]) p=0.5+(t-vec[0])/(vec[1]-vec[0]);
    else p=0.0;
  } else if(i<bb-1) {
    if(vec[i+1]>vec[i]) p=0.5+(double)i+(t-vec[i])/(vec[i+1]-vec[i]);
    else p=0.5+(double)i; /* <- should never happen */
  } else {
    if(vec[bb-1]-vec[bb-2]>0) p=n-0.5+(t-vec[bb-1])/(vec[bb-1]-vec[bb-2]);
    else p=n;
  }
  if(p>n) p=n; else if(p<0.0) p=0.0;
  return p;
}

void doWeightedLeastSquare(int n, double *w, double *a, double *b, double *c, double &x, double &y, double &se) {
    int k;
    double BC = 0.0, AB = 0.0, AC = 0.0, A2 = 0.0, B2 = 0.0;
    double denom;
    for (k = 0; k < n; k++) {
        double wa = w[k]*a[k];
        double wb = w[k]*b[k];
        AB += wa*b[k];
        BC += wb*c[k];
        AC += wa*c[k];
        A2 += wa*a[k];
        B2 += wb*b[k];
    }
    denom = 1.0/(AB*AB - A2*B2);
    x = (BC*AB - AC*B2) * denom;
    y = (AC*AB - BC*A2) * denom;
    
    se = -denom*(B2+A2+2*AB);
    assert(se >= 0.0);
}

template<class T1, class T2>
void quicksort(T1* arr, int left, int right, T2* arr2 = NULL) {
  int i = left, j = right;
  T1 pivot = arr[(left + right) / 2];

  /* partition */
  while (i <= j) {
    while (arr[i] < pivot)
      i++;
    while (arr[j] > pivot)
      j--;
    if (i <= j) {
      T1 tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;
      if (arr2) {
        T2 tmp2 = arr2[i];
        arr2[i] = arr2[j];
        arr2[j] = tmp2;
      }
      i++;
      j--;
    }
  };

  /* recursion */
  if (left < j)
    quicksort(arr, left, j, arr2);
  if (i < right)
    quicksort(arr, i, right, arr2);
}




void performAUTest(std::vector<Likelihoods> likelihoods,
    size_t nboot,
    std::vector<double> &au_pvalues) {

  size_t ntrees = likelihoods.size();
  assert(ntrees != 0);
  size_t nelems = likelihoods[0].size();
  au_pvalues.resize(ntrees);
  size_t ran_seed = 42;    
  srand(ran_seed); 

  /* STEP 1: specify scale factors */
  size_t nscales = 10;
  double r[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
  double rr[] = {sqrt(0.5), sqrt(0.6), sqrt(0.7), sqrt(0.8), sqrt(0.9), 1.0, 
        sqrt(1.1), sqrt(1.2), sqrt(1.3), sqrt(1.4)};
    double rr_inv[] = {sqrt(1/0.5), sqrt(1/0.6), sqrt(1/0.7), sqrt(1/0.8), sqrt(1/0.9), 1.0, 
        sqrt(1/1.1), sqrt(1/1.2), sqrt(1/1.3), sqrt(1/1.4)};
        
    /* STEP 2: compute bootstrap proportion */
    
    size_t nptn = likelihoods[0].size();
    size_t maxnptn = get_safe_upper_limit(nptn);
    
    double *treelhs;
    treelhs = new double[ntrees*nscales*nboot];
    


    int *rstream = nullptr; //randstream;
    size_t boot;
    std::vector<int> boot_sample(maxnptn, 0);
    std::vector<double> boot_sample_dbl(maxnptn);
    for (size_t k = 0; k < nscales; k++) {
		    for (boot = 0; boot < nboot; boot++) {
            if (r[k] == 1.0 && boot == 0) {
                // 2018-10-23: get one of the bootstrap sample as the original alignment
              std::fill(boot_sample.begin(), 
                  boot_sample.end(), 
                  1);
            } else {
              std::fill(boot_sample.begin(), 
                  boot_sample.end(), 
                  0);
              for (size_t i = 0; i < nelems; ++i) {
                boot_sample[rand() % nelems]++;
              }
            }
            for (size_t ptn = 0; ptn < maxnptn; ptn++) {
                boot_sample_dbl[ptn] = boot_sample[ptn];
            }
            double max_lh = -std::numeric_limits<double>::max();
            double second_max_lh = -std::numeric_limits<double>::max();
            int max_tid = -1;
            for (size_t tid = 0; tid < ntrees; tid++) {
                const auto &pattern_lh = likelihoods[tid];
                //double *pattern_lh = pattern_lhs + (tid*maxnptn);
                double tree_lh = 0.0;
                for (size_t ptn = 0; ptn < nptn; ptn++) {
                  tree_lh += pattern_lh[ptn] * boot_sample_dbl[ptn];
                }
                // rescale lh
                tree_lh /= r[k];
                
                // find the max and second max
                if (tree_lh > max_lh) {
                    second_max_lh = max_lh;
                    max_lh = tree_lh;
                    max_tid = tid;
                } else if (tree_lh > second_max_lh)
                    second_max_lh = tree_lh;
                    
                treelhs[(tid*nscales+k)*nboot + boot] = tree_lh; 
            }
            
            // compute difference from max_lh
            for (size_t tid = 0; tid < ntrees; tid++) { 
                if (tid != max_tid) {
                    treelhs[(tid*nscales+k)*nboot + boot] = max_lh - treelhs[(tid*nscales+k)*nboot + boot];
                } else {
                    treelhs[(tid*nscales+k)*nboot + boot] = second_max_lh - max_lh;
                }
            }
        } // for boot
        
        // sort the replicates
        for (size_t tid = 0; tid < ntrees; tid++) {
            quicksort<double,int>(treelhs + (tid*nscales+k)*nboot, 0, nboot-1);
        }
        
    } // for scale



    
    /* STEP 3: weighted least square fit */
    
    double *cc = new double[nscales];
    double *w = new double[nscales];
    double *this_bp = new double[nscales];
    std::cout << "TreeID\tAU\tRSS\td\tc" << std::endl;
    for (size_t tid = 0; tid < ntrees; tid++) {
        double *this_stat = treelhs + tid*nscales*nboot;
        double xn = this_stat[(nscales/2)*nboot + nboot/2], x;
        double c, d; // c, d in original paper
        int idf0 = -2;
        double z = 0.0, z0 = 0.0, thp = 0.0, th = 0.0, ze = 0.0, ze0 = 0.0;
        double pval, se;
        int df;
        double rss = 0.0;
        int step;
        const int max_step = 30;
        bool failed = false;
        for (step = 0; step < max_step; step++) {
            x = xn;
            int num_k = 0;
            for (size_t k = 0; k < nscales; k++) {
                this_bp[k] = cntdist3(this_stat + k*nboot, nboot, x) / nboot;
                if (this_bp[k] <= 0 || this_bp[k] >= 1) {
                    cc[k] = w[k] = 0.0;
                } else {
                    double bp_val = this_bp[k];
                    cc[k] = -gsl_cdf_ugaussian_Pinv(bp_val);
                    double bp_pdf = gsl_ran_ugaussian_pdf(cc[k]);
                    w[k] = bp_pdf*bp_pdf*nboot / (bp_val*(1.0-bp_val));
                    num_k++;
                }
            }
            df = num_k-2;
            if (num_k >= 2) {
                // first obtain d and c by weighted least square
                doWeightedLeastSquare(nscales, w, rr, rr_inv, cc, d, c, se);

                // maximum likelhood fit
                double coef0[2] = {d, c};
                int mlefail = mlecoef(this_bp, r, nboot, nscales, coef0, &rss, &df, &se);
                
                if (!mlefail) {
                    d = coef0[0];
                    c = coef0[1];
                }

                se = gsl_ran_ugaussian_pdf(d-c)*sqrt(se);
                
                // second, perform MLE estimate of d and c
    //            OptimizationAUTest mle(d, c, nscales, this_bp, rr, rr_inv);
    //            mle.optimizeDC();
    //            d = mle.d;
    //            c = mle.c;

                /* STEP 4: compute p-value according to Eq. 11 */
                pval = gsl_cdf_ugaussian_Q(d-c);
                z = -pval;
                ze = se;
                // compute sum of squared difference
                rss = 0.0;
                for (size_t k = 0; k < nscales; k++) {
                    double diff = cc[k] - (rr[k]*d + rr_inv[k]*c);
                    rss += w[k] * diff * diff;
                }
                
            } else {
                // not enough data for WLS
                int num0 = 0;
                for (size_t k = 0; k < nscales; k++)
                    if (this_bp[k] <= 0.0) num0++;
                if (num0 > nscales/2)
                    pval = 0.0;
                else
                    pval = 1.0;
                se = 0.0;
                d = c = 0.0;
                rss = 0.0;
            }

            
            if(df < 0 && idf0 < 0) { failed = true; break;} /* degenerated */
            
            if ((df < 0) || (idf0 >= 0 && (z-z0)*(x-thp) > 0.0 && fabs(z-z0)>0.1*ze0)) {
                th=x;
                xn=0.5*x+0.5*thp;
                continue;
            }
            if(idf0 >= 0 && (fabs(z-z0)<0.01*ze0)) {
                if(fabs(th)<1e-10) 
                    xn=th; 
                else th=x;
            } else 
                xn=0.5*th+0.5*x;
            au_pvalues[tid] = pval;
            thp=x; 
            z0=z;
            ze0=ze;
            idf0 = df;
            if(fabs(x-th)<1e-10) break;
        } // for step
        
        
        if (step == max_step) {
            failed = true;
        }
        
        double pchi2 = (failed) ? 0.0 : computePValueChiSquare(rss, df);
        std::cout << tid+1 << "\t" << au_pvalues[tid] << "\t" << rss << "\t" << d << "\t" << c;
        
        // warning if p-value of chi-square < 0.01 (rss too high)
        if (pchi2 < 0.01) 
            std::cout << " !!!";
        std::cout << std::endl;
    }
    
    delete [] this_bp;
    delete [] w;
    delete [] cc;
//    delete [] bp;
//
#ifdef PLOP
#endif
}

