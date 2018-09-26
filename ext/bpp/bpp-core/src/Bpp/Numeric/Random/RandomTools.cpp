//
// File RandomTools.cpp
// Author : Julien Dutheil
// Last modification : Friday Septembre 24 2004
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for numerical calculus.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "RandomTools.h"
#include "Uniform01K.h"
#include "../VectorTools.h"
#include "../NumConstants.h"

#include <iostream>

using namespace bpp;
using namespace std;

RandomFactory* RandomTools::DEFAULT_GENERATOR = new Uniform01K(time(NULL));

// Initiate random seed :
// RandomTools::RandInt RandomTools::r = time(NULL) ;

void RandomTools::setSeed(long seed)
{
  DEFAULT_GENERATOR->setSeed(seed);
}

// Method to get a double random value (between 0 and specified range)
// Note : the number you get is between 0 and entry not including entry !
double RandomTools::giveRandomNumberBetweenZeroAndEntry(double entry, const RandomFactory& generator)
{
  // double tm = r.drawFloatNumber();
  double tm = generator.drawNumber();
  return tm * entry;
}

// Method to get a boolean random value
bool RandomTools::flipCoin(const RandomFactory& generator)
{
  return (RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator) - 0.5) > 0;
}

double RandomTools::randGaussian(double mean, double variance, const RandomFactory& generator)
{
  return RandomTools::qNorm(generator.drawNumber(), mean, sqrt(variance));
}

double RandomTools::randGamma(double dblAlpha, const RandomFactory& generator)
{
  assert(dblAlpha > 0.0);
  if (dblAlpha < 1.0)
    return RandomTools::DblGammaLessThanOne(dblAlpha, generator);
  else if (dblAlpha > 1.0)
    return RandomTools::DblGammaGreaterThanOne(dblAlpha, generator);
  return -log(RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator));
}

double RandomTools::randGamma(double alpha, double beta, const RandomFactory& generator)
{
  double x = RandomTools::randGamma(alpha, generator) / beta;
  return x;
}

double RandomTools::randExponential(double mean, const RandomFactory& generator)
{
  return -mean* log(RandomTools::giveRandomNumberBetweenZeroAndEntry(1, generator));
}

std::vector<size_t> RandomTools::randMultinomial(size_t n, const std::vector<double>& probs)
{
  double s = VectorTools::sum(probs);
  double r;
  double cumprob;
  vector<size_t> sample(n);
  for (unsigned int i = 0; i < n; i++)
  {
    r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
    cumprob = 0;
    bool test = true;
    for (unsigned int j = 0; test &(j < probs.size()); j++)
    {
      cumprob += probs[j] / s;
      if (r <= cumprob)
      {
        sample[i] = j;
        test = false;
      }
    }
    // This test should never be true if probs sum to one:
    if (test)
      sample[i] = probs.size();
  }
  return sample;
}

// ------------------------------------------------------------------------------


double RandomTools::DblGammaGreaterThanOne(double dblAlpha, const RandomFactory& generator)
{
  // Code adopted from David Heckerman
  // -----------------------------------------------------------
  //  DblGammaGreaterThanOne(dblAlpha)
  //
  //  routine to generate a gamma random variable with unit scale and
  //      alpha > 1
  //  reference: Ripley, Stochastic Simulation, p.90
  //  Chang and Feast, Appl.Stat. (28) p.290
  // -----------------------------------------------------------
  double rgdbl[6];

  rgdbl[1] = dblAlpha - 1.0;
  rgdbl[2] = (dblAlpha - (1.0 / (6.0 * dblAlpha))) / rgdbl[1];
  rgdbl[3] = 2.0 / rgdbl[1];
  rgdbl[4] = rgdbl[3] + 2.0;
  rgdbl[5] = 1.0 / sqrt(dblAlpha);

  for ( ; ; )
  {
    double dblRand1;
    double dblRand2;
    do
    {
      dblRand1 = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator);
      dblRand2 = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator);
      if (dblAlpha > 2.5)
        dblRand1 = dblRand2 + rgdbl[5] * (1.0 - 1.86 * dblRand1);
    }
    while (!(0.0 < dblRand1 && dblRand1 < 1.0));

    double dblTemp = rgdbl[2] * dblRand2 / dblRand1;

    if (rgdbl[3] * dblRand1 + dblTemp + 1.0 / dblTemp <= rgdbl[4] ||
        rgdbl[3] * log(dblRand1) + dblTemp - log(dblTemp) < 1.0)
    {
      return dblTemp * rgdbl[1];
    }
  }
  assert(false);
  return 0.0;
}

double RandomTools::DblGammaLessThanOne(double dblAlpha, const RandomFactory& generator)
{
  // routine to generate a gamma random variable with
  // unit scale and alpha < 1
  // reference: Ripley, Stochastic Simulation, p.88
  double dblTemp;
  const double dblexp = exp(1.0);
  for ( ; ; )
  {
    double dblRand0 = giveRandomNumberBetweenZeroAndEntry(1.0, generator);
    double dblRand1 = giveRandomNumberBetweenZeroAndEntry(1.0, generator);
    if (dblRand0 <= (dblexp / (dblAlpha + dblexp)))
    {
      dblTemp = pow(((dblAlpha + dblexp) * dblRand0) /
                    dblexp, 1.0 / dblAlpha);
      if (dblRand1 <= exp(-1.0 * dblTemp))
        return dblTemp;
    }
    else
    {
      dblTemp = -1.0 * log((dblAlpha + dblexp) * (1.0 - dblRand0) / (dblAlpha * dblexp));
      if (dblRand1 <= pow(dblTemp, dblAlpha - 1.0))
        return dblTemp;
    }
  }
  assert(false);
  return 0.0;
}

/******************************************************************************/

// From Yang's PAML package:

/******************************************************************************/

double RandomTools::qNorm(double prob)
{
  double a0 = -.322232431088, a1 = -1, a2 = -.342242088547, a3 = -.0204231210245;
  double a4 = -.453642210148e-4, b0 = .0993484626060, b1 = .588581570495;
  double b2 = .531103462366, b3 = .103537752850, b4 = .0038560700634;
  double y, z = 0, p = prob, p1;

  p1 = (p < 0.5 ? p : 1 - p);
  if (p1 < 1e-20)
    return -9999;

  y = sqrt (log(1 / (p1 * p1)));
  z = y + ((((y * a4 + a3) * y + a2) * y + a1) * y + a0) / ((((y * b4 + b3) * y + b2) * y + b1) * y + b0);
  return p < 0.5 ? -z : z;
}

double RandomTools::qNorm(double prob, double mu, double sigma)
{
  return RandomTools::qNorm(prob) * sigma + mu;
}

double RandomTools::incompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
  size_t i;
  double p = alpha, g = ln_gamma_alpha;
  double accurate = 1e-8, overflow = 1e30;
  double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0;
  vector<double> pn(6);

  if (x == 0)
    return 0;
  if (x < 0 || p <= 0)
    return -1;

  factor = exp(p * log(x) - x - g);
  if (x > 1 && x >= p)
    goto l30;
  /* (1) series expansion */
  gin = 1;  term = 1;  rn = p;
l20:
  rn++;
  term *= x / rn;   gin += term;

  if (term > accurate)
    goto l20;
  gin *= factor / p;
  goto l50;
l30:
  /* (2) continued fraction */
  a = 1 - p;   b = a + x + 1;  term = 0;
  pn[0] = 1;  pn[1] = x;  pn[2] = x + 1;  pn[3] = x * b;
  gin = pn[2] / pn[3];
l32:
  a++;  b += 2;  term++;   an = a * term;
  for (i = 0; i < 2; i++)
  {
    pn[i + 4] = b * pn[i + 2] - an * pn[i];
  }
  if (pn[5] == 0)
    goto l35;
  rn = pn[4] / pn[5];   dif = fabs(gin - rn);
  if (dif > accurate)
    goto l34;
  if (dif <= accurate * rn)
    goto l42;
l34:
  gin = rn;
l35:
  for (i = 0; i < 4; i++)
  {
    pn[i] = pn[i + 2];
  }
  if (fabs(pn[4]) < overflow)
    goto l32;
  for (i = 0; i < 4; i++)
  {
    pn[i] /= overflow;
  }
  goto l32;
l42:
  gin = 1 - factor * gin;

l50:

  return gin;
}


double RandomTools::qChisq(double prob, double v)
{
  double e = .5e-6, aa = .6931471805, p = prob, g;
  double xx, c, ch, a = 0, q = 0, p1 = 0, p2 = 0, t = 0, x = 0, b = 0, s1, s2, s3, s4, s5, s6;

  if (p < .000002 || p > .999998 || v <= 0)
    return -1;

  g = lnGamma (v / 2);
  xx = v / 2;   c = xx - 1;
  if (v >= -1.24 * log(p))
    goto l1;

  ch = pow((p * xx * exp(g + xx * aa)), 1 / xx);
  if (ch - e < 0)
    return ch;
  goto l4;
l1:
  if (v > .32)
    goto l3;
  ch = 0.4;   a = log(1 - p);
l2:
  q = ch;  p1 = 1 + ch * (4.67 + ch);  p2 = ch * (6.73 + ch * (6.66 + ch));
  t = -0.5 + (4.67 + 2 * ch) / p1 - (6.73 + ch * (13.32 + 3 * ch)) / p2;
  ch -= (1 - exp(a + g + .5 * ch + c * aa) * p2 / p1) / t;
  if (fabs(q / ch - 1) - .01 <= 0)
    goto l4;
  else
    goto l2;

l3:
  x = qNorm (p);
  p1 = 0.222222 / v;   ch = v * pow((x * sqrt(p1) + 1 - p1), 3.0);
  if (ch > 2.2 * v + 6)
    ch = -2 * (log(1 - p) - c * log(.5 * ch) + g);
l4:
  q = ch;   p1 = .5 * ch;
  if ((t = incompleteGamma (p1, xx, g)) < 0)
  {
    std::cerr << "err IncompleteGamma" << std::endl;
    return -1;
  }
  p2 = p - t;
  t = p2 * exp(xx * aa + g + p1 - c * log(ch));
  b = t / ch;  a = 0.5 * t - b * c;

  s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420;
  s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520;
  s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520;
  s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040;
  s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520;
  s6 = (120 + c * (346 + 127 * c)) / 5040;
  ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
  if (fabs(q / ch - 1) > e)
    goto l4;

  return ch;
}


double RandomTools::pNorm(double x, double mu, double sigma)
{
  return RandomTools::pNorm((x - mu) / sigma);
}


double RandomTools::pNorm(double x)
{
  const static double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const static double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };
  const static double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const static double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };
  const static double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const static double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  double xden, xnum, temp, del, eps, xsq, y, cum;
  int i;

  eps = 1e-20;

  y = fabs(x);
  if (y <= 0.67448975)   /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
  {
    if (y > eps)
    {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (i = 0; i < 3; ++i)
      {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
      }
    }
    else
      xnum = xden = 0.0;

    temp = x * (xnum + a[3]) / (xden + b[3]);
    cum = 0.5 + temp;
  }
  else if (y <= sqrt(32))
  {
    /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i)
    {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

    xsq = trunc(y * 16) / 16;
    del = (y - xsq) * (y + xsq);
    cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;

    if (x > 0.)
      cum = 1 - cum;
  }
  else if (-37.5193 < x  &&  x < 8.2924)
  {
    xsq = 1.0 / (x * x);
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i)
    {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (1 / sqrt(2 * M_PI) - temp) / y;

    xsq = trunc(x * 16) / 16;
    del = (x - xsq) * (x + xsq);

    cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;

    if (x > 0.)
      cum = 1. - cum;
  }
  else   /* no log_p , large x such that probs are 0 or 1 */
  {
    if (x > 0)
      cum = 1.;
    else
      cum = 0.;
  }

  return cum;
}

double RandomTools::lnBeta(double alpha, double beta)
{
  return lnGamma(alpha) + lnGamma(beta) - lnGamma(alpha + beta);
}

double RandomTools::randBeta(double alpha, double beta, const RandomFactory& generator)
{
  return RandomTools::qBeta(generator.drawNumber(), alpha, beta);
}


double RandomTools::qBeta(double prob, double alpha, double beta)
{
  double lower = NumConstants::VERY_TINY();
  double upper = 1 - NumConstants::VERY_TINY();
  double const1 = 2.30753;
  double const2 = 0.27061;
  double const3 = 0.99229;
  double const4 = 0.04481;


  int swap_tail, i_pb, i_inn;
  double a, adj, logbeta, g, h, pp, prev, qq, r, s, t, tx, w, y, yprev;
  double acu;
  volatile double xinbta;

  if (alpha <= 0. || beta < 0.)
    throw ("RandomTools::qBeta wih non positive parameters");

  if (prob < 0. || prob > 1.)
    throw ("RandomTools::qBeta wih bad probability");

  /* initialize */
  logbeta = lnBeta(alpha, beta);

  /* change tail if necessary;  afterwards   0 < a <= 1/2   */
  if (prob <= 0.5)
  {
    a = prob;  pp = alpha; qq = beta; swap_tail = 0;
  }
  else   /* change tail, swap  alpha <-> beta :*/
  {
    a = 1 - prob;
    pp = beta; qq = alpha; swap_tail = 1;
  }

  /* calculate the initial approximation */

  /* y := {fast approximation of} qnorm(1 - a) :*/
  r = sqrt(-2 * log(a));
  y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);
  if (pp > 1 && qq > 1)
  {
    r = (y * y - 3.) / 6.;
    s = 1. / (pp + pp - 1.);
    t = 1. / (qq + qq - 1.);
    h = 2. / (s + t);
    w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));
    xinbta = pp / (pp + qq * exp(w + w));
  }
  else
  {
    r = qq + qq;
    t = 1. / (9. * qq);
    t = r * pow(1. - t + y * sqrt(t), 3.0);
    if (t <= 0.)
      xinbta = 1. - exp((log1p(-a) + log(qq) + logbeta) / qq);
    else
    {
      t = (4. * pp + r - 2.) / t;
      if (t <= 1.)
        xinbta = exp((log(a * pp) + logbeta) / pp);
      else
        xinbta = 1. - 2. / (t + 1.);
    }
  }

  /* solve for x by a modified newton-raphson method, */
  /* using the function pbeta_raw */

  r = 1 - pp;
  t = 1 - qq;
  yprev = 0.;
  adj = 1;
  /* Sometimes the approximation is negative! */
  if (xinbta < lower)
    xinbta = 0.5;
  else if (xinbta > upper)
    xinbta = 0.5;

  /* Desired accuracy should depend on  (a,p)
   * This is from Remark .. on AS 109, adapted.
   * However, it's not clear if this is "optimal" for IEEE double prec.

   * acu = fmax2(lower, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));

   * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
   * ---- i.e.,  "new acu" = sqrt(old acu)

   */
  double po = pow(10., -13 - 2.5 / (pp * pp) - 0.5 / (a * a));
  acu = (lower > po) ? lower : po;

  tx = prev = 0.;  /* keep -Wall happy */

  for (i_pb = 0; i_pb < 1000; i_pb++)
  {
    y = incompleteBeta(xinbta, pp, qq);
// #ifdef IEEE_754
//     if(!R_FINITE(y))
// #else
//       if (errno)
// #endif
//         ML_ERR_return_NAN;

    y = (y - a) *
        exp(logbeta + r * log(xinbta) + t * log1p(-xinbta));
    if (y * yprev <= 0.)
      prev = (fabs(adj) > lower) ? fabs(adj) : lower;
    g = 1;
    for (i_inn = 0; i_inn < 1000; i_inn++)
    {
      adj = g * y;
      if (fabs(adj) < prev)
      {
        tx = xinbta - adj; /* trial new x */
        if (tx >= 0. && tx <= 1)
        {
          if ((prev <= acu) || (fabs(y) <= acu))
            return swap_tail ? 1 - xinbta : xinbta;
          if (tx != 0. && tx != 1)
            break;
        }
      }
      g /= 3;
    }
    if (fabs(tx - xinbta) < 1e-15 * xinbta)
      return swap_tail ? 1 - xinbta : xinbta;

    xinbta = tx;
    yprev = y;
  }
  // throw Exception("Bad precision in RandomTools::qBeta");

  return swap_tail ? 1 - xinbta : xinbta;
}

double RandomTools::incompleteBeta(double x, double alpha, double beta)
{
  double t;
  double xc;
  double w;
  double y;
  int flag;
  double big;
  double biginv;
  double maxgam;
  double minlog;
  double maxlog;

  big = 4.503599627370496e15;
  biginv = 2.22044604925031308085e-16;
  maxgam = 171.624376956302725;
  minlog = log(NumConstants::VERY_TINY());
  maxlog = log(NumConstants::VERY_BIG());

  if ((alpha <= 0) || (beta <= 0))
    throw Exception("RandomTools::incompleteBeta not valid with non-positive parameters");

  if ((x < 0) || (x > 1))
    throw Exception("RandomTools::incompleteBeta out of bounds limit");

  if (x == 0)
    return 0;

  if (x == 1)
    return 1;

  flag = 0;
  if ((beta * x <= 1.0) && (x <= 0.95))
  {
    return incompletebetaps(alpha, beta, x, maxgam);
  }
  w = 1.0 - x;

  if (x > alpha / (alpha + beta))
  {
    flag = 1;
    t = alpha;
    alpha = beta;
    beta = t;
    xc = x;
    x = w;
  }
  else
  {
    xc = w;
  }
  if (flag == 1 && (beta * x <= 1.0) && (x <= 0.95) )
  {
    t = incompletebetaps(alpha, beta, x, maxgam);
    if (t <= NumConstants::VERY_TINY())
      return 1.0 - NumConstants::VERY_TINY();
    else
      return 1.0 - t;
  }

  y = x * (alpha + beta - 2.0) - (alpha - 1.0);
  if (y < 0.0)
  {
    w = incompletebetafe(alpha, beta, x, big, biginv);
  }
  else
  {
    w = incompletebetafe2(alpha, beta, x, big, biginv) / xc;
  }
  y = alpha * log(x);
  t = beta * log(xc);
  if ( (alpha + beta < maxgam) && (fabs(y) < maxlog) && (fabs(t) < maxlog) )
  {
    t = pow(xc, beta);
    t = t * pow(x, alpha);
    t = t / alpha;
    t = t * w;
    t = t * exp(lnGamma(alpha + beta) - (lnGamma(alpha) + lnGamma(beta)));
    if (flag == 1)
    {
      if (t < NumConstants::VERY_TINY())
        return 1.0 - NumConstants::VERY_TINY();
      else
        return 1.0 - t;
    }
    else
      return t;
  }
  y = y + t + lnGamma(alpha + beta) - lnGamma(alpha) - lnGamma(beta);
  y = y + log(w / alpha);
  if (y < minlog)
  {
    t = 0.0;
  }
  else
  {
    t = exp(y);
  }
  if (flag == 1)
  {
    if (t < NumConstants::VERY_TINY())
      t = 1.0 - NumConstants::VERY_TINY();
    else
      t = 1.0 - t;
  }
  return t;
}


/**********************************************/

double RandomTools::incompletebetafe(double a,
                                     double b,
                                     double x,
                                     double big,
                                     double biginv)
{
  double result;
  double xk;
  double pk;
  double pkm1;
  double pkm2;
  double qk;
  double qkm1;
  double qkm2;
  double k1;
  double k2;
  double k3;
  double k4;
  double k5;
  double k6;
  double k7;
  double k8;
  double r;
  double t;
  double ans;
  double thresh;
  int n;

  k1 = a;
  k2 = a + b;
  k3 = a;
  k4 = a + 1.0;
  k5 = 1.0;
  k6 = b - 1.0;
  k7 = k4;
  k8 = a + 2.0;
  pkm2 = 0.0;
  qkm2 = 1.0;
  pkm1 = 1.0;
  qkm1 = 1.0;
  ans = 1.0;
  r = 1.0;
  n = 0;
  thresh = 3.0 * NumConstants::VERY_TINY();
  do
  {
    xk = -x * k1 * k2 / (k3 * k4);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    xk = x * k5 * k6 / (k7 * k8);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if (qk != 0)
    {
      r = pk / qk;
    }
    if (r != 0)
    {
      t = fabs((ans - r) / r);
      ans = r;
    }
    else
    {
      t = 1.0;
    }
    if (t < thresh)
    {
      break;
    }
    k1 = k1 + 1.0;
    k2 = k2 + 1.0;
    k3 = k3 + 2.0;
    k4 = k4 + 2.0;
    k5 = k5 + 1.0;
    k6 = k6 - 1.0;
    k7 = k7 + 2.0;
    k8 = k8 + 2.0;
    if (fabs(qk) + fabs(pk) > big)
    {
      pkm2 = pkm2 * biginv;
      pkm1 = pkm1 * biginv;
      qkm2 = qkm2 * biginv;
      qkm1 = qkm1 * biginv;
    }
    if ((fabs(qk) < biginv) || (fabs(pk) < biginv))
    {
      pkm2 = pkm2 * big;
      pkm1 = pkm1 * big;
      qkm2 = qkm2 * big;
      qkm1 = qkm1 * big;
    }
    n = n + 1;
  }
  while (n != 300);
  result = ans;
  return result;
}


/*************************************************************************
   Continued fraction expansion #2
   for incomplete beta integral

   Cephes Math Library, Release 2.8:  June, 2000
   Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double RandomTools::incompletebetafe2(double a,
                                      double b,
                                      double x,
                                      double big,
                                      double biginv)
{
  double result;
  double xk;
  double pk;
  double pkm1;
  double pkm2;
  double qk;
  double qkm1;
  double qkm2;
  double k1;
  double k2;
  double k3;
  double k4;
  double k5;
  double k6;
  double k7;
  double k8;
  double r;
  double t;
  double ans;
  double z;
  double thresh;
  int n;

  k1 = a;
  k2 = b - 1.0;
  k3 = a;
  k4 = a + 1.0;
  k5 = 1.0;
  k6 = a + b;
  k7 = a + 1.0;
  k8 = a + 2.0;
  pkm2 = 0.0;
  qkm2 = 1.0;
  pkm1 = 1.0;
  qkm1 = 1.0;
  z = x / (1.0 - x);
  ans = 1.0;
  r = 1.0;
  n = 0;
  thresh = 3.0 * NumConstants::VERY_TINY();
  do
  {
    xk = -z * k1 * k2 / (k3 * k4);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    xk = z * k5 * k6 / (k7 * k8);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if (qk != 0)
    {
      r = pk / qk;
    }
    if (r != 0)
    {
      t = fabs((ans - r) / r);
      ans = r;
    }
    else
    {
      t = 1.0;
    }
    if (t < thresh)
    {
      break;
    }
    k1 = k1 + 1.0;
    k2 = k2 - 1.0;
    k3 = k3 + 2.0;
    k4 = k4 + 2.0;
    k5 = k5 + 1.0;
    k6 = k6 + 1.0;
    k7 = k7 + 2.0;
    k8 = k8 + 2.0;
    if (fabs(qk) + fabs(pk) > big)
    {
      pkm2 = pkm2 * biginv;
      pkm1 = pkm1 * biginv;
      qkm2 = qkm2 * biginv;
      qkm1 = qkm1 * biginv;
    }
    if ((fabs(qk) < biginv) || (fabs(pk) < biginv))
    {
      pkm2 = pkm2 * big;
      pkm1 = pkm1 * big;
      qkm2 = qkm2 * big;
      qkm1 = qkm1 * big;
    }
    n = n + 1;
  }
  while (n != 300);
  result = ans;
  return result;
}


/*************************************************************************
   Power series for incomplete beta integral.
   Use when b*x is small and x not too close to 1.

   Cephes Math Library, Release 2.8:  June, 2000
   Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double RandomTools::incompletebetaps(double a, double b, double x, double maxgam)
{
  double result;
  double s;
  double t;
  double u;
  double v;
  double n;
  double t1;
  double z;
  double ai;

  ai = 1.0 / a;
  u = (1.0 - b) * x;
  v = u / (a + 1.0);
  t1 = v;
  t = u;
  n = 2.0;
  s = 0.0;
  z = NumConstants::VERY_TINY() * ai;
  while (fabs(v) > z)
  {
    u = (n - b) * x / n;
    t = t * u;
    v = t / (a + n);
    s = s + v;
    n = n + 1.0;
  }
  s = s + t1;
  s = s + ai;
  u = a * log(x);
  if ((a + b < maxgam) && (fabs(u) < log(NumConstants::VERY_BIG())))
  {
    t = exp(lnGamma(a + b) - (lnGamma(a) + lnGamma(b)));
    s = s * t * pow(x, a);
  }
  else
  {
    t = lnGamma(a + b) - lnGamma(a) - lnGamma(b) + u + log(s);
    if (t < log(NumConstants::VERY_TINY()))
    {
      s = 0.0;
    }
    else
    {
      s = exp(t);
    }
  }
  result = s;
  return result;
}

/**************************************************************************/

