// Program to list irreducible cubic forms associated to elliptic
// curves of given conductor with no 2-torsion

#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>

#include "cubic_utils.h"

//#define DEBUG
//#define VERBOSE

int main ()
{
  initprimes("PRIMES");
  bigint N, N0, N1, D;
  int n, alpha, beta;
  // while(cerr << "Enter conductor (0 to stop): ", cin >> n, n>0)
  while(cin >> n, n>0)
    {
#ifdef VERBOSE
      cout << "N = " << n << endl;
#endif
      if (!is_valid_conductor(n))
        continue;

      N0 = N = BIGINT(n);
      vector<bigint> badprimes = pdivs(N0);
      vector<bigint> badprimes_not_2 = badprimes;
      vector<bigint> badprimes_not_2_or_3 = badprimes;
      alpha = beta = 0;
      if (div(2,N0))
        {
          alpha = divide_out(N0,2);
          badprimes_not_2 = badprimes_not_2_or_3 = pdivs(N0);
        }
      if (div(3,N0))
        {
          beta = divide_out(N0,3);
          badprimes_not_2_or_3 = pdivs(N0);
        }

      // we have a local test at odd primes p dividing N once, sp we collect these here
      vector<bigint> local_test_primes;
      for (auto pi = badprimes_not_2.begin(); pi!=badprimes_not_2.end(); ++pi)
        if (val(*pi,N)==1)
          local_test_primes.push_back(*pi);

#ifdef VERBOSE
      cout << "alpha, beta, N0 = " << alpha << ", " << beta << ", " << N0 << endl;
#endif
      vector<int> alpha0s = alpha0list(alpha);
      vector<int> beta0s = beta0list(beta);
      vector<int> signs = {1,-1};
      vector<bigint> N1s = posdivs(N0);
#ifdef DEBUG
      cout << "alpha0 in "<<alpha0s<<endl;
      cout << "beta0  in "<<beta0s<<endl;
      cout << "N1     in "<<N1s<<endl;
#endif
      for (auto N1i  = N1s.begin(); N1i!=N1s.end(); ++N1i)
        {
          bigint D23 = *N1i;
#ifdef DEBUG
          cout << "D23 = " << D23 << endl;
#endif
          for (auto alpha0  = alpha0s.begin(); alpha0!=alpha0s.end(); ++alpha0)
            {
              bigint D3 = D23*pow(2,*alpha0);
#ifdef DEBUG
              cout << "alpha0 = " << (*alpha0) << ", D3 = " << D3 << endl;
#endif
              for (auto beta0  = beta0s.begin(); beta0!=beta0s.end(); ++beta0)
                {
                  bigint absD = D3*pow(3,*beta0);
                  int mod3condition = ((*beta0) >=3);
#ifdef DEBUG
                  cout << "beta0 = " << (*beta0) << ", |D| = " << absD << endl;
#endif
                  for (auto s = signs.begin(); s!=signs.end(); s++)
                    {
                      D = (*s) * absD;
#ifdef DEBUG
                      cout << "D = " << D << endl;
#endif
                      vector<cubic> glist = reduced_cubics(D,
                                                           0, // 0 to exclude reducibles
                                                           1, // 1 for GL2-equivalence
                                                           0); // verbosity level
                      for (auto gi = glist.begin(); gi!=glist.end(); ++gi)
                        {
                          int locally_ok = 1;
                          if (mod3condition)
                            locally_ok = (div(3,gi->b()) && (div(3,gi->c())));
                          for (auto pi=local_test_primes.begin(); pi!=local_test_primes.end() && locally_ok; ++pi)
                            locally_ok = local_test(*gi, D, *pi);
                          if (locally_ok)
                            cout << N << " " << D << " " << (*gi) << endl;
#ifdef DEBUG
                          else
                            cout << (*gi) << " fails a local test" << endl;
#endif
                        }
                    }
                }
            }
        }
    }

}
