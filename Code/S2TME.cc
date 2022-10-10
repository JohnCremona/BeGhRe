// Program to list irreducible cubic forms associated to elliptic
// curves of conductor supported on finite sets of primes {2,3} U S
// (with no 2-torsion)

#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>

#include "cubic_utils.h"

//#define DEBUG

int main ()
{
  initprimes("PRIMES");

  int nS;
  cerr << "Number of primes in S (excluding 2 and 3)?" << endl;
  cin >> nS;
  cerr << "Enter " << nS << " primes: " << endl;
  vector<bigint> S(nS);
  bigint M(1);
  for (int i=0; i<nS; i++)
    {
      cin >> S[i];
      M *= S[i];
    }

  vector<int> alpha0s = {0,2,3,4};
  vector<int> beta0s = {0,1,3,4,5};
  vector<int> signs = {1,-1};
  vector<bigint> N1s = posdivs(M*M);
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
#ifdef DEBUG
              cout << "beta0 = " << (*beta0) << ", |D| = " << absD << endl;
#endif
              for (auto s = signs.begin(); s!=signs.end(); s++)
                {
                  bigint D = (*s) * absD;
#ifdef DEBUG
                  cout << "D = " << D << endl;
#endif
                  vector<cubic> glist = reduced_cubics(D,
                                                       0, // 0 to exclude reducibles
                                                       1, // 1 for GL2-equivalence
                                                       0); // verbosity level
                  for (auto gi = glist.begin(); gi!=glist.end(); ++gi)
                    {
                      if (gi->is_primitive())
                        {
                          cout << D << " ";
                          vector<int> data = {*s, *alpha0, *beta0};
                          cout << "[" << data << ", " << D23 << "] ";
                          cout << (*gi) << endl;
                        }
#ifdef DEBUG
                      else
                        cout << (*gi) << " is not primitive" << endl;
#endif
                    }
                }
            }
        }
    }

}
