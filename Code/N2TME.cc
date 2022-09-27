// Program to list irreducible cubic forms associated to elliptic
// curves of given conductor with no 2-torsion

#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>

#include "TME.h"

//#define DEBUG
//#define VERBOSE

int main ()
{
  initprimes("PRIMES");
  bigint N, N0, N1, D;
  int n;
#ifdef VERBOSE
  while(cerr << "Enter conductor (0 to stop): ", cin >> n, n>0)
#else
  while(cin >> n, n>0)
#endif
    {
#ifdef VERBOSE
      cout << "N = " << n << endl;
#endif
      if (!is_valid_conductor(n))
        continue;

      Ndata NN(n);
#ifdef VERBOSE
      cout << "alpha, beta, N0 = " << NN.alpha << ", " << NN.beta << ", " << NN.N0 << endl;
#endif

      vector<TM_eqn> TM_eqns = get_TMeqnsN(NN);
#ifdef VERBOSE
      int neqns = TM_eqns.size();
      if (neqns)
        cout << neqns;
      else
        cout <<"No";
      cout << " TM equations for conductor " << n << endl;
#endif

      for (auto T = TM_eqns.begin(); T!=TM_eqns.end(); ++T)
        {
          cout << n << "," << string(*T) << endl;
        }
    }
}
