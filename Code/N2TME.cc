// Program to list irreducible cubic forms associated to elliptic
// curves of given conductor with no 2-torsion

// Usage:
// ./N2TME <N>  # for one conductor N
// ./N2TME <N1> <N2>  # for all conductors from N1 to N2 inclusive

#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>

#include "TME.h"

#define VERBOSE 0

int main (int argc, char *argv[])
{
  if ( (argc < 2) || (argc > 3) )
    {
      cerr << "Usage: N2TME N or N2TME N1 N2" <<endl;
      return 0;
    }
  long n1 = atoi(argv[1]), n2;
  if (argc == 3)
    {
      n2 = atoi(argv[2]);
    }
  else
    {
      n2 = n1;
    }

#if VERBOSE
  if (argc==2)
    cerr << "TM equations for conductor " << n1 << endl;
  else
    cerr << "TM equations for conductors from " << n1 << " to " << n2 << endl;
#endif
  initprimes("PRIMES");

  for (long n=n1; n<=n2; n++)
    {
#if VERBOSE
      cout << "N = " << n << ": " <<flush;
#endif
      if (!is_valid_conductor(n))
        {
#if VERBOSE
          cout << " -- not a valid conductor" << endl;
#endif
          continue;
        }
      vector<TM_eqn> TM_eqns = get_TMeqnsN(n);

#if VERBOSE
      int neqns = TM_eqns.size();
      if (neqns)
        cout << neqns;
      else
        cout <<"No";
      cout << " TM equations for conductor " << n << endl;
#endif

      for (auto T = TM_eqns.begin(); T!=TM_eqns.end(); ++T)
        {
          string s(*T);
          cout << s << endl;
        }
    }
}
