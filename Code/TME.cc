#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>
#include <assert.h>

#include "cubic_utils.h"
#include "TME.h"

vector<int> signs = {1,-1};
vector<vector<int>> alpha0list = {{2}, {2,3}, {2,4}, {2,3,4}, {2,3,4}, {2,3}, {2,3,4}, {3,4}, {3}};
vector<vector<int>> beta0list = {{0}, {0,1}, {0,1,3}, {3}, {4}, {5}};
vector<int> powersof2 = {1,2,4,8,16,32,64,128,256};
vector<int> powersof3 = {1,3,9,27,81,243};

void Ndata::init()
{
  support = pdivs(N);
  N0 = N;
  alpha = divide_out(N0, 2);
  beta  = divide_out(N0, 3);
  for (auto pi = support.begin(); pi!=support.end(); ++pi)
    {
      bigint p = *pi;
      if (p>3)
        {
          if (val(p,N)==1)
            Mprimes.push_back(p);
          else
            Aprimes.push_back(p);
        }
    }
}

Ddata::Ddata(const Ndata& Ndat, const bigint& D23, int al, int be, int sg)
  :NN(Ndat), D0(D23), alpha(al), beta(be), s(sg)
{
  D = D0 * powersof2[alpha] * powersof3[beta] * s;
}

vector<Ddata> get_discriminants(const Ndata& NN)
{
  vector<int> alpha0s = alpha0list[NN.alpha];
  vector<int> beta0s = beta0list[NN.beta];
  vector<bigint> N1s = posdivs(NN.N0);

  vector<Ddata> Dlist;
  for (auto N1i  = N1s.begin(); N1i!=N1s.end(); ++N1i)
    for (auto alpha0  = alpha0s.begin(); alpha0!=alpha0s.end(); ++alpha0)
      for (auto beta0  = beta0s.begin(); beta0!=beta0s.end(); ++beta0)
        for (auto s = signs.begin(); s!=signs.end(); s++)
          Dlist.push_back(Ddata(NN, *N1i, *alpha0, *beta0, *s));
  return Dlist;
}

// these two maps are not yet used in the code which follows

map<pair<int,int>, vector<int>> alpha_map = { {{0,2}, {0,3}},
                                              {{2,2}, {1}},
                                              {{2,4}, {0,1}},
                                              {{3,2}, {1,2}},
                                              {{3,3}, {2}},
                                              {{3,4}, {0,1}},
                                              {{4,4}, {0,1}},
                                              {{5,2}, {0}},
                                              {{5,3}, {1}},
                                              {{6,4}, {0,1}},
                                              {{7,3}, {0}},
                                              {{7,3}, {0}},
                                              {{8,3}, {1}}};
map<pair<int,int>, vector<int>> beta_map = { {{0,0}, {0}},
                                             {{2,3}, {0}},
                                             {{3,3}, {0,1}},
                                             {{4,4}, {0,1}},
                                             {{5,5}, {0,1}}};

vector<TM_RHS> get_RHS(const Ddata& D)
{
  bigint one(1), two(2), three(3);
  vector<bigint> alist;
  alist.push_back(one);
  vector<bigint> plist;

  // p=2
  int aN = D.NN.alpha, aD = D.alpha;
  switch(aN*10+aD) {
  case 2:
    alist = multiply_list_by_powers(two, {0,3}, alist);
    break;
  case 22:
  case 53:
  case 83:
    alist = multiply_list_by_powers(two, {1}, alist);
    break;
  case 24:
  case 34:
  case 44:
  case 64:
    alist = multiply_list_by_powers(two, {0,1}, alist);
    break;
  case 32:
    alist = multiply_list_by_powers(two, {1,2}, alist);
    break;
  case 33:
    alist = multiply_list_by_powers(two, {2}, alist);
    break;
  case 52:
  case 73:
  case 74:
    // alist = multiply_list_by_powers(two, {0}, alist);
    break;
  default:
    plist.push_back(two);
  }
  // p=3
  int bN = D.NN.beta, bD = D.beta;
  switch(bN*10+bD) {
  case 0:
  case 23:
    // alist = multiply_list_by_powers(3, {0}, alist);
    break;
  default:
    if (bN>=3)
      alist = multiply_list_by_powers(three, {0,1}, alist);
    else
      plist.push_back(three);
  }

  // p>3, multiplicative
  for (auto pi=D.NN.Mprimes.begin(); pi!=D.NN.Mprimes.end(); ++pi)
    plist.push_back(*pi);

  // p>3, additive
  for (auto pi=D.NN.Aprimes.begin(); pi!=D.NN.Aprimes.end(); ++pi)
    {
      bigint p = *pi;
      if (val(p,D.D)==2)
        alist = multiply_list_by_powers(p, {0,1}, alist);
      else
        plist.push_back(p);
    }

  vector<TM_RHS> RHSs;
  for (auto ai=alist.begin(); ai!=alist.end(); ++ai)
    RHSs.push_back(TM_RHS(*ai, plist));
  return RHSs;
}

TM_RHS::operator string() const
{
  ostringstream s;
  s << a << ",[";
  for (auto pi=plist.begin(); pi!=plist.end(); ++pi)
    {
      if (pi!=plist.begin())
        s << ",";
      s << *pi;
    }
  s << "]";
  return s.str();
}

vector<cubic> get_cubics(const Ddata& DD)
{
  return reduced_cubics(DD.D,
                        0, // 0 to exclude reducibles
                        1, // 1 for GL2-equivalence
                        0); // verbosity level
}

//#define DEBUG_LOCAL_TEST

int TM_eqn::local_test()
{
#ifdef DEBUG_LOCAL_TEST
  cout<<"Local testing of "<<(string)(*this)<<endl;
#endif
  if (!modaCheck(F, RHS.a))
    {
#ifdef DEBUG_LOCAL_TEST
      cout<<" - fails modaCheck"<<endl;
#endif
        return 0;
    }

  vector<bigint> newprimes;

  for (auto pi=RHS.plist.begin(); pi!=RHS.plist.end(); ++pi)
    {
      bigint p = *pi;
      if (has_roots_mod(F, p)) // then F mod p roots, so we keep p as an RHS prime
        {
#ifdef DEBUG_LOCAL_TEST
          cout<<" - p = "<<p<<" is kept"<<endl;
#endif
          newprimes.push_back(p);
        }
      else // F mod p has no roots...
        {
          // check (13) is satisfied, F(u,v) is not 0 mod p
          if ((p>=3) && (val(p,DD.NN.N)==1) && !div(p,DD.D))
            {
#ifdef DEBUG_LOCAL_TEST
              cout<<" - p = "<<p<<" fails condition (13)"<<endl;
#endif
              return 0;
            }
          // Now the equation has no solutions with a positive power
          // of p on RHS, and we check whether it is satisfiable mod p
          // with p not dividing the RHS:
          if (!modpCheck(F, RHS.a, RHS.plist, p))
            {
#ifdef DEBUG_LOCAL_TEST
              cout<<" - p = "<<p<<" fails"<<endl;
#endif
              return 0; // no it is not, discard this RHS
            }
          // otherwise it is, we do *not* keep this prime and hence
          // only look for TME solutions with no p dividing the RHS
#ifdef DEBUG_LOCAL_TEST
          cout<<" - p = "<<p<<" passes and is removed"<<endl;
#endif
        }
    }
  RHS.plist = newprimes;
  return 1;
}

TM_eqn::operator string() const
{
  ostringstream s;
  F.output(s);
  return s.str() + "," + (string)RHS;
}

vector<TM_eqn> get_TMeqnsD(const Ddata& DD)
{
  vector<TM_eqn> TMeqns;
  vector<cubic> Flist = get_cubics(DD);
  vector<TM_RHS> RHSlist = get_RHS(DD);
  for (auto Fi=Flist.begin(); Fi!=Flist.end(); ++Fi)
    {
      for (auto RHSi=RHSlist.begin(); RHSi!=RHSlist.end(); ++RHSi)
        {
          TM_eqn T(DD, *Fi, *RHSi);
          if (T.local_test())
            TMeqns.push_back(T);
        }
    }
  return TMeqns;
}

vector<TM_eqn> get_TMeqnsN(const Ndata& NN)
{
  vector<TM_eqn> TMeqns;
  vector<Ddata> Dlist = get_discriminants(NN);
  for (auto Di=Dlist.begin(); Di!=Dlist.end(); ++Di)
    {
      vector<TM_eqn> TMeqnsD = get_TMeqnsD(*Di);
      TMeqns.insert(TMeqns.end(), TMeqnsD.begin(), TMeqnsD.end());
    }
  return TMeqns;
}

// for p||N and p not dividing D=disc(F) we require that F(u,v)=0 (mod
// p) has a nontrivial solution:
int local_test(const cubic& F, const Ddata& DD, const bigint& p)
{
  return (div(p,DD.D) || has_roots_mod(F, p));
}

