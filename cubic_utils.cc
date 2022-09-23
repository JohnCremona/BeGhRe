#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>
#include <assert.h>

#include "cubic_utils.h"

// Some simple utilities which should really be in eclib's cubic functions:

bigint content(const cubic& F)
{
  return gcd(gcd(gcd(F.a(), F.b()), F.c()), F.d());
}

int is_primitive(const cubic& F)
{
  return content(F)==1;
}

// affine roots of F mod q, assuming leading coefficient F.a() is
// nonzero:
vector<bigint> roots_mod(const cubic& F, const bigint& q)
{
  vector<bigint> coeffs = {F.a(), F.b(), F.c(), F.d()};
  return rootsmod(coeffs, q);
}

int has_roots_mod(const cubic& F, const bigint& q)
{
  return div(q,F.a()) || roots_mod(F,q).size() > 0;
}

// Now some more serious functions used in the algorithm

vector<int> alpha0list(int alpha)
{
  switch (alpha)
    {
    case 0:
      return {2};
    case 1: case 5:
      return {2,3};
    case 2:
      return {2,4};
    case 3: case 4: case 6:
      return {2,3,4};
    case 7:
      return {3,4};
    case 8:
      return {3};
    default: // will not happen
      return {};
    }
}

vector<int> beta0list(int beta)
{
  switch (beta)
    {
    case 0:
      return {0};
    case 1:
      return {0,1};
    case 2:
      return {0,1,3};
    default:
      return {beta};
    }
}

// for p||N and p not dividing D=disc(F) we require that F(u,v)=0 (mod
// p) has a nontrivial solution:
int local_test(const cubic& F, const bigint& D, const bigint& p)
{
  if (div(p,D))
    return 1; // no condition to pass
  return has_roots_mod(F, p);
}

int is_cube(const bigint& a, const bigint& q)
{
  if (div(q,a) || div(q,a-1) || div(3,q+1))
    return 1;
  bigint b;
  power_mod(b, a, (q-1)/3, q); // b = a^(q-1)/3 mod q
  return div(q,b);
}

// for q prime, returns a list of representatives of the values of
// F(u,v) mod q modulo cubes
vector<bigint> image_mod_cubes(const cubic& F, const bigint& q)
{
  vector<bigint> images;

  // first see if 0 is a value:

  vector<bigint> coeffs = {F.a(), F.b(), F.c(), F.d()};
  if (has_roots_mod(F, q))
    images.push_back(BIGINT(0));

  // if q=2 (mod 3) or q=3, then all nozero values occur and all are cubes:

  if (div(3,q+1) || q==3)
    {
      images.push_back(BIGINT(1));
      return images;
    }

  // Now q=1 (mod 3) and we must see which of the three nonzero cosets mod cubes are hit:

  // We keep track of the inverses of the images found to simplify the coset check
  vector<bigint> inverses;

  // check F(1,0):
  bigint v = F.a();
  if (v!=0)
    {
      images.push_back(v);
      inverses.push_back(invmod(v, q));
    }

  // check F(0,x) for all x mod q:
  for (bigint x(0); x<q; x++)
    {
      v = F.eval(x);
      if (v==0)
        continue;
      int repeat=0;
      for (auto w = inverses.begin(); !repeat && w!=inverses.end(); ++w)
        repeat = is_cube((*w) * v, q);
      if (!repeat) // we have a new coset
        {
          images.push_back(v);
          inverses.push_back(invmod(v, q));
        }
      if (inverses.size()==3)
        break;
    }
  return images;
}

// (function name as in AG's magma code) Return 1 iff there exist
// (u,v) not (0,0) mod q and exponents e such that F(u,v)=a*prod(p^e)
// mod q.
int modpCheck(const cubic& F, const bigint& a, const vector<bigint>& primes, const bigint& q)
{
  if (div(2,q+1)) // F takes one and hence all nonzero values since all are cubes
    return 1;

  for (auto pi = primes.begin(); pi!=primes.end(); ++pi)
    if (!is_cube(*pi,q)) // powers of p cover all cosets mod cubes
      return 1;

  // Now all p are cubes mod q so can be ignored
  vector<bigint> images = image_mod_cubes(F, q);
  bigint b = invmod(a,q);
  for (auto ci = images.begin(); ci!=images.end(); ++ci)
    {
      if (is_zero(*ci)) // if 0 is a value, ignore it0 as a value (if it is)
        continue;
      if (is_cube(b*(*ci),q))
        return 1;
    }
  return 0;
}

// similar to AG's modpcheckDivRHS. Return 1 iff there exists
// primitive (u,v) such that F(u,v)=0 (mod a).
int modaCheck(const cubic& F, const bigint& a)
{
  if (a==1)
    return 1;
  vector<bigint> plist = pdivs(a);

  if (plist.size()>1) // use CRT
    {
      for (auto pi = plist.begin(); pi!=plist.end(); ++pi)
        {
          bigint q = pow(*pi, val(*pi, a));
          if (!modaCheck(F, q))
            return 0;
        }
      return 1;
    }

  // Now a is a prime power p^e
  bigint p = plist[0];
  if (!has_roots_mod(F,p))
    return 0;
  if (val(p,a)==1) // a=p, nothing more to do
    return 1;

  // Now we must see if any projective root mod p lifts to p^e. Note
  // that in our application we expect p to divide the discriminant of
  // F, so roots mod p will not be simple, otherwise Hensel would make
  // this redundant.

  bigint b = a/p; // = p^{e-1}
  bigint one(1);
  // Test for roots above (1:0) mod p:
  if (div(p,F.a()))
    {
      // Test (u,v) = (1, p*w) for w mod p^{e-1}
      for (bigint w(0); w<b; w++)
        if (div(a,F.eval(one, p*w)))
          return 1;
    }
  // find affine roots r mod p:
  vector<bigint> roots = roots_mod(F,p);
  for (auto r = roots.begin(); r!=roots.end(); ++r)
    {
      // Test (u,v) = (r+p*w, 1) for w mod p^{e-1}
      for (bigint w(0); w<b; w++)
        if (div(a,F.eval(p*w + *r, one)))
          return 1;
    }
  // None of the roots did lift, so we fail
  return 0;
}
