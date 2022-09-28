// Some simple utilities which should really be in eclib's cubic functions:
bigint content(const cubic& F);
int is_primitive(const cubic& F);
// divide a cubic by a constant factor (which should divide all the coefficients)
cubic divide_out(const cubic& F, const bigint& g);

// sl2/gl2-reduce a cubic
void sl2_reduce(cubic& F);
void gl2_reduce(cubic& F);
// for an sl2-reduced cubic, normalise w.r.t. diag(1,-1), i.e. [a,b,c,d] --> [a,-b,c,-d]
void gl2_normalise(cubic& F);

// Tests for sl2/gl2-equivalence and equality:
int identical(const cubic& F, const cubic& G);
int sl2_equivalent(const cubic& F, const cubic& G);
int gl2_equivalent(const cubic& F, const cubic& G);

// affine roots of F mod q, assuming leading coefficient F.a() is
// nonzero:
vector<bigint> roots_mod(const cubic& F, const bigint& q);

// Return 1 iff F has a projective root mod q:
int has_roots_mod(const cubic& F, const bigint& q);

// multiply all integers in a list by a constant:
vector<bigint> multiply_list(const bigint& a, const vector<bigint>& L);

// multiply all integers in L1 by all in L2:
vector<bigint> multiply_lists(const vector<bigint>& L1, const vector<bigint>& L2);

// multiply all integers in L by p^e for e in exponents:
vector<bigint> multiply_list_by_powers(const bigint& p, const vector<int>& exponents, const vector<bigint>& L);

// for q prime > 3, returns a list of representatives of the values of
// F(u,v) mod q modulo cubes
vector<bigint> image_mod_cubes(const cubic& F, const bigint& q);

// Similar to AG's magma function of the same name. Return 1 iff there
// exist (u,v) not (0,0) mod q and exponents e such that
// F(u,v)=a*prod(p^e) mod q.
int modpCheck(const cubic& F, const bigint& a, const vector<bigint>& primes, const bigint& q);

// similar to AG's modpcheckDivRHS. Return 1 iff there exists
// primitive (u,v) such that F(u,v)=0 (mod a).
int modaCheck(const cubic& F, const bigint& a);
