
bigint content(const cubic& F);
int is_primitive(const cubic& F);

// affine roots of F mod q, assuming leading coefficient F.a() is
// nonzero:
vector<bigint> roots_mod(const cubic& F, const bigint& q);

// Return 1 iff F has a projective root mod q:
int has_roots_mod(const cubic& F, const bigint& q);

vector<int> alpha0list(int alpha);
vector<int> beta0list(int beta);

// eclib has a function is_valid_conductor(N) for a long int n

// for p||N and p not dividing D=disc(F) we require that F(u,v)=0 (mod
// p) has a nontrivial solution:
int local_test(const cubic& F, const bigint& D, const bigint& p);

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
