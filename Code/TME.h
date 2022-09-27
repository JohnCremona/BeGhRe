// class for holding a conductor and associated factorization data:

class Ndata {
public:
  bigint N;
  bigint N0; // prime-to-6 part of N
  int alpha; // valuation of 2
  int beta;  // valuation of 3
  vector<bigint> support; // all prime factors
  vector<bigint> Mprimes; // prime factors > 3 with valuation 2
  vector<bigint> Aprimes; // prime factors > 3 with valuation 1

  // constructors
  Ndata(int conductor)
    :N(conductor)
  {
    init();
  }
  Ndata(const bigint& conductor)
    :N(conductor)
  {
    init();
  }
  void init();
};

// class for holding a discriminant and associated factorization data

class Ddata {
public:
  Ndata NN; // pointer to conductor data
  bigint D;
  bigint D0; // prime-to-6 part of |D|
  int alpha; // valuation of 2
  int beta;  // valuation of 3
  int s;     // sign
  // constructor
  Ddata(const Ndata& Ndat, const bigint& D23, int al, int be, int sg);
};

// class for holding RHS data for a TM-equation

class TM_RHS {
public:
  bigint a;
  vector<bigint> plist;
  // constructor
  TM_RHS(const bigint& aa, const vector<bigint>& pl)
    :a(aa), plist(pl)
  {;}
  // for output:
  operator string() const;
};

class TM_eqn {
public:
  Ddata DD; // discriminant data (which also contains conductor data)
  cubic F;
  TM_RHS RHS;

  TM_eqn(const Ddata& dd, const cubic& f, const TM_RHS& rhs)
    :DD(dd), F(f), RHS(rhs)
  {;}

  // local test: return 0 if impossible, else 1 (and the RHS may have changed)
  int local_test();

  // for output:
  operator string() const;
};

// Return a list of discriminants for one conductor
vector<Ddata> get_discriminants(const Ndata& NN);

// Return a list of RHSs (a, primes) for one discriminant
vector<TM_RHS> get_RHS(const Ddata& D);

// Return a list of irreducible cubic forms up to GL(2,Z)-equivalence) for one discriminant
vector<cubic> get_cubics(const Ddata& DD);

// for p||N and p not dividing D=disc(F) we require that F(u,v)=0 (mod
// p) has a nontrivial solution:
int local_test(const cubic& F, const Ddata& DD, const bigint& p);

// Return all TM equations for one discriminant
vector<TM_eqn> get_TMeqnsD(const Ddata& DD);

// Return all TM equations for one conductor
vector<TM_eqn> get_TMeqnsN(const Ndata& NN);
