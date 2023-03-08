#define TMB_LIB_INIT R_init_shmm

#undef  CPPAD_FRAMEWORK
#define TMBAD_FRAMEWORK

// Spatial hidden Markov model
#include <TMB.hpp>

template<class Type>
struct parameter_list {
  vector<Type> sd;
  vector<Type> bias;
  Type time;
  // From C++14 use 'auto' for return type
  tmbutils::Concat<vector<Type>, vector<Type>, Type>
  unlist() { return tmbutils::concat(sd, bias, time); }
};

template<class Type>
struct data_likelihood {
  SEXP map;  // List of maps by time
  SEXP map_keep;
  SEXP obs;  // List of observations by time
  SEXP obs_keep;
  SEXP con;
  int ntime;
  parameter_list<Type> pl;
  data_likelihood(SEXP data, parameter_list<Type> pl) : pl(pl) {
    map = getListElement(data, "map"); // list
    obs = getListElement(data, "obs"); // list
    map_keep = getListElement(data, "map_keep"); // list
    obs_keep = getListElement(data, "obs_keep"); // list
    if (LENGTH(map) != LENGTH(obs)) Rf_error("length(maps/obs) mismatch");
    con = getListElement(data, "con"); // matrix
    ntime = LENGTH(map);
  }
  // Get list elt
  vector<Type> getData(SEXP x, Type i) {
    Type t1 = atomic::dynamic_data::sexp_to_double(x);
    Type t2 = atomic::dynamic_data::list_lookup_by_index(t1, i);
    return atomic::dynamic_data::sexp_to_vector(t2);
  }
  // Get matrix
  vector<Type> getData(SEXP x) {
    Type t1 = atomic::dynamic_data::sexp_to_double(x);
    return atomic::dynamic_data::sexp_to_vector(t1);
  }
  vector<Type> operator()(vector<Type> parms) {
    typedef vectorize::vector<Type> vec;
    // Load parms to parameter list
    pl.unlist() = parms;
    vector<Type> sd = pl.sd;
    vector<Type> bias = pl.bias;
    Type time = pl.time;
    // Load data
    matrix<Type> m = getData(map, time);
    vector<Type> o = getData(obs, time);
    o = o - bias;
    int ntyp = o.size();
    m.resize(m.size() / ntyp, ntyp);
    matrix<Type> m_keep = getData(map_keep, time);
    vector<Type> o_keep = getData(obs_keep, time);
    m_keep.resize(m.size() / ntyp, ntyp);
    vec ans;
    for (int i=0; i<ntyp; i++) {
      vec m_col(m.col(i));
      vec m_keep_col(m_keep.col(i));
      vec tmp = dnorm(m_col, o[i], sd[i], true) * o_keep[i] * m_keep_col;
      if (i==0)
        ans = tmp;
      else
        ans += tmp;
    }
    ans = exp(ans);
    vec c = getData(con);
    ans = ans * c;
    return ans;
  }
};

template<class Type>
struct separable_matrix_exponential {
  matrix<Type> EXPGx;
  matrix<Type> EXPGy;
  separable_matrix_exponential(matrix<Type> Gx,
                               matrix<Type> Gy) {
    EXPGx = atomic::expm(Gx);
    EXPGy = atomic::expm(Gy); // Transpose here if not symmetric!
  }
  vector<Type> operator()(vector<Type> u) {
    matrix<Type> um(u);
    int nx = EXPGx.rows();
    int ny = EXPGy.rows();
    um.resize(nx, ny);
    um = TMBad::matmul(EXPGx, um);
    um = TMBad::matmul(um, EXPGy);
    return um.vec();
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(logDx);        // log diffusion in east-west (x) direction
  PARAMETER(logDy);        // log diffusion in north-south (y) direction
  //  DATA_INTEGER(doviterbi); // Type of solver to use
  DATA_MATRIX(datlik);     // Data likelihood
  //  DATA_INTEGER(solvetype); // Type of solver to use
  DATA_INTEGER(ns);        // Number of time steps of solver
  DATA_IVECTOR(iobs);      // Indices to which observations correspond
  DATA_SPARSE_MATRIX(I);   // Identity matrix
  DATA_SCALAR(dt);         // Time step
  //  DATA_INTEGER(m);         // Number of iterations of uniformization
  DATA_SPARSE_MATRIX(Sns); // North-south generator skeleton
  DATA_SPARSE_MATRIX(Sew); // East-west generator skeleton
  DATA_MATRIX(Sx); // Alternative representation
  DATA_MATRIX(Sy); // Alternative representation

  Type Dx = exp(logDx);
  Type Dy = exp(logDy);

  DATA_INTEGER(rawdat); // Rawdata available?

  if (rawdat) {
    /*
      ## Raw Data layout:
      ## - list (by time) of 3D arrays 'map' each indexed by (x, y, datatype)
      ## - list (by time) of 1D arrays 'obs' each indexed by (datatype)
      ## - 2D array 'con' indexed by (x, y)
    */
    PARAMETER_VECTOR(logsd); // obs
    PARAMETER_VECTOR(bias); // obs
    vector<Type> sd = exp(logsd);
    parameter_list<Type> pl = {sd, bias, Type(0)};
    data_likelihood<Type> dlik(this->data, pl);
    atomic::AtomicLocal<data_likelihood<Type> > F(dlik);
    // Keep the old datlik begin and end for now (FIXME)
    for (int i=1; i<dlik.ntime-1; i++) {
      pl.time = Type(i);
      vector<Type> parms = pl.unlist();
      datlik.col(i).array() = F(parms);
    }
    // Store the old datlik begin and end
  }

  REPORT(datlik); // For debugging

  // Matrix exponential (maps vector to vector)
  DATA_STRING(solvetype);
  std::function< vector<Type>( vector<Type> ) > ExpGdt;
  // Sparse case:
  if (solvetype == "sparse") {
    Eigen::SparseMatrix<Type> G = Dx*Sew + Dy*Sns;  // Make generator
    Eigen::SparseMatrix<Type> Gdt = G * dt;
    ExpGdt = sparse_matrix_exponential::expm_generator<Type>(Gdt);
  } else if (solvetype == "separable") {
    // Dense separable case
    matrix<Type> Gx = (Dx*dt)*Sx;
    matrix<Type> Gy = (Dy*dt)*Sy;
    ExpGdt = separable_matrix_exponential<Type>(Gx, Gy);
  } else {
    Rf_error( ("solvetype='" + solvetype + "'" + " not valid").c_str() );
  }

  // ==================================== HMM FILTER
  typedef vectorize::vector<Type> vec; // For faster vectorized math
  matrix<Type> pred;
  matrix<Type> phi;
  vector<Type> psi;
  matrix<Type> smoo; // output from smoothing
  vector<int> track;
  int nobs = datlik.cols();
  int n = datlik.rows();
  // Initialise HMM grids
  pred.resize(n, ns);
  phi.resize(n, ns);
  psi.resize(nobs - 1);
  // First state is at time of first observation
  pred.col(0) = datlik.col(0) / datlik.col(0).sum();
  phi.col(0) = pred.col(0);
  // Filter loop
  for(int t=1; t<ns; t++) {
    // Time update using uniformization algorithm
    vector<Type> svec = phi.col(t-1);
    //vector<Type> predtmp = svec;
    vector<Type> predtmp = ExpGdt(svec);
    pred.col(t) = predtmp; // Store prediction
    if (iobs(t) > 0) {
      // Data update
      int ind = (iobs(t)-1);
      vector<Type> d = datlik.col(ind);
      vec post = vec(d) * vec(predtmp);
      psi(ind - 1) = post.sum();
      vec post_normalized = (post / (psi(ind - 1) + 1e-20)); // Add small value to avoid division by zero
      phi.col(t).array() = post_normalized.array();
    } else {
      // No data update
      phi.col(t) = pred.col(t);
    }
  }
  // Negative log likelihood
  Type ans = -sum(log(psi));

  /* branch is only entered by obj$report(). */
  if(isDouble<Type>::value) {
    smoo.resize(n, ns);
    // Smoothing loop
    smoo.col(ns-1) = phi.col(ns-1);
    vector<Type> rho = datlik.col(ns-1);
    for(int t=1; t < ns; t++) {
      int tt = ns - t - 1;
      vector<Type> tmp = ExpGdt(rho); // FIXME: Transpose G if not symmetric!
      vector<Type> phat = tmp * vector<Type>(phi.col(tt).array()); phat = phat / sum(phat);
      smoo.col(tt) = phat;
      rho = tmp * vector<Type>(datlik.col(tt).array()); rho = rho / sum(rho);
    }
    // Reports
    matrix<Type> predout = pred;
    matrix<Type> phiout = phi;
    matrix<Type> psiout = psi;
    matrix<Type> smooout = smoo;
    REPORT(predout);
    REPORT(phiout);
    REPORT(psi);
    REPORT(smooout);
  }
  
  return ans;

}
