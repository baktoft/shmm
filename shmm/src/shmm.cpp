#define TMB_LIB_INIT R_init_shmm

#undef  CPPAD_FRAMEWORK
#define TMBAD_FRAMEWORK

// Spatial hidden Markov model
// 07.11.2016 - Added tinyad
#include <TMB.hpp>

// /* Class to build generator and project state one step forward  */
// namespace shmm {
//   // Input that is constant (do not depend on parameters)
//   int m;
//   Eigen::SparseMatrix<double> I;
//   Eigen::SparseMatrix<double> Sns;
//   Eigen::SparseMatrix<double> Sew;
//   double dt;
//   vector<double> lgam;
//   matrix<double> datlik;     // Data likelihood
//   int solvetype;             // Type of solver to use
//   int ns;                    // Number of time steps of solver
//   vector<int> iobs;          // Indices to which observations correspond

//   /* 'filter' class calculates the entire likelihood (not just one
//      step). We call the type 'Float' to emphasize that this is a
//      tiny_ad type */
//   template<class Float>
//   struct filter_t {
//     /* Objects we only calculate once during initialization */
//     //    Eigen::SparseMatrix<Float> G;
//     Eigen::SparseMatrix<Float> FPdt; // For uniformization
//     Float F;
//     Eigen::SimplicialLLT<Eigen::SparseMatrix<Float> > solver; // For implicit

//     // HMM grids (used by both 'eval' and 'smoothing')
//     matrix<Float> pred;
//     matrix<Float> phi;
//     vector<Float> psi;
//     matrix<Float> smoo; // output from smoothing
//     vector<int> track;

//     /* Initialize the generator */
//     void initialize(Float Dx, Float Dy) {
//       // Cast required 'double' objects to 'Float'
//       Eigen::SparseMatrix<Float>   I = shmm::  I.cast<Float>();
//       Eigen::SparseMatrix<Float> Sns = shmm::Sns.cast<Float>();
//       Eigen::SparseMatrix<Float> Sew = shmm::Sew.cast<Float>();
//       Float                       dt = shmm::dt;
//       vector<Float>             lgam = shmm::lgam.cast<Float>();
//       // Build generator
//       Eigen::SparseMatrix<Float> G = Dx*Sew + Dy*Sns;  // Make generator
//       if (solvetype == 1) { // Uniformisation
// 	F = 2 * (Dx + Dy); // Absolute largest jump rate, max(abs(diag(G)))
// 	Eigen::SparseMatrix<Float> P = G/F + I; // Sub-stochastic matrix
// 	FPdt = F * P * Float(dt);
//       }
//       else if (solvetype == 2) { // Implicit solving
// 	Eigen::SparseMatrix<Float> A = I - G*dt;
// 	solver.analyzePattern(A);
// 	solver.factorize(A);
//       }
//       else error("'solvetype' can by 1 or 2");
//     }

//     // 1. Uniformization: project state one step forward
//     matrix<Float> forwardProject(matrix<Float> svec){
//       matrix<Float> predtmp = svec; // Initialise
//       for(int i=0; i<m; i++){
// 	svec = svec * FPdt; // Vector x Matrix (this can be optimised?)
// 	predtmp = predtmp + svec/exp(lgam(i)); // exp(lgamma()) is factorial
//       }
//       predtmp = predtmp * exp(-F*dt);
//       predtmp = predtmp / predtmp.sum(); // Ensure total probability mass is 1, should be a minor correction
//       return predtmp;
//     }

//     // 2. Implicit: project state one step forward
//     matrix<Float> forwardProjecti(matrix<Float> svec){
//       matrix<Float> predtmp = solver.solve(svec.transpose());
//       predtmp = predtmp / predtmp.sum(); // Ensure total probability mass is 1, should be a minor correction
//       return predtmp.transpose();
//     }

//     /* Run filter loop */
//     Float eval(Float Dx, Float Dy, bool do_smoothing = false, bool do_viterbi = false) {
//       matrix<Float> datlik = shmm::datlik.cast<Float>();
//       int nobs = datlik.rows();
//       int n = datlik.cols();
//       // Initialize solver or FPdt
//       initialize(Dx, Dy); // FIXME: Put in constructor (?)
//       // Initialise HMM grids
//       pred.resize(ns, n);
//       phi.resize(ns, n);
//       psi.resize(nobs - 1);
//       // First state is at time of first observation
//       pred.row(0) = datlik.row(0) / datlik.row(0).sum();
//       phi.row(0) = pred.row(0);
//       // Filter loop
//       for(int t=1; t<ns; t++) {
//       	// Time update using uniformization algorithm
//       	matrix<Float> svec = phi.row(t-1);
//       	matrix<Float> predtmp = svec;
//         if (solvetype == 1) {
//       	  predtmp = forwardProject(svec);
//         }
//         else if (solvetype == 2) {
//       	  predtmp = forwardProjecti(svec);
//         }
//         else error("'solvetype' can by 1 or 2");
//       	pred.row(t) = predtmp; // Store prediction
//         if (iobs(t) > 0) {
//       	  // Data update
//       	  int ind = (iobs(t)-1);
//       	  matrix<Float> post = pred.row(t).cwiseProduct(datlik.row(ind)); // Element-wise product
//       	  psi(ind - 1) = post.sum(); 
//       	  phi.row(t) = post / (psi(ind - 1) + 1e-20); // Add small value to avoid division by zero
//       	} else {
// 	  // No data update
//       	  phi.row(t) = pred.row(t);
//       	}
//       }
//       // Negative log likelihood
//       Float ans = -sum(log(psi));

//       if(do_smoothing) {
// 	smoo.resize(ns, n);
// 	// Smoothing loop
// 	smoo.row(ns-1) = phi.row(ns-1);
// 	for(int t=1; t < ns; t++) {
// 	  int tt = ns - t;
// 	  // Time update using uniformization algorithm
// 	  matrix<Float> predrow = pred.row(tt);
// 	  for (int i=0; i < n; i++) {
// 	    predrow(0, i) += 1e-10;
// 	  }
// 	  matrix<Float> ratio = smoo.row(tt).cwiseQuotient(predrow);
// 	  //matrix<double> asd(1, n);
// 	  //for (int i=0; i < n; i++){
// 	  //asd(0, i) = std::isnan(asDouble(ratio(0, i)));
// 	  //cout 
// 	  //}
// 	  //matrix<double> asd = std::isnan(asDouble(ratio));
// 	  matrix<Float> ratiotmp = ratio;
// 	  if (solvetype == 1){
// 	    ratiotmp = forwardProject(ratio);
// 	  } else if (solvetype == 2) {
// 	    ratiotmp = forwardProjecti(ratio);
// 	  } else error("'solvetype' can by 1 or 2");
// 	  matrix<Float> post = phi.row(tt-1).cwiseProduct(ratiotmp);
// 	  post = post / (post.sum() + 1e-20);
// 	  smoo.row(tt-1) = post;
// 	}	
//       }


//       if(do_viterbi) {
// 	// === Viterbi ===
// 	// --- Forward sweep ---
// 	track.resize(ns);

// 	matrix<Float> xi(ns, n);
// 	xi.row(0) = phi.row(0);
// 	matrix<Float> zerovec(1, n);
// 	for(int i=0; i < n; i++){
// 	  zerovec(0, i) = 0.0;
// 	}

// 	for(int t=1; t < ns; t++){

// 	  vector<Float> xitmp(n);
// 	  for(int i=0; i < n; i++){
// 	    // Get transition probabilities
// 	    matrix<Float> vec = zerovec;
// 	    vec(0, i) = 1.0;
// 	    matrix<Float> transprob = vec;

// 	    if (solvetype == 1){
// 	      transprob = forwardProject(vec);
// 	    } else if (solvetype == 2) {
// 	      transprob = forwardProjecti(vec);
// 	    } else error("'solvetype' can by 1 or 2");
// 	    vector<Float> tmp = xi.row(t-1).cwiseProduct(transprob);
// 	    // Find max value
// 	    Float mx = tmp(0);
// 	    for(int i=1; i < n; i++){
// 	      if (tmp(i) > mx){
// 		mx = xitmp(i);
// 	      }
// 	    }
// 	    //Float mx = max(tmp.cast<double>).cast<Float>;
// 	    //Float mx = thismax;
// 	    if (iobs(t) > 0){
// 	      // Data update
// 	      int ind = (iobs(t)-1);
// 	      xitmp(i) = mx * datlik(ind, i);
// 	    } else {
// 	      // No data update
// 	      xitmp(i) = mx;
// 	    }
// 	  }
// 	  xitmp = xitmp / xitmp.sum();
// 	  xi.row(t) = xitmp;
// 	}

// 	// --- Backard sweep to find track ---
// 	// Start with final step
// 	vector<Float> xitmp = xi.row(ns-1);
// 	int j = 0;
// 	Float thismax = xitmp(0);
// 	for(int i=1; i < n; i++){
// 	  if (xitmp(i) > thismax){
// 	    thismax = xitmp(i);
// 	    j = i;
// 	  }
// 	}
// 	track(ns-1) = j + 1; 
// 	for(int t=1; t < ns; t++){
// 	  int tt = ns - t;
// 	  // Get transition probabilities
// 	  matrix<Float> vec = zerovec;
// 	  //int ind = CppAD::Integer(track(tt)-1);
// 	  int ind = track(tt)-1;
// 	  vec(0, ind) = 1.0;
// 	  matrix<Float> transprob = vec;
// 	  if (solvetype == 1){
// 	    transprob = forwardProject(vec);
// 	  } else if (solvetype == 2) {
// 	    transprob = forwardProjecti(vec);
// 	  } else error("'solvetype' can by 1 or 2");
// 	  vector<Float> tmp = xi.row(tt-1).cwiseProduct(transprob);
// 	  j = 0;
// 	  thismax = tmp(0);
// 	  for(int i=1; i < n; i++){
// 	    if (tmp(i) > thismax){
// 	      thismax = tmp(i);
// 	      j = i;
// 	    }
// 	  }
// 	  track(tt-1) = j + 1;
// 	}
//       }  

//       return ans;
//     }
    
//   };

  
//   // ****** How to use it in TMB:
//   // 1. Create an evaluator 'eval' for previous class
//   template<class Float>
//   Float eval(Float Dx, Float Dy) {
//     filter_t<Float> f;
//     return f.eval(Dx, Dy);
//   }
//   // 2. Run 'eval' through tiny_ad and obtain an atomic function
//   //    'func'.  The '11' tells tiny_ad that we need
//   //    derivatives wrt. both Dx and Dy.
//   TMB_BIND_ATOMIC(func, 11, eval(x[0], x[1]))
//   // 3. Create a more user-friendly version ('func' takes vector
//   //    arguments and there's a final invisible argument that
//   //    corresponds to the derivative order)
//   template<class Type>
//   Type hmm_nll(Type Dx, Type Dy) {
//     vector<Type> args(3); // Last index reserved for derivative order
//     args << Dx, Dy, 0;
//     return shmm::func(CppAD::vector<Type>(args))[0];
//   }
// }

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
    um = atomic::matmul(EXPGx, um);
    um = atomic::matmul(um, EXPGy);
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
      ## - 4D array 'map' indexed by (x, y, time, datatype)
      ## - 2D array 'obs' indexed by (time, datatype)
      ## - 2D array 'con' indexed by (x, y)
    */
    DATA_ARRAY(map);
    DATA_ARRAY(obs);
    DATA_VECTOR(con);
    PARAMETER_VECTOR(logsd); // obs
    vector<Type> sd = exp(logsd);
    array<Type> keep(map.dim);
    for (int i=0; i<map.size(); i++) {
      if (R_FINITE(asDouble(map[i]))) {
        keep[i] = 1;
      } else {
        keep[i] = 0; map[i] = 0;
      }
    }
    /* Manually calculate dnorm to keep memory low:

       The default implementation does this:

       x = (x - mean) / sd;
       logres = -log(Type(sqrt(2*M_PI)) * sd) - Type(.5) * x * x;

       We need the special case where mean=0 and the x vector is constant:

       -log(Type(sqrt(2*M_PI)) - log(sd) - Type(.5) * x * x / sd;
       (-.5 * x * x) / sd - log(sd) - log(sqrt(2*M_PI)

       in addition we need to multiply by 'keep':

       (-.5 * x * x * keep) / sd - ( log(sd) + log(sqrt(2*M_PI) ) * keep

       We notice that '(-.5 * x * x * keep)' is constant, so can be
       calculated without adding anything to the tape.

    */


    /* Create partially calculated dnorm array (the constant):  */
    array<Type> res(map.dim);
    for (int k=0; k<map.cols(); k++) { // loop across k=datatype
      array<Type> submap = map.col(k); // (x, y, time)
      for (int j=0; j<submap.cols(); j++) { // loop across j=time
        res.col(k).col(j) = obs(j, k) - submap.col(j); // x-mu
        res.col(k).col(j) = res.col(k).col(j) * res.col(k).col(j); // (x-mu)*(x-mu)
        res.col(k).col(j) = Type(-.5)*res.col(k).col(j); // -.5*(x-mu)*(x-mu)
        res.col(k).col(j) = res.col(k).col(j) * keep.col(k).col(j);
      }
    }
    /* Sum up terms using vectorized operations */
    typedef vectorize::vector<Type> vec;
    vec tot;
    for (int k=0; k<map.cols(); k++) { // loop across k=datatype
      vec tmp = vec(res.col(k)) / sd[k] + ( -logsd[k] - log(sqrt(2.*M_PI)) ) * vec(keep.col(k));
      if (k==0)
        tot = tmp;
      else
        tot += tmp;
    }
    tot = exp(tot); // (x,y,time)
    // TODO: Add constraint
    vector<Type> vtot = tot;
    datlik.array() = vtot;
    vector<int> d = map.dim;
    datlik.resize(d[0]*d[1], d[2]);
    //tot = tot * vec();
  }

  // Matrix exponential

  // Sparse case:
  if (false) {
    Eigen::SparseMatrix<Type> G = Dx*Sew + Dy*Sns;  // Make generator
    Eigen::SparseMatrix<Type> Gdt = G * dt;
    sparse_matrix_exponential::expQ<Type> ExpGdt(Gdt);
  }

  // Dense separable case
  matrix<Type> Gx = (Dx*dt)*Sx;
  matrix<Type> Gy = (Dy*dt)*Sy;
  separable_matrix_exponential<Type> ExpGdt(Gx, Gy);
  
  // ==================================== HMM FILTER
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
    Rcout << "t=" << t << "\n";
    // Time update using uniformization algorithm
    vector<Type> svec = phi.col(t-1);
    //vector<Type> predtmp = svec;
    vector<Type> predtmp = ExpGdt(svec);
    pred.col(t) = predtmp; // Store prediction
    if (iobs(t) > 0) {
      // Data update
      int ind = (iobs(t)-1);
      vector<Type> d = datlik.col(ind);
      vectorize::vector<Type> post =
        vectorize::vector<Type>(d) * vectorize::vector<Type>(predtmp);
      //vectorize::vector<Type> post = pred.col(t).cwiseProduct(datlik.col(ind)); // Element-wise product
      psi(ind - 1) = post.sum();
      vectorize::vector<Type> post_normalized = (post / (psi(ind - 1) + 1e-20)); // Add small value to avoid division by zero
      vector<Type> qw = post_normalized;
      phi.col(t) = qw;
    } else {
      // No data update
      phi.col(t) = pred.col(t);
    }
  }
  // Negative log likelihood
  Type ans = -sum(log(psi));

  /* With this construct 'dosmoo' flag becomes redundant. The branch
     is only entered by obj$report(). */
  if(isDouble<Type>::value) {
    smoo.resize(n, ns);
    // Smoothing loop
    smoo.col(ns-1) = phi.col(ns-1);
    vector<Type> rho = datlik.col(ns-1);
    for(int t=1; t < ns; t++) {
      int tt = ns - t - 1;
      vector<Type> tmp = ExpGdt(rho); // FIXE: Transpose G if not symmetric!
      vector<Type> phat = tmp * vector<Type>(phi.col(tt).array()); phat = phat / sum(phat);
      smoo.col(tt) = phat;
      rho = tmp * vector<Type>(datlik.col(tt).array()); rho = rho / sum(rho);
    }
      // // Time update using uniformization algorithm
      // vector<Type> predcol = pred.col(tt);
      // for (int i=0; i < n; i++) {
      //   predcol(0, i) += 1e-10;
      // }
      
      // matrix<Float> ratio = smoo.col(tt).cwiseQuotient(predcol);
      //     //matrix<double> asd(1, n);
      //     //for (int i=0; i < n; i++){
      //     //asd(0, i) = std::isnan(asDouble(ratio(0, i)));
      //     //cout 
      //     //}
      //     //matrix<double> asd = std::isnan(asDouble(ratio));
      //     matrix<Float> ratiotmp = ratio;
      //     if (solvetype == 1){
      //       ratiotmp = forwardProject(ratio);
      //     } else if (solvetype == 2) {
      //       ratiotmp = forwardProjecti(ratio);
      //     } else error("'solvetype' can by 1 or 2");
      //     matrix<Float> post = phi.col(tt-1).cwiseProduct(ratiotmp);
      //     post = post / (post.sum() + 1e-20);
      //     smoo.col(tt-1) = post;
      //   }	
      // }


    
  //   shmm::filter_t<double> filter;
  //   filter.eval(asDouble(Dx), asDouble(Dy), true /* do_smoothing */, doviterbi == 1 /* do_viterbi */ );
  //   // Store a subset of distribution for output
  //   int nobs = shmm::datlik.cols();
  //   int n = shmm::datlik.rows();
  //   int ns = shmm::ns;
  //   vector<int> iobs = shmm::iobs;
  //   matrix<double> smooout(nobs, n);
  //   matrix<double> phiout(nobs, n);
  //   matrix<double> predout(nobs, n);
  //   for(int t=0; t<ns; t++){
  //     if (iobs(t) > 0){
  //       int ind = (iobs(t)-1);
  //       smooout.col(ind) = filter.smoo.col(t);
  //       phiout.col(ind) = filter.phi.col(t);
  //       predout.col(ind) = filter.pred.col(t);
  //     }
  //   }
  //   vector<double> psi = filter.psi;
  //   vector<int> trackout = filter.track;
  //   // Reports
    matrix<Type> predout = pred;
    matrix<Type> phiout = phi;
    matrix<Type> psiout = psi;
    matrix<Type> smooout = smoo;
    REPORT(predout);
    REPORT(phiout);
    REPORT(psi);
    REPORT(smooout); 
  //   //REPORT(xi);
  //   REPORT(trackout);
  }
  
  return ans;

}

