#include <Rcpp.h>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <limits>

namespace {
using Eigen::MatrixXd;
using Eigen::VectorXd;

double softplus(double x) { return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x))); }
double sigmoid(double x) {
    if (x >= 0) return 1.0 / (1.0 + std::exp(-x));
    double z = std::exp(x); return z / (1.0 + z);
}

struct Evaluation {
    double value;
    VectorXd gradient, curvature, eta;
};

// All arithmetic here uses mean weighted loss. The R interface documents the
// equivalent summed-loss lambda and ridge parameters.
Evaluation evaluate(const MatrixXd& A, const MatrixXd& R, const VectorXd& y,
                    const VectorXd& w, const VectorXd& offset, const VectorXd& b,
                    const VectorXd& penalty, double mu, bool barrier, bool gaussian = false) {
    Evaluation e;
    e.eta = A*b + offset;
    VectorXd residual(y.size()); e.curvature.resize(y.size());
    // Compensated accumulation keeps the objective accurate enough for
    // near-converged Newton line searches on large Bernoulli samples.
    // Otherwise O(n) rounding can exceed the true Newton decrease and
    // repeatedly reject a valid step while stationarity still exceeds tol.
    // Evaluate the ridge as ||R b||^2 / 2, and its gradient as R'(R b).
    // Expanding this into b'(R'R)b loses small differences between nearly
    // parallel spline columns and can make a valid Newton step look uphill.
    VectorXd ridge_residual=R*b;
    double total=0.5*ridge_residual.squaredNorm(), correction=0;
    auto add=[&](double value) {
        double next=total+value;
        correction+=std::abs(total)>=std::abs(value) ? (total-next)+value : (value-next)+total;
        total=next;
    };
    for (int i=0; i<y.size(); ++i) {
        if (gaussian) {
            double difference = e.eta[i] - y[i];
            add(0.5 * w[i] * difference * difference);
            residual[i] = w[i] * difference;
            e.curvature[i] = w[i];
            continue;
        }
        double p = sigmoid(e.eta[i]);
        // This form avoids cancellation for y=1 and very large positive eta.
        add(w[i] * ((1-y[i])*softplus(e.eta[i]) + y[i]*softplus(-e.eta[i])));
        residual[i] = w[i] * (p-y[i]);
        double z = std::exp(-std::abs(e.eta[i]));
        e.curvature[i] = w[i] * z / ((1+z)*(1+z));
    }
    e.gradient = A.transpose()*residual + R.transpose()*ridge_residual;
    if (barrier) {
        for (int j=0; j<b.size(); ++j) if (penalty[j]>0) {
            double s = std::hypot(mu, penalty[j]*b[j]);
            add(s - mu*std::log(mu+s));
            e.gradient[j] += penalty[j] * (penalty[j]*b[j]/(mu+s));
        }
    }
    e.value=total+correction;
    return e;
}

// Compare a trial with the current iterate without subtracting two rounded
// objectives or recomputing two large, cancelling A*b products. In a late
// barrier stage, a meaningful dual correction can have a likelihood decrease
// far below the rounding error of either absolute objective. Use the actual
// representable coefficient change, including any rounding in b + step*d.
double objective_change(const MatrixXd& A, const MatrixXd& R,
                        const VectorXd& y, const VectorXd& w,
                        const VectorXd& b, const VectorXd& change,
                        const VectorXd& eta, const VectorXd& penalty,
                        double mu, bool gaussian) {
    VectorXd h=A*change, rb=R*b, rh=R*change;
    double total=rb.dot(rh)+0.5*rh.squaredNorm(), correction=0;
    auto add=[&](double value) {
        double next=total+value;
        correction+=std::abs(total)>=std::abs(value) ? (total-next)+value : (value-next)+total;
        total=next;
    };
    for(int i=0;i<y.size();++i) {
        if(gaussian) add(w[i]*h[i]*(eta[i]-y[i]+0.5*h[i]));
        else if(std::abs(h[i])<=1) {
            add(w[i]*((1-y[i])*std::log1p(sigmoid(eta[i])*std::expm1(h[i]))+
                y[i]*std::log1p(sigmoid(-eta[i])*std::expm1(-h[i]))));
        } else {
            add(w[i]*((1-y[i])*(softplus(eta[i]+h[i])-softplus(eta[i]))+
                y[i]*(softplus(-eta[i]-h[i])-softplus(-eta[i]))));
        }
    }
    for(int j=0;j<b.size();++j) if(penalty[j]>0) {
        double a=penalty[j]*b[j], delta=penalty[j]*change[j];
        double s=std::hypot(mu,a), next=std::hypot(mu,a+delta);
        double ds=delta*((a+delta)/(s+next)+a/(s+next));
        // A large decrease can round ds/(mu+s) to -1. Use the positive
        // endpoints there; log1p is needed only for the small-change case.
        double logratio=ds < -0.5*(mu+s) ? std::log(mu+next)-std::log(mu+s) : std::log1p(ds/(mu+s));
        add(ds-mu*logratio);
    }
    return total+correction;
}

// Search only directions with zero L1 and ridge penalty. A verified logistic
// recession direction establishes nonexistence of a finite minimizer. Failure
// to find one is not a proof that no separation exists.
VectorXd recession_direction(const MatrixXd& A, const MatrixXd& R,
                            const VectorXd& y, const VectorXd& w,
                            const VectorXd& off, const VectorXd& penalty) {
    std::vector<int> free;
    for(int j=0;j<A.cols();++j) if(penalty[j]==0 && R.col(j).squaredNorm()==0) free.push_back(j);
    if(free.empty()) return VectorXd::Zero(A.cols());
    MatrixXd B(A.rows(),free.size());
    for(int j=0;j<B.cols();++j) B.col(j)=A.col(free[j]);
    VectorXd b=VectorXd::Zero(B.cols()),zero=VectorXd::Zero(B.cols());
    MatrixXd zeroQ=MatrixXd::Zero(B.cols(),B.cols());
    for(int it=0;it<60;++it) {
        Evaluation e=evaluate(B,zeroQ,y,w,off,b,zero,0,false);
        MatrixXd h=B.transpose()*(B.array().colwise()*e.curvature.array()).matrix();
        h.diagonal().array()+=1e-14;
        VectorXd step=h.ldlt().solve(-e.gradient);
        if(!step.allFinite()) break;
        VectorXd v=B*step;
        double size=v.cwiseAbs().maxCoeff(),positive=0,violation=0;
        if(size>1e-8) {
            for(int i=0;i<y.size();++i) {
                if(y[i]==0 || y[i]==1) {
                    double margin=(2*y[i]-1)*v[i]/size;
                    violation=std::max(violation,-margin);
                    positive=std::max(positive,margin);
                } else violation=std::max(violation,std::abs(v[i])/size);
            }
            if(violation<=1e-12 && positive>1e-6) {
                VectorXd out=VectorXd::Zero(A.cols());
                for(int j=0;j<B.cols();++j) out[free[j]]=step[j]/size;
                return out;
            }
        }
        double length=1;
        while(length>1e-8 && evaluate(B,zeroQ,y,w,off,b+length*step,zero,0,false).value>e.value)
            length*=0.5;
        b+=length*step;
        if(e.gradient.lpNorm<Eigen::Infinity>()<1e-13) break;
    }
    return VectorXd::Zero(A.cols());
}
}

// Shared loss-family implementation; exported wrappers keep the binary API stable.
Rcpp::List univariate_tf_native(Rcpp::NumericMatrix design, Rcpp::NumericMatrix ridge_factor,
                             Rcpp::NumericVector response, Rcpp::NumericVector weights,
                             Rcpp::NumericVector offsets, Rcpp::NumericVector penalty_scale,
                             Rcpp::NumericVector lambdas, double tolerance,
                             double gap_tolerance, int max_iterations, bool gaussian) {
    const int n=design.nrow(), p=design.ncol();
    Eigen::Map<MatrixXd> fullA(design.begin(), n, p);
    Eigen::Map<MatrixXd> fullR(ridge_factor.begin(), ridge_factor.nrow(), p);
    Eigen::Map<VectorXd> y(response.begin(),n), raww(weights.begin(),n), off(offsets.begin(),n);
    Eigen::Map<VectorXd> scale(penalty_scale.begin(),p);
    double mass=raww.sum(); VectorXd w=raww/mass;
    VectorXd warm=VectorXd::Zero(p);
    double prevalence=w.dot(y);
    warm[0]=(gaussian ? prevalence : std::log(prevalence/(1-prevalence)))-w.dot(off);
    Rcpp::List fits(lambdas.size());
    for (int li=0; li<lambdas.size(); ++li) {
        Rcpp::checkUserInterrupt();
        const double lambda=lambdas[li];
        std::vector<int> keep;
        for (int j=0;j<p;++j) if (std::isfinite(lambda) || scale[j]==0) keep.push_back(j);
        int q=keep.size(), np=0;
        MatrixXd A(n,q),R(fullR.rows(),q); VectorXd b(q),penalty(q);
        for (int j=0;j<q;++j) {
            A.col(j)=fullA.col(keep[j]); b[j]=warm[keep[j]];
            penalty[j]=std::isfinite(lambda) ? lambda*scale[keep[j]]/mass : 0;
            if (penalty[j]>0) ++np;
            R.col(j)=fullR.col(keep[j])/std::sqrt(mass);
        }
        MatrixXd Q=R.transpose()*R;
        double target_mu=gap_tolerance/std::max(1,np);
        double mu=np ? std::max(0.01,target_mu) : 0;
        double stationarity=R_PosInf, complementarity=R_PosInf;
        double objective=R_PosInf, max_eta=0, max_damping=0;
        int iterations=0, stages=0, backtracks=0;
        bool converged=false;
        std::string status="iteration_limit";
        VectorXd recession = gaussian ? VectorXd::Zero(q).eval() : recession_direction(A,R,y,w,off,penalty);
        bool separated=recession.cwiseAbs().maxCoeff()>0;
        if(separated) status="separation_detected";
        while(!separated && iterations<max_iterations) {
            ++iterations;
            if (iterations%10==0) Rcpp::checkUserInterrupt();
            Evaluation e=evaluate(A,R,y,w,off,b,penalty,mu,np>0,gaussian);
            if (!std::isfinite(e.value) || !e.gradient.allFinite() || !b.allFinite()) {
                status="nonfinite_arithmetic"; break;
            }
            stationarity=e.gradient.lpNorm<Eigen::Infinity>();
            complementarity=0;
            for(int j=0;j<q;++j) if(penalty[j]>0) {
                double d=penalty[j]*(penalty[j]*b[j]/(mu+std::hypot(mu,penalty[j]*b[j])));
                complementarity+=penalty[j]*std::abs(b[j])-d*b[j];
            }
            if(stationarity<=tolerance && complementarity<=gap_tolerance) {
                converged=true; status="converged"; break;
            }
            if (np && stationarity<=std::max(tolerance,mu*0.05) && mu>target_mu) {
                mu=std::max(target_mu,mu*0.1); ++stages; continue;
            }
            MatrixXd weighted=A.array().colwise()*e.curvature.array();
            MatrixXd H=A.transpose()*weighted+Q;
            for(int j=0;j<q;++j) if(penalty[j]>0) {
                double s=std::hypot(mu,penalty[j]*b[j]);
                H(j,j)+=penalty[j]*(penalty[j]/(mu+s))*(mu/s);
            }
            VectorXd hs=H.diagonal().cwiseMax(1e-30).cwiseSqrt().cwiseInverse();
            MatrixXd scaled=hs.asDiagonal()*H*hs.asDiagonal();
            VectorXd direction(q); bool solved=false;
            double damping=0;
            for(int attempt=0;attempt<10;++attempt) {
                MatrixXd system=scaled;
                system.diagonal().array()+=damping;
                Eigen::LDLT<MatrixXd> factor(system);
                if(factor.info()==Eigen::Success && factor.isPositive()) {
                    direction=hs.asDiagonal()*factor.solve(-(hs.array()*e.gradient.array()).matrix());
                    if(direction.allFinite() && e.gradient.dot(direction)<0) { solved=true; break; }
                }
                damping=damping==0 ? 1e-12 : damping*10;
            }
            max_damping=std::max(max_damping,damping);
            if(!solved) { status="singular_newton_system"; break; }
            double step=1, slope=e.gradient.dot(direction);
            bool accepted=false;
            for(int ls=0;ls<60;++ls) {
                VectorXd candidate=b+step*direction;
                double change=objective_change(A,R,y,w,b,candidate-b,e.eta,penalty,mu,gaussian);
                // Allow only rounding-level slack, and still require the KKT
                // checks above before declaring convergence.
                if(std::isfinite(change) && change<=1e-4*step*slope+1e-15*std::max(1.0,std::abs(e.value))) {
                    b=candidate; accepted=true; break;
                }
                step*=0.5; ++backtracks;
            }
            if(!accepted) { status="line_search_failed"; break; }
            if(b.cwiseAbs().maxCoeff()>1e12) { status="unbounded_or_ill_conditioned"; break; }
        }
        Evaluation final=evaluate(A,R,y,w,off,b,penalty,mu,false,gaussian);
        objective=final.value;
        VectorXd dual=VectorXd::Zero(q);
        complementarity=0;
        for(int j=0;j<q;++j) if(penalty[j]>0) {
            dual[j]=penalty[j]*b[j]/(mu+std::hypot(mu,penalty[j]*b[j]));
            final.gradient[j]+=penalty[j]*dual[j];
            complementarity+=penalty[j]*(std::abs(b[j])-dual[j]*b[j]);
            objective+=penalty[j]*std::abs(b[j]);
        }
        stationarity=final.gradient.lpNorm<Eigen::Infinity>();
        max_eta=final.eta.cwiseAbs().maxCoeff();
        converged=converged && std::isfinite(objective) && final.gradient.allFinite() &&
            dual.allFinite() && dual.cwiseAbs().maxCoeff()<=1 &&
            stationarity<=tolerance && complementarity<=gap_tolerance;
        VectorXd probability(n);
        for(int i=0;i<n;++i) probability[i]=gaussian ? final.eta[i] : sigmoid(final.eta[i]);
        VectorXd full=VectorXd::Zero(p),fulldual=VectorXd::Zero(p);
        for(int j=0;j<q;++j) { full[keep[j]]=b[j]; fulldual[keep[j]]=dual[j]; }
        if(converged) warm=full;
        else if(status=="converged") status="diagnostic_check_failed";
        Rcpp::NumericVector coeff(p),duals(p),eta(n),prob(n),witness(p);
        std::copy(full.data(),full.data()+p,coeff.begin());
        std::copy(fulldual.data(),fulldual.data()+p,duals.begin());
        std::copy(final.eta.data(),final.eta.data()+n,eta.begin());
        std::copy(probability.data(),probability.data()+n,prob.begin());
        for(int j=0;j<q;++j) witness[keep[j]]=recession[j];
        fits[li]=Rcpp::List::create(Rcpp::_["coefficients"]=coeff,
            Rcpp::_["dual"]=duals, Rcpp::_["linear.predictors"]=eta,
            Rcpp::_["separation.direction"]=witness,
            Rcpp::_["fitted.values"]=prob,Rcpp::_["objective"]=objective*mass,
            Rcpp::_["converged"]=converged,Rcpp::_["status"]=status,
            Rcpp::_["stationarity"]=stationarity,
            Rcpp::_["complementarity"]=complementarity,
            Rcpp::_["dual.feasibility"]=std::max(0.0,dual.cwiseAbs().maxCoeff()-1),
            Rcpp::_["iterations"]=iterations,Rcpp::_["barrier.stages"]=stages,
            Rcpp::_["barrier.parameter"]=mu,Rcpp::_["backtracks"]=backtracks,
            Rcpp::_["max.damping"]=max_damping,Rcpp::_["max.abs.eta"]=max_eta);
    }
    return fits;
}

// [[Rcpp::export]]
Rcpp::List logistic_tf_native(Rcpp::NumericMatrix design, Rcpp::NumericMatrix ridge_factor,
    Rcpp::NumericVector response, Rcpp::NumericVector weights, Rcpp::NumericVector offsets,
    Rcpp::NumericVector penalty_scale, Rcpp::NumericVector lambdas,
    double tolerance, double gap_tolerance, int max_iterations) {
    return univariate_tf_native(design, ridge_factor, response, weights, offsets,
        penalty_scale, lambdas, tolerance, gap_tolerance, max_iterations, false);
}

// [[Rcpp::export]]
Rcpp::List gaussian_tf_native(Rcpp::NumericMatrix design, Rcpp::NumericMatrix ridge_factor,
    Rcpp::NumericVector response, Rcpp::NumericVector weights, Rcpp::NumericVector offsets,
    Rcpp::NumericVector penalty_scale, Rcpp::NumericVector lambdas,
    double tolerance, double gap_tolerance, int max_iterations) {
    return univariate_tf_native(design, ridge_factor, response, weights, offsets,
        penalty_scale, lambdas, tolerance, gap_tolerance, max_iterations, true);
}
