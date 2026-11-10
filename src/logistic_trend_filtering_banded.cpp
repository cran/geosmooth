#include <Rcpp.h>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <limits>
#include <vector>

// Small unpenalized quadratic initialization uses the existing, independently
// validated solver. No observation-by-full-knot dense matrix is constructed.
Rcpp::List logistic_tf_native(Rcpp::NumericMatrix, Rcpp::NumericMatrix,
    Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector,
    Rcpp::NumericVector, Rcpp::NumericVector, double, double, int);

Rcpp::List gaussian_tf_native(Rcpp::NumericMatrix, Rcpp::NumericMatrix,
    Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector,
    Rcpp::NumericVector, Rcpp::NumericVector, double, double, int);

namespace ltf_banded {
using Real = long double;
using Vec = Eigen::Matrix<Real,Eigen::Dynamic,1>;
using Mat = Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic>;
using Row = std::vector<std::pair<int,Real>>;
double prob(double a) { if(a>=0) return 1/(1+std::exp(-a)); double e=std::exp(a); return e/(1+e); }
Real sp(Real a) { return std::max(a,Real(0))+std::log1p(std::exp(-std::abs(a))); }
Real infnorm(const Vec& a) { return a.size() ? a.cwiseAbs().maxCoeff() : 0; }
struct Sum {
    Real total=0,correction=0;
    void add(Real value) {
        Real next=total+value;
        correction+=std::abs(total)>=std::abs(value) ? (total-next)+value : (value-next)+total;
        total=next;
    }
    Real value() const {return total+correction;}
};

// Shifted suffix moments avoid subtracting nearly equal powers near x=1.
Mat moments(const std::vector<double>& x, const Vec& v, const std::vector<double>& knots) {
    static const int choose[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},{1,3,3,1,0},{1,4,6,4,1}};
    std::vector<int> order(x.size()); std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){return x[a]<x[b];});
    Mat out(knots.size(),5); std::array<Real,5> m{};
    int pos=order.size()-1; Real origin=knots.back();
    for(int j=int(knots.size())-1;j>=0;--j) {
        Real delta=origin-Real(knots[j]); std::array<Real,5> next{};
        for(int k=0;k<5;++k) for(int l=0;l<=k;++l) next[k]+=choose[k][l]*std::pow(delta,k-l)*m[l];
        m=next; origin=knots[j];
        while(pos>=0 && x[order[pos]]>=knots[j]) {
            int i=order[pos--]; Real t=Real(x[i])-origin,power=1;
            for(int k=0;k<5;++k) {m[k]+=v[i]*power;power*=t;}
        }
        for(int k=0;k<5;++k) out(j,k)=m[k];
    }
    return out;
}

struct Geometry {
    std::vector<double> x,knots;
    Vec center,scale,w;
    Geometry(const std::vector<double>& xx,const std::vector<double>& kk,const Vec& ww):x(xx),knots(kk),w(ww) {
        int m=knots.size(); center=Vec::Zero(m-1);scale=Vec::Zero(m-1);
        Mat s=moments(x,w,knots);
        for(int j=0;j<2;++j) {
            Real mean=0,square=0;
            for(int i=0;i<int(x.size());++i) {Real a=std::pow(Real(x[i]),j+1);mean+=w[i]*a;square+=w[i]*a*a;}
            center[j]=mean;scale[j]=std::sqrt(std::max(Real(0),square-mean*mean));
        }
        for(int j=2;j<m-1;++j) {
            Real h=Real(knots[j])-knots[j-1];
            Real mean=(s(j,2)+h*s(j,1))/2;
            Real square=(s(j,4)+2*h*s(j,3)+h*h*s(j,2))/4;
            center[j]=mean;scale[j]=std::sqrt(std::max(Real(0),square-mean*mean));
        }
        if(!center.allFinite() || !scale.allFinite() || (scale.array()<=0).any())
            Rcpp::stop("Degenerate structured falling-factorial geometry.");
    }
    Vec gradient(const Vec& r) const {
        Mat s=moments(x,r,knots); Vec g(knots.size());g.setZero();
        for(int i=0;i<int(x.size());++i) {g[0]+=r[i];g[1]+=r[i]*x[i];g[2]+=r[i]*x[i]*x[i];}
        for(int j=2;j<int(knots.size())-1;++j)
            g[j+1]=(s(j,2)+(Real(knots[j])-knots[j-1])*s(j,1))/2;
        return g;
    }
};

// A symmetric matrix stored as its diagonal and three lower diagonals.
struct Band {
    int n; std::vector<std::array<Real,4>> a,l; Vec d,sc;
    explicit Band(int size):n(size),a(size),l(size),d(size),sc(size) {}
    void add(int i,int j,Real v) {
        if(i<j) std::swap(i,j);
        if(i-j>3) {if(std::abs(v)>1e-25) Rcpp::stop("Structured basis exceeded its bandwidth.");return;}
        a[i][i-j]+=v;
    }
    bool factor(Real damping=0) {
        for(int i=0;i<n;++i) {sc[i]=std::sqrt(std::max(a[i][0],Real(1e-30)));l[i].fill(0);}
        for(int i=0;i<n;++i) {
            for(int j=std::max(0,i-3);j<i;++j) {
                Real v=a[i][i-j]/(sc[i]*sc[j]);
                for(int k=std::max(0,i-3);k<j;++k) if(j-k<=3) v-=l[i][i-k]*d[k]*l[j][j-k];
                l[i][i-j]=v/d[j];
            }
            Real v=a[i][0]/(sc[i]*sc[i])+damping;
            for(int k=std::max(0,i-3);k<i;++k) v-=l[i][i-k]*l[i][i-k]*d[k];
            if(!(v>0) || !std::isfinite(v)) return false;
            d[i]=v;
        }
        return true;
    }
    Vec solve(const Vec& rhs) const {
        Vec z=rhs.cwiseQuotient(sc);
        for(int i=0;i<n;++i) for(int k=std::max(0,i-3);k<i;++k) z[i]-=l[i][i-k]*z[k];
        z=z.cwiseQuotient(d);
        for(int i=n-1;i>=0;--i) for(int k=i+1;k<std::min(n,i+4);++k) z[i]-=l[k][k-i]*z[k];
        return z.cwiseQuotient(sc);
    }
};

// Compact columns are exact linear combinations of four consecutive
// falling-factorial functions whose polynomial tails cancel. Hence each
// observation touches at most three columns and the Hessian is banded.
struct Basis {
    std::vector<double> knots;
    std::vector<int> active;
    std::vector<std::array<Real,4>> transform,polynomial0,polynomial1,polynomial2;
    std::vector<Real> origin,span;
    std::vector<Row> rows;
    int d;
    Real start(int j) const {return j<3 ? 0 : knots[active[j-3]];}
    std::array<Real,3> polynomial(int j,Real s,Real h) const {
        if(j==0) return {1,0,0};
        if(j==1) return {s,h,0};
        if(j==2) return {s*s,2*s*h,h*h};
        Real a=knots[active[j-3]-1],b=knots[active[j-3]];
        return {(s-a)*(s-b)/2,h*(2*s-a-b)/2,h*h/2};
    }
    Real eval(int j,Real u) const {
        if(u<start(j) || (j<d-3 && u>=start(j+3))) return 0;
        Real t=(u-origin[j])/span[j],v=0;
        int count=std::min(4,d-j);
        for(int k=0;k<count;++k) if(j+k<3 || u>start(j+k))
            v+=transform[j][k]*(polynomial0[j][k]+t*(polynomial1[j][k]+t*polynomial2[j][k]));
        return v;
    }
    Basis(const std::vector<double>& kk,const std::vector<int>& aa,const std::vector<double>& x):knots(kk),active(aa),d(aa.size()+3) {
        transform.resize(d);polynomial0.resize(d);polynomial1.resize(d);polynomial2.resize(d);
        origin.resize(d);span.resize(d);rows.resize(x.size());
        std::vector<int> order(x.size());std::iota(order.begin(),order.end(),0);
        std::stable_sort(order.begin(),order.end(),[&](int a,int b){return x[a]<x[b];});
        for(int j=0;j<d;++j) {
            origin[j]=start(j);span[j]=(j<d-3 ? start(j+3) : Real(1))-origin[j];
            if(!(span[j]>0)) Rcpp::stop("Degenerate compact basis interval.");
            int count=std::min(4,d-j);
            for(int k=0;k<count;++k) {auto p=polynomial(j+k,origin[j],span[j]);polynomial0[j][k]=p[0];polynomial1[j][k]=p[1];polynomial2[j][k]=p[2];}
            if(j<d-3) {
                Eigen::Matrix<Real,3,4> P;
                std::array<Real,4> norm;
                for(int k=0;k<4;++k) {
                    norm[k]=std::max({std::abs(polynomial0[j][k]),std::abs(polynomial1[j][k]),std::abs(polynomial2[j][k])});
                    if(norm[k]==0) Rcpp::stop("Degenerate compact basis column.");
                    P.col(k)<<polynomial0[j][k]/norm[k],polynomial1[j][k]/norm[k],polynomial2[j][k]/norm[k];
                }
                for(int k=0;k<4;++k) {
                    Eigen::Matrix<Real,3,3> minor;int at=0;
                    for(int h=0;h<4;++h) if(h!=k) minor.col(at++)=P.col(h);
                    transform[j][k]=(k%2 ? -1 : 1)*minor.determinant()/norm[k];
                }
                Real largest=0;for(Real v:transform[j]) largest=std::max(largest,std::abs(v));
                if(!(largest>0)) Rcpp::stop("Singular compact basis transform.");
                for(Real& v:transform[j]) v/=largest;
            } else transform[j][0]=1;
            Real peak=0;
            std::vector<Real> breaks{origin[j],origin[j]+span[j]};
            for(int k=1;k<count;++k) breaks.push_back(start(j+k));
            std::sort(breaks.begin(),breaks.end());
            for(int k=0;k<int(breaks.size());++k) {
                peak=std::max(peak,std::abs(eval(j,breaks[k])));
                if(k) peak=std::max(peak,std::abs(eval(j,(breaks[k]+breaks[k-1])/2)));
            }
            if(!(peak>0) || !std::isfinite(peak)) Rcpp::stop("Unresolved compact basis normalization.");
            for(Real& v:transform[j]) v/=peak;
            auto first=std::lower_bound(order.begin(),order.end(),double(start(j)),[&](int i,double a){return x[i]<a;});
            for(auto it=first;it!=order.end();++it) {
                int i=*it;
                if(j<d-3 && x[i]>=start(j+3)) break;
                Real v=eval(j,x[i]);if(v!=0) rows[i].push_back({j,v});
            }
        }
    }
    Vec values(const Vec& c) const {
        Vec v(rows.size());v.setZero();
        for(int i=0;i<int(rows.size());++i) for(auto a:rows[i]) v[i]+=a.second*c[a.first];
        return v;
    }
    Vec transpose(const Vec& v) const {
        Vec out=Vec::Zero(d);for(int i=0;i<int(rows.size());++i) for(auto a:rows[i]) out[a.first]+=a.second*v[i];return out;
    }
    Vec coefficients(const Vec& c) const {
        Vec out=Vec::Zero(d);for(int j=0;j<d;++j) for(int k=0;k<std::min(4,d-j);++k) out[j+k]+=transform[j][k]*c[j];return out;
    }
    Vec transpose_transform(const Vec& v) const {
        Vec out=Vec::Zero(d);for(int j=0;j<d;++j) for(int k=0;k<std::min(4,d-j);++k) out[j]+=transform[j][k]*v[j+k];return out;
    }
    Vec project(const Vec& values,const Vec& weights) const {
        // The trailing truncated quadratics can be nearly parallel when active
        // knots cluster. Normal equations square that condition number and
        // can produce nonpositive LDL pivots even for a full-rank design.
        // Scaled, row-wise Givens QR solves the same weighted least-squares
        // projection without a Gram matrix, ridge, or column deletion.
        // Visit rows in nondecreasing right support so rotations cannot fill
        // outside the four stored upper diagonals (O(n+d) workspace).
        std::vector<std::array<Real,4>> r(d);
        Vec rhs=Vec::Zero(d),scale=Vec::Zero(d);
        for(int i=0;i<int(rows.size());++i) for(auto a:rows[i]) scale[a.first]+=weights[i]*a.second*a.second;
        for(int j=0;j<d;++j) {
            scale[j]=std::sqrt(scale[j]);
            if(!(scale[j]>0)) Rcpp::stop("Singular compact-basis projection.");
        }
        std::vector<int> order(rows.size());std::iota(order.begin(),order.end(),0);
        std::stable_sort(order.begin(),order.end(),[&](int a,int b) {
            int left=rows[a].empty() ? -1 : rows[a].back().first;
            int right=rows[b].empty() ? -1 : rows[b].back().first;
            return left<right;
        });
        for(int i:order) {
            if(rows[i].empty()) continue;
            int first=rows[i].front().first,last=rows[i].back().first;
            if(last-first>3) Rcpp::stop("Structured projection exceeded its bandwidth.");
            std::array<Real,4> v{};Real root=std::sqrt(weights[i]),target=root*values[i];
            for(auto a:rows[i]) v[a.first-first]=root*a.second/scale[a.first];
            for(int j=first;j<=last;++j) {
                Real entry=v[j-first];if(entry==0)continue;
                Real norm=std::hypot(r[j][0],entry),co=r[j][0]/norm,si=entry/norm;
                r[j][0]=norm;
                for(int k=j+1;k<=last;++k) {
                    Real before=r[j][k-j];
                    r[j][k-j]=co*before+si*v[k-first];
                    v[k-first]=-si*before+co*v[k-first];
                }
                Real before=rhs[j];rhs[j]=co*before+si*target;target=-si*before+co*target;
            }
        }
        Vec c=rhs;
        for(int j=d-1;j>=0;--j) {
            if(!(r[j][0]>0) || !std::isfinite(r[j][0])) Rcpp::stop("Singular compact-basis projection.");
            for(int k=j+1;k<std::min(d,j+4);++k)c[j]-=r[j][k-j]*c[k];
            c[j]/=r[j][0];
        }
        return c.cwiseQuotient(scale);
    }
};

struct Eval { Real objective,mean; Vec eta,probability,residual,curvegradient,curvature; };
Eval evaluate(const Vec& curve,const Mat& Z,const Vec& gamma,const Vec& offset,
              const Vec& y,const Vec& w,Real nr,Real cr,const Vec& theta,Real lambda,bool gaussian) {
    Eval e; int n=y.size();e.eta=curve+Z*gamma+offset;e.mean=w.dot(curve);
    e.probability.resize(n);e.residual.resize(n);e.curvature.resize(n);
    Sum objective;objective.add(nr*gamma.squaredNorm()/2);
    for(int i=0;i<n;++i) {
        if (gaussian) {
            e.probability[i]=e.eta[i];
            Real difference=e.eta[i]-y[i];
            e.residual[i]=w[i]*difference;e.curvature[i]=w[i];
            objective.add(w[i]*difference*difference/2);
        } else {
            e.probability[i]=prob(double(e.eta[i]));e.residual[i]=w[i]*(e.probability[i]-y[i]);
            Real z=std::exp(-std::abs(e.eta[i]));e.curvature[i]=w[i]*z/((1+z)*(1+z));
            objective.add(w[i]*((1-y[i])*sp(e.eta[i])+y[i]*sp(-e.eta[i])));
        }
        objective.add(cr*w[i]*(curve[i]-e.mean)*(curve[i]-e.mean)/2);
    }
    e.curvegradient=e.residual+cr*w.cwiseProduct((curve.array()-e.mean).matrix());
    if(std::isfinite(lambda)) for(int j=3;j<theta.size();++j) objective.add(lambda*std::abs(theta[j]));
    e.objective=objective.value();
    return e;
}

struct Diagnostic { Real stationarity=0,active=0,complementarity=0;int candidate=-1; Vec dual; };
Diagnostic diagnose(const Geometry& geometry,const Vec& g,const Vec& zg,const Vec& theta,
                    const std::vector<int>& active,const std::vector<int>& signs,Real lambda) {
    int m=geometry.knots.size();Diagnostic out;out.dual=Vec::Zero(m);
    out.stationarity=out.active=std::max(std::abs(g[0]),infnorm(zg));
    for(int j=1;j<3;++j) {
        Real v=std::abs((g[j]-geometry.center[j-1]*g[0])/geometry.scale[j-1]);
        out.stationarity=std::max(out.stationarity,v);out.active=std::max(out.active,v);
    }
    int at=0;Real worst=0;
    for(int j=2;j<m-1;++j) {
        bool included=at<int(active.size()) && active[at]==j;
        // The reported KKT equations use centered, scaled columns. Project
        // the dual in those same coordinates: projecting the raw moment g
        // instead leaves a spurious center*g[0] residual at inactive knots.
        // Near convergence that can admit a feasible knot and cycle forever
        // between adding it and removing it at a zero-length Newton step.
        Real centered=g[j+1]-geometry.center[j]*g[0];
        Real dual=0;
        if(lambda>0 && std::isfinite(lambda)) dual=included ? signs[at] : std::max(Real(-1),std::min(Real(1),-centered/lambda));
        out.dual[j+1]=dual;
        Real residual=std::isfinite(lambda) ? std::abs((centered+lambda*dual)/geometry.scale[j]) : 0;
        out.stationarity=std::max(out.stationarity,residual);
        if(included) {
            out.active=std::max(out.active,residual);
            out.complementarity+=lambda*(std::abs(theta[at+3])-dual*theta[at+3]); ++at;
        } else if(residual>worst) {worst=residual;out.candidate=j;}
    }
    return out;
}

bool direction(const Basis& B,const Eval& e,const Mat& Z,const Vec& w,
               const Vec& gb,const Vec& zg,Real nr,Real cr,Vec& dc,Vec& dg,Real& maxdamping) {
    int q=Z.cols(),s=q+(cr>0);Band h(B.d);Mat cross=Mat::Zero(B.d,s),small=Mat::Zero(s,s);
    Vec rhs=Vec::Zero(s);if(q) rhs.head(q)=-zg;
    Vec v=B.transpose(w);
    for(int i=0;i<int(B.rows.size());++i) {
        Real wt=e.curvature[i]+cr*w[i];
        for(auto a:B.rows[i]) {
            for(auto b:B.rows[i]) if(a.first>=b.first) h.add(a.first,b.first,wt*a.second*b.second);
            for(int j=0;j<q;++j) cross(a.first,j)+=a.second*e.curvature[i]*Z(i,j);
        }
        for(int j=0;j<q;++j) for(int k=0;k<=j;++k) small(j,k)+=e.curvature[i]*Z(i,j)*Z(i,k);
    }
    for(int j=0;j<q;++j) {small(j,j)+=nr;for(int k=0;k<j;++k) small(k,j)=small(j,k);}
    if(cr>0) {cross.col(q)=-cr*v;small(q,q)=cr;}
    Real damping=0;
    for(int attempt=0;attempt<10;++attempt) {
        if(h.factor(damping)) {
            Vec base=h.solve(-gb);Mat solved(B.d,s);
            for(int j=0;j<s;++j) solved.col(j)=h.solve(cross.col(j));
            Vec extra=Vec::Zero(s);
            bool okay=true;
            if(s) {
                Mat schur=small-cross.transpose()*solved;
                for(int j=0;j<s;++j) schur(j,j)+=damping*std::max(std::abs(small(j,j)),Real(1e-20));
                Eigen::LDLT<Mat> factor(schur);
                okay=factor.info()==Eigen::Success && factor.isPositive();
                if(okay) extra=factor.solve(rhs-cross.transpose()*base);
            }
            if(okay) {
                dc=base-solved*extra;dg=extra.head(q);
                if(dc.allFinite() && dg.allFinite() && gb.dot(dc)+zg.dot(dg)<0) {maxdamping=std::max(maxdamping,damping);return true;}
            }
        }
        damping=damping==0 ? 1e-12 : damping*10;
    }
    maxdamping=std::max(maxdamping,damping);return false;
}
Rcpp::NumericVector numeric(const Vec& a) {Rcpp::NumericVector out(a.size());for(int i=0;i<a.size();++i) out[i]=double(a[i]);return out;}
}

// [[Rcpp::export]]
Rcpp::List logistic_tf_structured_geometry(Rcpp::NumericVector u,Rcpp::NumericVector knots,Rcpp::NumericVector weights) {
    using namespace ltf_banded;Vec w(weights.size());Real mass=std::accumulate(weights.begin(),weights.end(),Real(0));
    for(int i=0;i<w.size();++i) w[i]=weights[i]/mass;
    Geometry g(Rcpp::as<std::vector<double>>(u),Rcpp::as<std::vector<double>>(knots),w);
    return Rcpp::List::create(Rcpp::_["center"]=numeric(g.center),Rcpp::_["scales"]=numeric(g.scale));
}

// [[Rcpp::export]]
Rcpp::NumericVector logistic_tf_banded_predict(Rcpp::NumericVector knots,Rcpp::IntegerVector active,
                                              Rcpp::NumericVector coefficients,Rcpp::NumericVector u) {
    using namespace ltf_banded;std::vector<double> kk=Rcpp::as<std::vector<double>>(knots),xx=Rcpp::as<std::vector<double>>(u);
    std::vector<int> aa=Rcpp::as<std::vector<int>>(active);Basis B(kk,aa,xx);
    Vec c(coefficients.size());for(int j=0;j<c.size();++j)c[j]=coefficients[j];return numeric(B.values(c));
}

Rcpp::List univariate_tf_banded_native(Rcpp::NumericVector u,Rcpp::NumericVector knots,Rcpp::NumericMatrix nuisance,
    Rcpp::NumericVector response,Rcpp::NumericVector weights,Rcpp::NumericVector offsets,Rcpp::NumericVector lambdas,
    double nuisance_ridge,double curve_ridge,double tolerance,double gap_tolerance,int max_iterations,bool gaussian) {
    using namespace ltf_banded;
    int n=u.size(),m=knots.size(),q=nuisance.ncol(),p=m+q;
    std::vector<double> xx=Rcpp::as<std::vector<double>>(u),kk=Rcpp::as<std::vector<double>>(knots);
    Vec y(n),w(n),off(n);Mat Z(n,q);Real mass=std::accumulate(weights.begin(),weights.end(),Real(0));
    for(int i=0;i<n;++i) {y[i]=response[i];w[i]=weights[i]/mass;off[i]=offsets[i];for(int j=0;j<q;++j)Z(i,j)=nuisance(i,j);}
    Geometry geometry(xx,kk,w);
    Rcpp::NumericMatrix initialA(n,q+3),initialR(q+(curve_ridge>0 ? n : 0),q+3);
    for(int i=0;i<n;++i) {
        initialA(i,0)=1;
        for(int j=0;j<q;++j)initialA(i,j+1)=double(Z(i,j));
        for(int j=0;j<2;++j)initialA(i,q+1+j)=double((std::pow(Real(xx[i]),j+1)-geometry.center[j])/geometry.scale[j]);
    }
    for(int j=0;j<q;++j)initialR(j,j+1)=std::sqrt(nuisance_ridge);
    if(curve_ridge>0)for(int j=0;j<2;++j)for(int i=0;i<n;++i)
        initialR(q+i,q+1+j)=std::sqrt(curve_ridge*double(w[i]))*initialA(i,q+1+j);
    auto initialize = gaussian ? gaussian_tf_native : logistic_tf_native;
    Rcpp::List initlist=initialize(initialA,initialR,response,weights,offsets,
        Rcpp::NumericVector(q+3),Rcpp::NumericVector::create(R_PosInf),tolerance,gap_tolerance,max_iterations);
    Rcpp::List init=initlist[0];Rcpp::NumericVector initb=init["coefficients"];
    Vec curve(n),gamma(q);
    for(int j=0;j<q;++j)gamma[j]=initb[j+1];
    for(int i=0;i<n;++i)curve[i]=initb[0]+initialA(i,q+1)*initb[q+1]+initialA(i,q+2)*initb[q+2];
    std::vector<int> active,signs;
    Rcpp::List fits(lambdas.size());
    for(int li=0;li<lambdas.size();++li) {
        Rcpp::checkUserInterrupt();Real lambda=std::isfinite(lambdas[li]) ? lambdas[li]/mass : R_PosInf;
        if(!std::isfinite(lambda)) {active.clear();signs.clear();}
        if(lambda==0) {active.clear();for(int j=2;j<m-1;++j)active.push_back(j);signs.assign(active.size(),1);}
        int iterations=0,backtracks=0,additions=0,removals=0;Real damping=0;
        Real laststep=0,lastslope=0;
        bool canfit=Rcpp::as<bool>(init["converged"]);
        std::string status=canfit ? "iteration_limit" : Rcpp::as<std::string>(init["status"]);
        int separated_knot=-1;
        if(!gaussian && lambda==0 && curve_ridge==0) {
            std::vector<Real> lo(m,1),hi(m,0);
            for(int i=0;i<n;++i) {
                int at=std::lower_bound(kk.begin(),kk.end(),xx[i])-kk.begin();
                lo[at]=std::min(lo[at],y[i]);hi[at]=std::max(hi[at],y[i]);
            }
            for(int j=0;j<m;++j) if(hi[j]==0 || lo[j]==1) {separated_knot=j;break;}
            if(separated_knot>=0) {canfit=false;status="separation_detected";}
        }
        bool converged=false,rebuild=true;Vec c,theta;Basis B(kk,active,xx);c=B.project(curve,w);
        Eval e=evaluate(B.values(c),Z,gamma,off,y,w,nuisance_ridge/mass,curve_ridge/mass,B.coefficients(c),lambda,gaussian);
        Diagnostic diag;
        for(;iterations<max_iterations && canfit;++iterations) {
            if(iterations%10==0)Rcpp::checkUserInterrupt();
            if(rebuild) {B=Basis(kk,active,xx);c=B.project(curve,w);rebuild=false;}
            theta=B.coefficients(c);curve=B.values(c);
            e=evaluate(curve,Z,gamma,off,y,w,nuisance_ridge/mass,curve_ridge/mass,theta,lambda,gaussian);
            Vec fullg=geometry.gradient(e.curvegradient),zg=Z.transpose()*e.residual+(nuisance_ridge/mass)*gamma;
            diag=diagnose(geometry,fullg,zg,theta,active,signs,lambda);
            if(!std::isfinite(e.objective) || !fullg.allFinite() || !zg.allFinite()) {status="nonfinite_arithmetic";break;}
            if(diag.stationarity<=tolerance && diag.complementarity<=gap_tolerance) {converged=true;status="converged";break;}
            // Do not oversolve a restricted subproblem to floating-point
            // stagnation while omitted knots have much larger violations.
            // This changes admission timing, never the final full-knot gate.
            if(diag.active<=std::max(Real(tolerance),Real(.001)*diag.stationarity) && diag.candidate>=0 && lambda>0 && std::isfinite(lambda)) {
                int at=std::lower_bound(active.begin(),active.end(),diag.candidate)-active.begin();
                signs.insert(signs.begin()+at,fullg[diag.candidate+1]-geometry.center[diag.candidate]*fullg[0]>0 ? -1 : 1);
                active.insert(active.begin()+at,diag.candidate);++additions;rebuild=true;continue;
            }
            Vec penalty=Vec::Zero(B.d);
            if(std::isfinite(lambda))for(int j=0;j<int(active.size());++j)penalty[j+3]=lambda*signs[j];
            Vec gb=B.transpose(e.curvegradient)+B.transpose_transform(penalty),dc,dg;
            if(!direction(B,e,Z,w,gb,zg,nuisance_ridge/mass,curve_ridge/mass,dc,dg,damping)) {
                if(diag.candidate>=0 && diag.active<100*tolerance && lambda>0 && std::isfinite(lambda)) {
                    int at=std::lower_bound(active.begin(),active.end(),diag.candidate)-active.begin();
                    signs.insert(signs.begin()+at,fullg[diag.candidate+1]-geometry.center[diag.candidate]*fullg[0]>0 ? -1 : 1);active.insert(active.begin()+at,diag.candidate);++additions;rebuild=true;continue;
                }
                status="singular_banded_system";break;
            }
            Vec dt=B.coefficients(dc);Real step=1;int drop=-1;
            if(lambda>0 && std::isfinite(lambda))for(int j=0;j<int(active.size());++j) {
                Real a=theta[j+3]*signs[j],b=dt[j+3]*signs[j];
                if(b<0 && -std::max(a,Real(0))/b<step) {step=-std::max(a,Real(0))/b;drop=j;}
            }
            if(step<1e-12 && drop>=0) {
                active.erase(active.begin()+drop);signs.erase(signs.begin()+drop);++removals;rebuild=true;continue;
            }
            Real initialstep=step,slope=gb.dot(dc)+zg.dot(dg);bool accepted=false;
            lastslope=slope;
            for(int ls=0;ls<60;++ls) {
                Vec nextc=c+step*dc,nextgamma=gamma+step*dg;
                Eval next=evaluate(B.values(nextc),Z,nextgamma,off,y,w,nuisance_ridge/mass,curve_ridge/mass,B.coefficients(nextc),lambda,gaussian);
                Real rounding=16*std::numeric_limits<Real>::epsilon()*std::max(Real(1),std::abs(e.objective));
                if(std::isfinite(next.objective) && next.objective<=e.objective+1e-4*step*slope+rounding) {
                    c=nextc;gamma=nextgamma;curve=B.values(c);accepted=true;break;
                }
                step*=.5;++backtracks;
            }
            laststep=step;
            if(!accepted) {status="line_search_failed";break;}
            if(drop>=0 && step==initialstep) {active.erase(active.begin()+drop);signs.erase(signs.begin()+drop);++removals;rebuild=true;}
        }
        if(rebuild) {B=Basis(kk,active,xx);c=B.project(curve,w);}
        theta=B.coefficients(c);curve=B.values(c);
        e=evaluate(curve,Z,gamma,off,y,w,nuisance_ridge/mass,curve_ridge/mass,theta,lambda,gaussian);
        Vec fullg=geometry.gradient(e.curvegradient),zg=Z.transpose()*e.residual+(nuisance_ridge/mass)*gamma;
        diag=diagnose(geometry,fullg,zg,theta,active,signs,lambda);
        converged=converged && e.eta.allFinite() && std::isfinite(e.objective) && diag.stationarity<=tolerance && diag.complementarity<=gap_tolerance;
        if(!converged && status=="converged")status="diagnostic_check_failed";
        Vec coeff=Vec::Zero(p),dual=Vec::Zero(p);coeff[0]=theta[0];
        for(int j=0;j<q;++j)coeff[j+1]=gamma[j];
        for(int j=0;j<2;++j){coeff[q+1+j]=theta[j+1]*geometry.scale[j];coeff[0]+=theta[j+1]*geometry.center[j];}
        for(int j=0;j<int(active.size());++j){int at=active[j];coeff[q+1+at]=theta[j+3]*geometry.scale[at];coeff[0]+=theta[j+3]*geometry.center[at];}
        for(int j=2;j<m-1;++j)dual[q+1+j]=diag.dual[j+1];
        Real dualfeasibility=std::max(Real(0),infnorm(dual)-1);
        converged=converged && dual.allFinite() && dualfeasibility==0;
        if(!converged && status=="converged")status="diagnostic_check_failed";
        Rcpp::IntegerVector aa=Rcpp::wrap(active);
        Rcpp::NumericVector witness(p),curvewitness(n);
        if(!Rcpp::as<bool>(init["converged"])) {
            Rcpp::NumericVector initialwitness=init["separation.direction"];
            std::copy(initialwitness.begin(),initialwitness.end(),witness.begin());
        }
        if(separated_knot>=0) for(int i=0;i<n;++i) if(xx[i]==kk[separated_knot]) curvewitness[i]=y[i]==0 ? -1 : 1;
        fits[li]=Rcpp::List::create(Rcpp::_["coefficients"]=numeric(coeff),Rcpp::_["dual"]=numeric(dual),
            Rcpp::_["linear.predictors"]=numeric(e.eta),Rcpp::_["fitted.values"]=numeric(e.probability),
            Rcpp::_["separation.direction"]=witness,Rcpp::_["separation.curve.direction"]=curvewitness,
            Rcpp::_["objective"]=double(e.objective*mass),
            Rcpp::_["converged"]=converged,Rcpp::_["status"]=status,Rcpp::_["stationarity"]=double(diag.stationarity),
            Rcpp::_["complementarity"]=double(diag.complementarity),Rcpp::_["dual.feasibility"]=double(dualfeasibility),
            Rcpp::_["iterations"]=iterations,Rcpp::_["barrier.stages"]=0,Rcpp::_["barrier.parameter"]=0.,
            Rcpp::_["backtracks"]=backtracks,Rcpp::_["max.damping"]=double(damping),Rcpp::_["max.abs.eta"]=double(infnorm(e.eta)),
            Rcpp::_["active.knots"]=aa,Rcpp::_["compact.coefficients"]=numeric(c),Rcpp::_["nuisance.coefficients"]=numeric(gamma),
            Rcpp::_["working.set.additions"]=additions,Rcpp::_["working.set.removals"]=removals,
            Rcpp::_["initialization.iterations"]=init["iterations"],
            Rcpp::_["last.line.search.step"]=double(laststep),Rcpp::_["last.directional.derivative"]=double(lastslope));
    }
    return fits;
}

// [[Rcpp::export]]
Rcpp::List logistic_tf_banded_native(Rcpp::NumericVector u, Rcpp::NumericVector knots,
    Rcpp::NumericMatrix nuisance, Rcpp::NumericVector response,
    Rcpp::NumericVector weights, Rcpp::NumericVector offsets, Rcpp::NumericVector lambdas,
    double nuisance_ridge, double curve_ridge, double tolerance,
    double gap_tolerance, int max_iterations) {
    return univariate_tf_banded_native(u, knots, nuisance, response, weights, offsets,
        lambdas, nuisance_ridge, curve_ridge, tolerance, gap_tolerance, max_iterations, false);
}

// [[Rcpp::export]]
Rcpp::List gaussian_tf_banded_native(Rcpp::NumericVector u, Rcpp::NumericVector knots,
    Rcpp::NumericMatrix nuisance, Rcpp::NumericVector response,
    Rcpp::NumericVector weights, Rcpp::NumericVector offsets, Rcpp::NumericVector lambdas,
    double nuisance_ridge, double curve_ridge, double tolerance,
    double gap_tolerance, int max_iterations) {
    return univariate_tf_banded_native(u, knots, nuisance, response, weights, offsets,
        lambdas, nuisance_ridge, curve_ridge, tolerance, gap_tolerance, max_iterations, true);
}
