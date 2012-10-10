/*
 * File:   HmcSampler.h
 * Author: aripakman
 *
 * Created on July 4, 2012, 10:44 AM
 */

#ifndef HMCSAMPLER_H
#define HMCSAMPLER_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <tr1/random>
#include <vector>
#include <Eigen/Dense>

//using namespace Eigen;
using namespace std;
using namespace std::tr1;

struct LinearConstraint{
  Eigen::VectorXd f;
  double  g;
};

struct QuadraticConstraint{
  Eigen::MatrixXd A;
  Eigen::VectorXd B;
  double  C;
};


class HmcSampler   {
public:

  HmcSampler();
  HmcSampler(const int & d, const int & seed);

  void setInitialValue(const Eigen::VectorXd & initial);
  void addLinearConstraint(const Eigen::VectorXd & f, const double & g);
  void addQuadraticConstraint(const Eigen::MatrixXd & A, const Eigen::VectorXd & B, const double & C);
  Eigen::MatrixXd sampleNext(bool returnTrace = false);

  static Eigen::MatrixXd rtnorm(const Eigen::MatrixXd& b,
				const Eigen::MatrixXd& P,
				const Eigen::MatrixXd& F,
				const Eigen::VectorXd& g,
				const Eigen::VectorXd& iv,
				int samp=1,
				int burn=0,
				bool kern=false,
				int seed=-1);

  double verifyConstraints(const Eigen::VectorXd &v);

private:
  int dim;
  Eigen::VectorXd lastSample;
  static const double min_t;
  vector<LinearConstraint> linearConstraints;
  vector<QuadraticConstraint> quadraticConstraints;

  ranlux64_base_01 eng1;
  //    mt19937 eng1; //to sample time and momenta
  uniform_real<> ud;
  normal_distribution<> nd;

  void _getNextLinearHitTime(const Eigen::VectorXd & a, const Eigen::VectorXd & b,  double & t, int & cn );
  void _getNextQuadraticHitTime(const Eigen::VectorXd & a, const Eigen::VectorXd & b, 
				double & t, int & cn, const bool );
  double _verifyConstraints(const Eigen::VectorXd &);
  void _updateTrace( Eigen::VectorXd const & a,  Eigen::VectorXd const & b, double const & tt, Eigen::MatrixXd & tracePoints);
};

#endif  /* HMCSAMPLER_H */

