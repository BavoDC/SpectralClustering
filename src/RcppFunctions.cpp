// [[Rcpp::depends(RcppEigen)]]
#define ARMA_64BIT_WORD 1
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXf;

// [[Rcpp::export]]
Eigen::VectorXd XtV(Eigen::SparseMatrix<double> X, const Eigen::Map<Eigen::VectorXd> V){
  Eigen::VectorXd C = X.transpose() * V;
  return C;
}

// [[Rcpp::export]]
Eigen::VectorXd XB(Eigen::SparseMatrix<double> X, const Eigen::Map<Eigen::VectorXd> B){
  Eigen::VectorXd C = X * B;
  return C;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> XXt(Eigen::SparseMatrix<double> X, Eigen::SparseMatrix<double> Xt){
  Eigen::SparseMatrix<double> C = X * Xt;
  return C;
}

// [[Rcpp::export]]
Eigen::MatrixXd GSM(const Eigen::Map<Eigen::MatrixXd>& A, const double sigma = 1) {
  int n = A.rows();
  //Eigen::MatrixXd S(n, n);
  Eigen::MatrixXd S = Eigen::VectorXd::Zero(n).asDiagonal();
  for(int i = 0; i < n; i++)
    for(int j = i + 1; j < n; j++)
      S(j, i) = S(i, j) = std::exp(1 / (2 * pow(sigma, 2)) *  - pow((A.row(i) - A.row(j)).norm(), 2));
  return S;
}

// [[Rcpp::export]]
Eigen::MatrixXd GSMLocalScaling(const Eigen::Map<Eigen::MatrixXd>& A, const Eigen::Map<Eigen::VectorXd>& sigma) {
  int n = A.rows();
  //Eigen::MatrixXd S(n, n);
  Eigen::MatrixXd S = Eigen::VectorXd::Zero(n).asDiagonal();
  for(int i = 0; i < n; i++)
    for(int j = i + 1; j < n; j++)
      S(j, i) = S(i, j) = std::exp(1 / (sigma[i] * sigma[j]) *  - pow((A.row(i) - A.row(j)).norm(), 2));
  return S;
}

// [[Rcpp::export]]
double Norm(const Eigen::Map<Eigen::VectorXd>& v) {
  double vn = v.norm();
  return vn;
}

// [[Rcpp::export]]
void redir(const std::string file){
  FILE* F = freopen(file.c_str(),"w+",stdout);
  }

// [[Rcpp::export]]
void resetredir(){
  FILE* F = freopen("CON","w+",stdout);
  }

// [[Rcpp::export]]
Eigen::MatrixXd CosineSimilarity(const Eigen::Map<Eigen::MatrixXd>& X) {
  int n = X.rows();
  //Eigen::MatrixXd S(n, n);
  Eigen::MatrixXd S = Eigen::VectorXd::Ones(n).asDiagonal();
  double num, denom;
  for(int i = 0; i < n; i++) {
    for(int j = i + 1; j < n; j++) {
      num   = (X.row(i) * X.row(j).transpose());
      denom = (X.row(i).norm() * X.row(j).norm());
      S(j, i) = S(i, j) = num / denom;
    }
  }
  return S;
}

// [[Rcpp::export]]
Eigen::MatrixXd Qm(const Eigen::Map<Eigen::MatrixXd>& A) {
  int n = A.rows();
  //Eigen::MatrixXd S(n, n);
  Eigen::MatrixXd Q = Eigen::VectorXd::Zero(n).asDiagonal();
  for(int i = 0; i < n; i++)
    for(int j = i + 1; j < n; j++)
      Q(j, i) = Q(i, j) = pow(1 + pow((A.row(i) - A.row(j)).norm(), 2), -1);
  return Q;
}