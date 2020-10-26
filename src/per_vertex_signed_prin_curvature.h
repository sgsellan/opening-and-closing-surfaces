#ifndef PC_H
#define PC_H
#include "per_vertex_signed_prin_curvature.cpp"


#include <Eigen/Core>

void per_vertex_signed_prin_curvature(Eigen::MatrixXd & V,Eigen::MatrixXi & F, Eigen::MatrixXd & PD1, Eigen::MatrixXd & PD2, Eigen::VectorXd & PC1, Eigen::VectorXd & PC2);


#endif
