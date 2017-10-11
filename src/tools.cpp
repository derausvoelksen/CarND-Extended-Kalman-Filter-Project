#include <iostream>
#include "tools.h"
#define eps 1e-4
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  //Taken from lesson
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  //checking for invalid input dimension
  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
  	cout << "Invalid estimation or ground_truth data" << endl;
  	return rmse;
  }
  
  //sum up squared residuals
  for(unsigned int i = 0; i < estimations.size(); ++i){
  	VectorXd residual = estimations[i] - ground_truth[i];
  	residual = residual.array() * residual.array();
  	rmse += residual;
  }
  
  // calculate mean
  rmse = rmse / estimations.size();
  
  // calculate square root
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);

    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    double ps = sqrt(px * px + py * py);
    double p2 = px * px + py * py;
    double p3 = ps * ps * ps;

    if (ps < eps) {
        return Hj;
    }
    //compute the Jacobian matrix
    Hj << px / ps, py / ps, 0, 0,
            -py / (p2), px / (p2), 0, 0,
            (py * (vx * py - vy * px)) / p3, (px * (vy * px - vx * py)) / p3, px / ps, py / ps;

    return Hj;
}

VectorXd Tools::FromCartesian(const VectorXd &x) {
    cout <<"X: " <<x <<endl;

    VectorXd z = VectorXd(3);
    double rho = sqrt(x[0] * x[0] + x[1] * x[1]);
    z << rho, atan2(x[1], x[0]), (x[0] * x[2] + x[1] * x[3]) / rho;
    cout <<"Z: " <<z <<endl;

    return z;
}