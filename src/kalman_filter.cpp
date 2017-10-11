#include "kalman_filter.h"
#include "tools.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd y = z - H_ * x_;
    UpdateInternal(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    H_ = Tools::CalculateJacobian(x_);
    VectorXd y = z - Tools::FromCartesian(x_);
    y[1] = atan2(sin(y[1]), cos(y[1]));
    UpdateInternal(y);
}

void KalmanFilter::UpdateInternal(const VectorXd &y) {
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}