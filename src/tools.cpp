#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * RMSE CALCULATION....!!
  */

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	if(estimations.size() == 0 || estimations.size() != ground_truth.size())
	{
		cout << "Invalid estimations or ground ground_truth values" << endl;
		return rmse;
	}

	for (int i =0; i < estimations.size(); i++)
	{
		VectorXd residual_error;
		residual_error = (estimations[i] - ground_truth[i]);

		residual_error = residual_error.array()*residual_error.array();

		rmse += residual_error;
	}

	rmse = rmse/estimations.size();

	rmse = rmse.array().sqrt();

	return rmse;
}