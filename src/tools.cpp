#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() : mse(4), 
    lastmse(4),
    rmse(4)
{
  resetRMSE();
}

Tools::~Tools() {}

/**
 * Calculate the RMSE
 * 
 * Logic:
 * 1. check input
 *    - check the size of the estimations and the ground truth,
 *      make sure they are the same. 
 * 2. cal the sum of all previous terms 
 * 3. cal the incremental term 
 * 4. add the incremental term to sum 
 * 5. cal the new rmse
 * 6. warn when the rmse larger than the limit
 * 
 **/
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{

    int number_estimations = estimations.size();
    int number_ground_truths = ground_truth.size();

    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    if (number_estimations < 1 || number_estimations != number_ground_truths) {
        std::cout << "Invalid length of input" << std::endl;
        return rmse;
    }

    for (int i = 0; i < number_estimations; ++i) {
        VectorXd residual_error = estimations[i] - ground_truth[i];
        residual_error = residual_error.array() * residual_error.array();
        rmse += residual_error;
    }

    return (rmse / number_estimations).array().sqrt();
}
