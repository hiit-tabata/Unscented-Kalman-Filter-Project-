#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

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
