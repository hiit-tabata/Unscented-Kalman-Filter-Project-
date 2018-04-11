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

void Tools::resetRMSE()
{
  mse.fill(0.);
  lastmse.fill(0.);
}

// Calculate the RMSE
/**
 * The function is calculating the RMSE.
 * 
 * Logic:
 * 1. check the size of the estimations and the ground truth,
 * make sure they are the same. 
 * 2. cal the last term and add cal the rmse
 * 
 **/
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
  // get idx
  float idx = estimations.size();

  // validate data
  if(idx == 0 ){
    cout << "estimations should not be 0 " << endl;
  } else if(idx != ground_truth.size()){
    cout << "estimations size not = to ground truth." << endl;
  }

  VectorXd tmp = mse * (idx-1);
  //get the last term error
  VectorXd nextTerm = estimations[idx-1] - ground_truth[idx-1]; 
  nextTerm = nextTerm.array() * nextTerm.array();
  mse = (lastmse + nextTerm)/estimations.size();
  lastmse = mse;

  // Calculate the RMSE
  rmse = mse.array().sqrt();

  if( rmse(0) > .09 ||
      rmse(1) > .10 ||
      rmse(2) > .40 ||
      rmse(3) > .30 )
    {
      cout << "Warning at timestep " << idx << ":  rmse = " 
          << rmse(0) << "  " << rmse(1) << "  " 
          << rmse(2) << "  " << rmse(3) << endl
          << " It learger than  .09, .10, .40, .30" << endl;
    }

  return rmse;
}
