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
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
    rmse << 0,0,0,0;

    if((estimations.size()==0) | (ground_truth.size()!=estimations.size())){
	        return rmse;
    }

    VectorXd squared_res(4);
    squared_res << 0,0,0,0;

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
        VectorXd diff=estimations[i]-ground_truth[i];
        VectorXd fact=diff.array().square();
		 squared_res=squared_res+fact;
	}

	//calculate the mean
    VectorXd mean = squared_res.array()/estimations.size();
	//calculate the squared root
    rmse=mean.array().sqrt();
    //cout << "RSME = " << rmse << endl;
	//return the result
	return rmse;
}
