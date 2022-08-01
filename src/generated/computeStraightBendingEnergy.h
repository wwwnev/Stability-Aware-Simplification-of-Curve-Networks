#ifndef CODEGEN_COMPUTESTRAIGHTBENDINGENERGY_H
#define CODEGEN_COMPUTESTRAIGHTBENDINGENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeStraightBendingEnergy(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, double& energy);
void computeStraightBendingEnergyGradient(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, Eigen::Matrix<double, 9, 1>& gradient);
void computeStraightBendingEnergyHessian(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, Eigen::Matrix<double, 9, 9>& hessian);
void computeStraightBendingEnergyHessianProductJacobian(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, const Eigen::Matrix<double,9,1> & v, Eigen::Matrix<double, 9, 9>& hessianProductJacobian);
void computeStraightBendingEnergyRestConfJac(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, Eigen::Matrix<double, 9, 0>& rest_conf_jac);
void computeStraightBendingEnergyParam_jac(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, Eigen::Matrix<double, 9, 3>& param_jac);
void computeStraightBendingEnergyParam_thirdOrderDerivatives(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, Eigen::Matrix<double, 9, 9>& d3_0, Eigen::Matrix<double, 9, 9>& d3_1, Eigen::Matrix<double, 9, 9>& d3_2);
void computeStraightBendingEnergyThirdOrderDerivatives(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, Eigen::Matrix<double, 9, 9>& d3_0, Eigen::Matrix<double, 9, 9>& d3_1, Eigen::Matrix<double, 9, 9>& d3_2, Eigen::Matrix<double, 9, 9>& d3_3, 
	Eigen::Matrix<double, 9, 9>& d3_4, Eigen::Matrix<double, 9, 9>& d3_5, Eigen::Matrix<double, 9, 9>& d3_6, Eigen::Matrix<double, 9, 9>& d3_7, Eigen::Matrix<double, 9, 9>& d3_8);
void computeStraightBendingEnergyRestConfThirdOrderDerivatives(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1);
void computeStraightBendingEnergyFourthOrderDerivatives(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, const Eigen::Matrix<double,9,1> & v1, const Eigen::Matrix<double,9,1> & v2, Eigen::Matrix<double, 9, 9>& fourthOrderDerivative);
void computeStraightBendingEnergyThirdOrderToRestConfFourthOrderDerivatives(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, const Eigen::Matrix<double,0,1> & v1, const Eigen::Matrix<double,9,1> & v2, Eigen::Matrix<double, 9, 9>& thirdOrderToRestConfFourthOrderDerivative);
void computeStraightBendingEnergyThirdOrderToParameterFourthOrderDerivatives(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double norm0, 
	double norm1, const Eigen::Matrix<double,3,1> & v1, const Eigen::Matrix<double,9,1> & v2, Eigen::Matrix<double, 9, 9>& thirdOrderToParameterFourthOrderDerivative);

 } 
#endif