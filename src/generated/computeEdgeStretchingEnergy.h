#ifndef CODEGEN_COMPUTEEDGESTRETCHINGENERGY_H
#define CODEGEN_COMPUTEEDGESTRETCHINGENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeEdgeStretchingEnergy(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, double& energy);
void computeEdgeStretchingEnergyGradient(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, Eigen::Matrix<double, 6, 1>& gradient);
void computeEdgeStretchingEnergyHessian(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, Eigen::Matrix<double, 6, 6>& hessian);
void computeEdgeStretchingEnergyHessianProductJacobian(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, const Eigen::Matrix<double,6,1> & v, 
	Eigen::Matrix<double, 6, 6>& hessianProductJacobian);
void computeEdgeStretchingEnergyRestConfJac(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, Eigen::Matrix<double, 6, 0>& rest_conf_jac);
void computeEdgeStretchingEnergyParam_jac(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, Eigen::Matrix<double, 6, 2>& param_jac);
void computeEdgeStretchingEnergyParam_thirdOrderDerivatives(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, Eigen::Matrix<double, 6, 6>& d3_0, 
	Eigen::Matrix<double, 6, 6>& d3_1);
void computeEdgeStretchingEnergyThirdOrderDerivatives(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, Eigen::Matrix<double, 6, 6>& d3_0, 
	Eigen::Matrix<double, 6, 6>& d3_1, Eigen::Matrix<double, 6, 6>& d3_2, Eigen::Matrix<double, 6, 6>& d3_3, Eigen::Matrix<double, 6, 6>& d3_4, Eigen::Matrix<double, 6, 6>& d3_5);
void computeEdgeStretchingEnergyRestConfThirdOrderDerivatives(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm);
void computeEdgeStretchingEnergyFourthOrderDerivatives(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, const Eigen::Matrix<double,6,1> & v1, 
	const Eigen::Matrix<double,6,1> & v2, Eigen::Matrix<double, 6, 6>& fourthOrderDerivative);
void computeEdgeStretchingEnergyThirdOrderToRestConfFourthOrderDerivatives(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, const Eigen::Matrix<double,0,1> & v1, 
	const Eigen::Matrix<double,6,1> & v2, Eigen::Matrix<double, 6, 6>& thirdOrderToRestConfFourthOrderDerivative);
void computeEdgeStretchingEnergyThirdOrderToParameterFourthOrderDerivatives(double stretchingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double init_e_norm, const Eigen::Matrix<double,2,1> & v1, 
	const Eigen::Matrix<double,6,1> & v2, Eigen::Matrix<double, 6, 6>& thirdOrderToParameterFourthOrderDerivative);

 } 
#endif