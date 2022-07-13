/**
   This file provides a mechanism to perform QP based Inverse Dynamics with
constraints.
   @author Victor Paredes (paredescauna.1@osu.edu)
*/

#ifndef IDQPcpp_HPP_
#define IDQPcpp_HPP_

#include "include/eiquadprog/eiquadprog-fast.hpp"
#include "Eigen/Dense"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace IDQP_type {
  // Variable Helpers to load many functions
  typedef std::vector< std::function< MatrixXd(VectorXd&)> > List_Jac;
  typedef std::vector< std::function< MatrixXd(VectorXd&, VectorXd&)> > List_Jac_dot;
  typedef std::vector< std::function< MatrixXd(VectorXd&)> > List_Ad;
  typedef std::vector< std::function< VectorXd(VectorXd&)> > List_Position;
  } // namespace IDQP_type

using namespace IDQP_type;

class IDQPcpp{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /**
   * Constructor to initialize controller.
   *  @param nStates Number of joints of the Robot, including base joints, i.e ne
   *  + nb
   *  @param nObj Number of Objectives of the virtual constraint
   *  @param nActuators Number of Actuated joints
   *  @param nContacts Number of contacts, we assume 6 constraints per contact.
   **/
  IDQPcpp(int nStates, int nObj, int nActuators, int nContacts);

  // ---------------- Dynamics ---------------- //
  // This a required function to have a well posed IDQP //
  // We call only one of these functions to specify the dynamics //

  /**
   * In this case Cmat = Coriolis.
   * @param Mmat Function pointer to the Mass Matrix that takes q
   * @param Cmat Function pointer to the Coriolis Matrix that takes (q,dq)
   * @param Gvec Function pointer to the Gravity Vector that takes q
   **/
  int addDynamicsCMat(std::function< MatrixXd (VectorXd&)> Mmat,
                  std::function< MatrixXd (VectorXd&, VectorXd&)> Cmat,
                  std::function< VectorXd (VectorXd&)> Gvec);

  /**
   * In this case Cvec = Coriolis * dq
   * @param Mmat Function pointer to the Mass Matrix that takes q
   * @param Cvec Function pointer to the Coriolis Vector that takes (q,dq)
   * @param Gvec Function pointer to the Gravity Vector that takes q
   **/  
  int addDynamicsCVec(std::function< MatrixXd (VectorXd&)> Mmat,
                  std::function< VectorXd (VectorXd&, VectorXd&)> Cvec,
                  std::function< VectorXd (VectorXd&)> Gvec);

  /**
   * In this case Hvec = Coriolis * dq + Gravity
   * @param Mmat Function pointer to the Mass Matrix that takes q
   * @param Hvec Function pointer to the Drift Vector thatf takes (q,dq)
   **/  
  int addDynamicsHvec(std::function< MatrixXd (VectorXd&)> Mmat,
                  std::function< VectorXd (VectorXd&, VectorXd&)> Hvec);

  int addDynamicsHvecNum(MatrixXd Mmat, VectorXd Hvec);  
  int updateDynamicsHvecNum(MatrixXd Mmat, VectorXd Hvec);
  // ---------------- End Dynamics ---------------- //  

  /**
   * Required. Specify a list of contact jacobians and its derivatives, to declare the list use the internal data types: List_Jac and List_Jac_dot that will store function pointers.
   * @param Jac List of function pointers to the contact jacobians with argument q
   * @param Jac_dot List of function pointers to the time-derivative of contact jacobians with argument (q,dq)
   **/    
  int addContactJacobians(List_Jac Jac,
                          List_Jac_dot Jac_dot);

  /**
   * Required. Specify mu (friction coefficient) in the friction cone.
   * @param m_friction Friction coefficient
   **/      
  int addFrictionCone(double mu_friction);

  /**
   * Optional. Specify the control limits on the controller.
   * @param uMin Vector of lower bounds for the control signal u
   * @param uMax Vector of upper bounds for the control signal u
   **/        
  int addControlLimits(VectorXd &uMin, VectorXd &uMax);

  int addPositionLimits(VectorXd &qMin, VectorXd &qMax, VectorXd &IndexMap);

  /**
   * Optional. Add the ZMP constraint indicating the positions of the contacts with the
ground and the Adjoint Matrices to project the wrenches to the World
Coordinates.
* All the following parameters are relative to the average foot position
   * @param LengthXmin Minimum X location of a feasible ZMP
   * @param LengthXmax Maximum X location of a feasible ZMP
   * @param LengthYmin Minimum Y location of a feasible ZMP
   * @param LengthYmax Maximum Y location of a feasible ZMP
   **/          
  int addZMP(double LengthXmin, double LengthXmax, double LengthYmin, double LengthYmax, List_Ad Ad_foot, List_Position Pos_foot);

  int addCustomEq(MatrixXd &Acustom, VectorXd &Bcustom);
  int addCustomIneq(MatrixXd &Acustom, VectorXd &Bcustom);

  int addTypeIconstraint(List_Jac M_I, int nConst_I);
  int addTypeIIconstraint(List_Jac M_II, List_Jac_dot M_II_dot, int nConst_II);    

  // ---------------- Virtual Constraints ---------------- //
  // Call only one of these functions

  /**
   * Required. Specify the actual and desired virtual constraints together with the
desired Kp and Kd gains for the those ouputs.
   * We assume that the desired output acceleration is Zero.
   * All the following parameters with except Kp and Kd are function pointers.
   * @param y Actual Output, depends on (q)
   * @param Jac_y Jacobian of the Output, depends on (q)
   * @param Jac_y_dot Time-derivative of the output jacobian, depends on (q,dq)
   * @param y_d Desired Output, depends on (t, q, dq, qdes, dqdes)
   * @param Dy_d Time-derivative of the desired Output, depends on (t, q, dq, qdes, dqdes)
   * @param Kp Vector of Proportional Gains of the Outputs
   * @param Kd Vector of Derivative Gains of the Outputs
   **/            
  int addObjectives(std::function< VectorXd (VectorXd&)> y,
                    std::function< MatrixXd (VectorXd&)> Jac_y,
                    std::function< MatrixXd (VectorXd&, VectorXd&)> Jac_y_dot,
                    std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> y_d,
                    std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> Dy_d,
                    VectorXd &Kp, VectorXd &Kd );

  /**
   * Required. Specify the actual and desired virtual constraints together with the
desired Kp and Kd gains for the those ouputs.
   * We DON'T assume that the desired output acceleration is Zero.
   * All the following parameters with except Kp and Kd are function pointers.
   * @param y Actual Output, depends on (q)
   * @param Jac_y Jacobian of the Output, depends on (q)
   * @param Jac_y_dot Time-derivative of the output jacobian, depends on (q,dq)
   * @param y_d Desired Output, depends on (t, q, dq, qdes, dqdes)
   * @param Dy_d Time-derivative of the desired Output, depends on (t, q, dq,
qdes, dqdes)
   * @param DDy_d Second Time-derivative of the desired Output, depends on (t, q, dq, qdes, dqdes)
   * @param Kp Vector of Proportional Gains of the Outputs
   * @param Kd Vector of Derivative Gains of the Outputs
   **/              
  int addObjectives(std::function< VectorXd (VectorXd&)> y,
                    std::function< MatrixXd (VectorXd&)> Jac_y,
                    std::function< MatrixXd (VectorXd&, VectorXd&)> Jac_y_dot,
                    std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> y_d,
                    std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> Dy_d,
                    std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> DDy_d,                    
                    VectorXd &Kp, VectorXd &Kd );
  // ---------------- End Virtual Constraints ---------------- //

  /**
   * Optional. (Default 1e-6) Specify the regularization terms for the joint accelerations, controls and
   wrenches.
   * For now, we set one term per class of variable: [joint accelerations,
   controls, wrenches]
   * @param q_ddot Specify the regularization term for joint accelerations
   * @param u_control Specify the regularization term for the control signals
   * @param lambda Specify the regularization term for the Wrenches
   **/                
  int setRegularizationTerm(double q_ddot, double u_control, double lambda);

  /**
   * Optional. Specify the joint acceleration bounds
   * @param qddMin Specify the Vector of lower bounds for the joint
accelerations.
   * @param qddMax Specify the Vector of upper bounds for the joint accelerations.
   **/                  
  int setAccelerationLimits(VectorXd &qddMin, VectorXd &qddMax);

  int specifyMotorTransmissions(VectorXd &transmissions);

  /**
   * Required. Builds the optimization problem. This must be called as the last command of the IDQP members. 
   **/                    
  int buildOptimization();
  //int buildCLFQP

  /**
   * Required. Computes the decision variable of the Optimization Problem.
   * The decision variable is [Qddot, u, Wrenches] with dimensions [nStates_,
   nActuators_, nContactcstr_] respectively.
   * Since the user constructs the desired outputs, we give them all the
arguments to have freedom on their construction.
   * @param t Actual time.
   accelerations.
   * @param q Actual joint coordinates
   * @param dq Actual joint velocities
   * @param qdes Desired joint coordinates
   * @param dqdes Desired joint velocities
   * @return The decision Variable X
   **/                    
  VectorXd solve(double t, VectorXd &q, VectorXd &dq, VectorXd &qdes, VectorXd &dqdes);

  int testEvaluation(VectorXd &q, VectorXd &dq);

  /**
   * Utility (Untested). Computes the current estimated ZMP
   * @param q Actual joint coordinates
   * @return The ZMP prediction as a VectorXd
   **/                      
  VectorXd getZMPprediction(VectorXd &q);

  MatrixXd solveRiccatiArimotoPotter();

  /**
   * Utility (Untested). Computes the current computed Vector of joint accelerations
   * @return The qdd solutions as a VectorXd
   **/                        
  VectorXd get_qdd();

  VectorXd get_torque();

private:
  int nContacts_;
  
  int nStates_;
  int nObj_;
  int nActuators_;
  int nContactCstr_;
  int nZMPCstr_;
public:
  int dimX;
private:

  double mu_friction;
  // Dynamics Matrix
  MatrixXd ADyn;
  VectorXd BDyn;

  // Contact Contraint
  MatrixXd ACont;
  VectorXd BCont;

  // Friction Cone Constraint
  MatrixXd ACone;
  VectorXd BCone;

  // ZMP Constraints
public:
  MatrixXd AZMP;
  MatrixXd AZMPr;
  VectorXd BZMP;
private:
  List_Ad List_Ad_contacts;
  List_Position List_Pos_contacts;

  // Custom Constraints
  MatrixXd Acustom_eq;
  VectorXd Bcustom_eq;
  MatrixXd Acustom_ineq;
  VectorXd Bcustom_ineq;

  // Custom TypeI & TypeII
  List_Jac List_MI;
  List_Jac List_MII;
  List_Jac_dot List_MII_dot;
  
  MatrixXd M_I_mat;
  VectorXd b_I_vec;
  MatrixXd M_II_mat;
  MatrixXd M_II_dot_mat;
  VectorXd b_II_vec;
  int nConstI;
  int nConstII;  

  // Torque Constraint
  MatrixXd Atorque;
  VectorXd Btorque;

  // Position Constraint
  MatrixXd Aposition;
  VectorXd Bposition;
  VectorXd PositionIndexMap;
  VectorXd qmin;
  VectorXd qmax;

  // Acceleration Constraint
  MatrixXd Aacc;
  VectorXd Bacc;

  // Pointer to functions
  // Dynamics
  std::function< MatrixXd(VectorXd&)> Mmat;
  std::function< MatrixXd(VectorXd&, VectorXd&)> Cmat;
  std::function< VectorXd(VectorXd&, VectorXd&)> Cvec;
  std::function< VectorXd(VectorXd&, VectorXd&)> Hvec;  
  std::function< VectorXd(VectorXd&)> Gvec;

  // Numeric Dynamics
  MatrixXd MmatNum;
  VectorXd HvecNum;

  // Contacts
  List_Jac List_Jac_cont;
  List_Jac_dot List_Jac_cont_dot;
  MatrixXd Jac_cont;
  MatrixXd Jac_cont_dot;
  MatrixXd Bmat;

  // Objectives
  std::function< VectorXd(VectorXd&)> ya;
  std::function< VectorXd(double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> yd;
  std::function< VectorXd(double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> Dyd;
  std::function< VectorXd(double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> DDyd; 
  std::function< MatrixXd(VectorXd&)> Jac_y;
  std::function< MatrixXd(VectorXd&, VectorXd&)> Jac_y_dot;
  VectorXd yStar;
  MatrixXd MatKp;
  MatrixXd MatKd;

public:
  VectorXd yOut;
  VectorXd DyOut;
  VectorXd ydOut;
  VectorXd DydOut;
  VectorXd yaOut;
  VectorXd DyaOut;
  VectorXd DDydOut;  
  MatrixXd Jacy;
  MatrixXd DJacy;  
private:
  // Transmissions
  VectorXd MotorjointTrans;

  // Control Barrier Functions (CBF)
  // Joint Limits
  double h_min;
  double h_max;
  double alpha_cbf;
  double k_cbf;
  VectorXd Breact;

  // Control Lyapunov Function (CLF)
  MatrixXd Flin;
  MatrixXd Glin;
  MatrixXd Pmat;
  MatrixXd Pemat;
  MatrixXd Qmat;
  MatrixXd Rmat;
  VectorXd eta;
  MatrixXd Aclf;
  VectorXd Bclf;
  MatrixXd Ie;
  double c_clf;
  double eps_clf;
  Eigen::EigenSolver<Eigen::MatrixXd> *Eigs;
  
  // Optimization
  // Aeq X = Beq
public:
  MatrixXd Aeq;
  VectorXd Beq;  

  int nEq;
  int EqCounter;
  double dt;

  // Aineq X \leq Bineq
  MatrixXd Aineq;
  VectorXd Bineq;
  int nIneq;
  int IneqCounter;
  int dyn_IneqCounter;

  // Cost Function X^T H X + f^T X

  MatrixXd Hopt;
  VectorXd fopt;

  VectorXd ZeroVec;
  int solutionStatus; // 0:Not solved, 1:solved  
private:
  // Regularization
  MatrixXd Reg;

  struct{
    bool dynamics;
    bool contacts;
    bool friction;
    bool zmp;
    bool torque;
    bool acceleration;
    bool velocity;
    bool position;
    bool customEq;
    bool customIneq;
    bool typeIEq;
    bool typeIIEq;    
  } activeConstraints;

  eiquadprog::solvers::EiquadprogFast qp;
  VectorXd Xsol;
  VectorXd x0;
};

#endif
