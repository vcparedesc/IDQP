#ifndef MOTION_GENERATOR_H_
#define MOTION_GENERATOR_H_

#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include "math_utils.hpp"
#include "robot_expressions.hpp"
#include "math_utils.hpp"
#include "Delta_Rule.hpp"
#include <memory>

using namespace Eigen;
typedef Matrix<double, 6, 1> Vector6d;

class LIP{
public:
  enum REF_FRAMES{
                  PELVIS_FRAME,
                  SUPPORT_FOOT_FRAME
  };
  
  LIP(double mass_, double height_, double Tssp_, double Tdsp_):
    mass(mass_), height(height_), Tssp(Tssp_), Tdsp(Tdsp_)
  {
    lambda = std::sqrt(gravity/height);
    sigma1 = lambda * 1/std::tanh(Tssp/2 * lambda);
    sigma2 = lambda * std::tanh(Tssp/2 * lambda);

    Matrix2d Assp; Assp << 0, 1, lambda*lambda, 0;
    Matrix2d Adsp; Adsp << 1, Tdsp, 0, 1;

    Vector2d B; B << -1, 0;

    Amat = (Assp * Tssp).exp() * Adsp;
    Bmat = (Assp * Tssp).exp() * B;

    Kdb << 1, Tdsp + 1/( lambda * std::tanh(Tssp * lambda));

    swing_xTraj = std::make_unique<BezierSupport>(bezierOrder);
    swing_yTraj = std::make_unique<BezierSupport>(bezierOrder);
    swing_zTraj = std::make_unique<BezierSupport>(bezierOrder);
    swing_thetaTraj = std::make_unique<BezierSupport>(bezierOrder);
    com_thetaTraj = std::make_unique<BezierSupport>(bezierOrder);
    bPoints.resize(bezierOrder + 1);

    task = VectorXd::Zero(14);
    twists.resize(2);
    twists[0] = VectorXd::Zero(6); // zero-position and zero-orientation
    twists[1] = VectorXd::Zero(6);

    twist_error.resize(2);
    twist_error[0] = VectorXd::Zero(6);
    twist_error[1] = VectorXd::Zero(6);

    thetaSwingState = VectorXd::Zero(2);
    thetaComState = VectorXd::Zero(2);
    yaw_reference = 0;

    adaptive_dx = std::make_unique<DeltaLearning>(50,2);
    adaptive_dy = std::make_unique<DeltaLearning>(40,4, 1, 5 * 1e-4); // Inputs y, ydot, Vx, Vy
  }

  Vector2d getTimedState(Vector2d &x0, double Time);
  int buildOrbits(double forward_velocity, double lateral_velocity, double u_Left2Right);
  Vector2d getFinalStatesP1();
  Vector2d getFinalStatesP2L();
  Vector2d getFinalStatesP2R();

  int updateModel(double mass_, double height_, double Tssp_, double Tdsp_);

  Vector2d computeStepOffsets(VectorXd &base_velocity, VectorXd &base_position, Vector3d &initialCoM, int support_leg, REF_FRAMES ref_frame, double time, int nSteps, double Yaw, VectorXd &base_ang_velocity);

  /** Computes the trajectory of the LIP CoM [pCoM, vCoM] a 6x1 vector given a sup
   * @Time: Time in seconds
   * @support_leg: 0 is left leg, 1 is right leg  
   * Not implemented, because comTraj is underactuated
   **/
  Vector6d comTraj(double Time, int support_leg, Vector3d comInit, Vector3d vcomInit, Vector3d comRPYInit);

  /** Computes the trajectory of the SwingFoot [pCoM, vCoM] a 6x1 vector, wrt the Pelvis frame, given a support_leg reference.
      @Time: Time in seconds
      @support_leg: 0 is left leg, 1 is right leg  
  **/  
  Vector6d swingTraj(double Time, Vector3d initialFootSupport, Vector3d initialRPY, double deltaX, double deltaY, int support_leg, double Yaw, REF_FRAMES ref_frame);

  double Tssp;
  double Tdsp;
  double lambda;
  double sigma1;
  double sigma2;  
  std::vector<VectorXd> twist_error;

  // Nominal step widths
  double uP1;
  double uP2_L;
  double uP2_R;

  // Final states of CoM
  Vector2d xfP1;
  Vector2d xfP2_L;
  Vector2d xfP2_R;
  Vector2d xfP2;

  // Initial states of CoM
  Vector2d x0P1;
  Vector2d x0P2;
  Vector2d x0P2_L;
  Vector2d x0P2_R;

  double yaw_reference;
  VectorXd thetaSwingState;
  VectorXd thetaComState;

  double height;  
private:
  // LIP parameters
  double mass;
  double gravity = 9.81;

  Matrix2d Amat;
  Matrix2d AmatHalf;
  Vector2d Bmat;
  Matrix2d BmatHalf;
  Vector2d Kdb;

  // LIP velocities
  double forward_velocity;
  double lateral_velocity;

  // Bezier support

  std::unique_ptr<BezierSupport> swing_xTraj;
  std::unique_ptr<BezierSupport> swing_yTraj;
  std::unique_ptr<BezierSupport> swing_zTraj;
  std::unique_ptr<BezierSupport> swing_thetaTraj;
  std::unique_ptr<BezierSupport> com_thetaTraj;    
  const int bezierOrder = 7;
  VectorXd bPoints;
  RobotExpressions rExpr;
  std::vector<VectorXd> twists;
  VectorXd task;
  double deltaTheta;
  std::unique_ptr<DeltaLearning> adaptive_dx;
  std::unique_ptr<DeltaLearning> adaptive_dy;
  

};

#endif
