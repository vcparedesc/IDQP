#include "motion_generator.hpp"
#include "math_utils.hpp"
#include "robot_expressions.hpp"

double sech(double x) { return 1 / std::cosh(x); }

int LIP::updateModel(double mass_, double height_, double Tssp_, double Tdsp_)
{
  mass = mass_;
  height = height_;
  Tssp = Tssp_;
  Tdsp = Tdsp_;
  
  lambda = std::sqrt(gravity/height);
  sigma1 = lambda * 1/std::tanh(Tssp/2 * lambda);
  sigma2 = lambda * std::tanh(Tssp/2 * lambda);

  Matrix2d Assp; Assp << 0, 1, lambda*lambda, 0;
  Matrix2d Adsp; Adsp << 1, Tdsp, 0, 1;

  Vector2d B; B << -1, 0;

  Amat = (Assp * Tssp).exp() * Adsp;
  Bmat = (Assp * Tssp).exp() * B;
  Kdb << 1, Tdsp + 1/( lambda * std::tanh(Tssp * lambda));
  
  return 0;
}

Vector2d LIP::getTimedState(Vector2d &x0, double Time)
{
  Vector2d state;
  double phi = lambda * Time;

  state << std::cosh(phi) * x0(0) + 1/lambda * std::sinh(phi) * x0(1),
    lambda * std::sinh(phi) * x0(0) + std::cosh(phi) * x0(1);

  return state;
}

Vector2d LIP::getFinalStatesP1()
{
  Vector2d state;
  state << uP1 / (2 + Tdsp*sigma1), ( uP1 / (2 + Tdsp*sigma1) ) * sigma1;
  return state;
}

Vector2d LIP::getFinalStatesP2L()
{
  Vector2d state;
  double D = (uP2_L + uP2_R)/2;
  double d2 = std::pow(lambda * sech(lambda/2 * Tssp), 2) * D /
    (lambda * lambda * Tdsp + 2 * sigma2);
  state << (uP2_L - Tdsp * d2)/(2 + Tdsp * sigma2),
    (uP2_L - Tdsp * d2)/(2 + Tdsp * sigma2) * sigma2 + d2;

  return state;
}

Vector2d LIP::getFinalStatesP2R()
{
  Vector2d state;
  double D = (uP2_L + uP2_R)/2;
  double d2 = std::pow(lambda * sech(lambda/2 * Tssp), 2) * D /
    (lambda * lambda * Tdsp + 2 * sigma2);
  state << (uP2_R - Tdsp * d2)/(2 + Tdsp * sigma2),
    (uP2_R - Tdsp * d2)/(2 + Tdsp * sigma2) * sigma2 + d2;

  return state;
}

int LIP::buildOrbits(double forward_velocity_, double lateral_velocity_, double u_Left2Right)
{
  forward_velocity = forward_velocity_;
  lateral_velocity = lateral_velocity_;

  uP1 = forward_velocity_ * (Tssp + Tdsp);
  uP2_L = u_Left2Right;
  uP2_R = 2 * lateral_velocity_ * (Tssp + Tdsp) - u_Left2Right;

  xfP1 = getFinalStatesP1();
  xfP2_L = getFinalStatesP2L();
  xfP2_R = getFinalStatesP2R();

  x0P1 << -xfP1(0), xfP1(1);
  x0P2_L = getTimedState(x0P2_L, -Tssp);
  x0P2_R = getTimedState(x0P2_R, -Tssp);
  
  return 0;
}

Vector6d LIP::comTraj(double Time, int support_leg, Vector3d comInit, Vector3d vcomInit, Vector3d comRPYInit) // Works for support trajectory as well
{
  double z0 = comInit(2);
  double Time_;
  Vector6d comState = Vector6d::Zero();
  Vector2d x0P1_actual;
  Vector2d x0P2_actual;

  double YawInit = comRPYInit(2);
  double thetaf = 0;
  
  Time_ = (Time > Tssp) ? Tssp : Time;
  Time_ = (Time_ < 0) ? 0 : Time_;

  double up_fast = 1 - exp(-5 * Time_/Tssp);
  double dt_up_fast = 5 * exp(-5 * Time_/Tssp);
  
  //double zt = z0 + (-height - z0) * Time_/Tssp;
  //double vzt = (-height - z0) * 1/Tssp;
  double w = 2 * 3.14156 / Tssp;
  double zt = z0 + (-height - z0) * Time_/Tssp;  
  double vzt = (-height - z0) * 1/Tssp;

  // WRT Support foot
  x0P1_actual << comInit(0), vcomInit(0);
  x0P2_actual << comInit(1), vcomInit(1);

  if(support_leg == 0) {
    x0P2 = x0P2_L;
  }else if(support_leg == 1) {
    x0P2 = x0P2_R;
  }

  Vector2d xt_vec =  getTimedState(x0P1, Time_);
  Vector2d yt_vec =  getTimedState(x0P2, Time_);

  bPoints << YawInit, YawInit, YawInit, YawInit, thetaf, thetaf, thetaf, thetaf;
  com_thetaTraj->updatePoints(bPoints);
  thetaComState(0) = com_thetaTraj->Bezier(Time_/Tssp);
  thetaComState(1) = com_thetaTraj->DBezier(Time_/Tssp) / Tssp;    

  if(Time_ >= Tssp) {
    vzt *= 0;
  }  
  comState << x0P1_actual(0) + up_fast * (xt_vec(0) - x0P1_actual(0)),
    x0P2_actual(0) + up_fast * (yt_vec(0) - x0P2_actual(0)),
    zt,
    x0P1_actual(1) + up_fast * (xt_vec(1) - x0P1_actual(1)),
    x0P2_actual(1) + up_fast * (yt_vec(1) - x0P2_actual(1)),    
    vzt;
  return comState;
}

// initialFootSupport: initial swing position
Vector6d LIP::swingTraj(double Time, Vector3d initialFootSupport, Vector3d initialRPY, double deltaX, double deltaY, int support_leg, double Yaw, REF_FRAMES ref_frame)
{
  Vector6d swingState;
  double swingHeight = 0.12;//0.05;
  double offset = 0.000;
  double w = (3.14159 + offset)/Tssp;
  double xf = 0;
  double yf = 0;
  double zf = 0;
  double thetaf = 0;
  Vector3d pswing;
  Vector3d vswing;
  VectorXd bPoints(8);
  double mult_frame = 0;
  static double YawInit = 0;

  double Time_;

  //Time = Time * 1/Tssp;
  Time_ = (Time > Tssp) ? Tssp : Time;
  Time_ = (Time_ < 0) ? 0 : Time_;  

  if(support_leg == 0) { // LEFT LEG
    if(ref_frame == REF_FRAMES::PELVIS_FRAME) {
      xf = uP1/2.0 + deltaX * (1 - exp(-3 * (Time_ - Tssp/2.0)));
      yf = uP2_L/2.0 + deltaY;
    }else if(ref_frame == REF_FRAMES::SUPPORT_FOOT_FRAME) {
      xf = uP1 + deltaX * (1 - exp(-3 * (Time_)));
      yf = uP2_L + deltaY  * (1 - exp(-3 * (Time_)));
    }
  }else if(support_leg == 1) { // RIGHT LEG
    if(ref_frame == REF_FRAMES::PELVIS_FRAME) {
      xf = uP1/2.0 + deltaX * (1 - exp(-3 * (Time_ - Tssp/2.0)));
      yf = uP2_R/2.0 + deltaY;
    }else if(ref_frame == REF_FRAMES::SUPPORT_FOOT_FRAME) {
      xf = uP1 + deltaX  * (1 - exp(-3 * (Time_)));
      yf = uP2_R + deltaY * (1 - exp(-3 * (Time_)));
    }
  }

  if(ref_frame == REF_FRAMES::PELVIS_FRAME){
    mult_frame = 1;
  }else if(ref_frame == REF_FRAMES::SUPPORT_FOOT_FRAME) {
    mult_frame = 0;
  }

  YawInit = initialRPY(2);

  thetaf = deltaTheta;

  bPoints << initialFootSupport(0), initialFootSupport(0), xf, xf, xf, xf, xf, xf;
  swing_xTraj->updatePoints(bPoints);
  pswing(0) = swing_xTraj->Bezier(Time_/Tssp);
  vswing(0) = swing_xTraj->DBezier(Time_/Tssp)/Tssp;

  bPoints << initialFootSupport(1), initialFootSupport(1), (yf-initialFootSupport(1))/2, (yf-initialFootSupport(1))/2, yf, yf, yf, yf;
  swing_yTraj->updatePoints(bPoints);
  pswing(1) = swing_yTraj->Bezier(Time_/Tssp);
  vswing(1) = swing_yTraj->DBezier(Time_/Tssp)/Tssp;

  //bPoints << initialFootSupport(2), relHeight, relHeight, relHeight,  relHeight, zf, zf, zf;
  bPoints << initialFootSupport(2), swingHeight/2, swingHeight, swingHeight, swingHeight/1.5, swingHeight/2,-offset,-offset;
  swing_zTraj->updatePoints(bPoints);
  pswing(2) = swing_zTraj->Bezier(Time_/Tssp);
  vswing(2) = swing_zTraj->DBezier(Time_/Tssp)/Tssp;

  bPoints << YawInit, YawInit, YawInit, YawInit, thetaf, thetaf, thetaf, thetaf;
  swing_thetaTraj->updatePoints(bPoints);
  thetaSwingState(0) = swing_thetaTraj->Bezier(Time_/Tssp);
  thetaSwingState(1) = swing_thetaTraj->DBezier(Time_/Tssp) / Tssp;
  

  if(Time_ >= Tssp) {
    vswing *= 0;
    thetaSwingState(1) = 0;
  }

  swingState.head(3) = pswing;
  swingState.tail(3) = vswing;

  return swingState;
}

// base_position: currentFootsupport
Vector2d LIP::computeStepOffsets(VectorXd &base_velocity, VectorXd &base_position, Vector3d &initialCoM, int support_leg, REF_FRAMES ref_frame,  double time, int nSteps, double Yaw, VectorXd &base_ang_velocity){
  Vector2d offset;
  Vector2d xCoM_state_init;
  Vector2d yCoM_state_init;
  Vector2d xCoM_state_end;
  Vector2d yCoM_state_end;

  static double YawInit = 0;
  static double vTheta_int_error = 0;

  double Ymax = 0.35;

  VectorXd dx_inputs(2);
  double dx_error;
  VectorXd dy_inputs(4);
  double dy_output;
  double dy_error;

  if(time < 0.02) {
    YawInit = Yaw;
  }

  // We need to compute the predicted pre-impact states of the CoM
  if( ref_frame == SUPPORT_FOOT_FRAME ) {
    xCoM_state_init << initialCoM(0), base_velocity(0);
    yCoM_state_init << initialCoM(1), base_velocity(1);

    xCoM_state_end = getTimedState(xCoM_state_init, (Tssp - time) * 0.95);
    yCoM_state_end = getTimedState(yCoM_state_init, (Tssp - time) * 0.95); 
  }else if( ref_frame == PELVIS_FRAME ) {
    // xCoM_state_init << -initialCoM(0), base_velocity(0);
    // yCoM_state_init << -initialCoM(1), base_velocity(1);
    xCoM_state_init << -base_position(0), base_velocity(0);
    yCoM_state_init << -base_position(1), base_velocity(1);

    xCoM_state_end = getTimedState(xCoM_state_init, (Tssp - time));    
    yCoM_state_end = getTimedState(yCoM_state_init, (Tssp - time));
  }



  dy_inputs.normalize();

  // WRT support foot
  if(support_leg == 0) { // Left leg
    x0P2 = x0P2_L;
    xfP2 = xfP2_L;
  }else if(support_leg == 1) {
    x0P2 = x0P2_R;
    xfP2 = xfP2_R;
  }

   Kdb <<0, 0.9;//0.4;
   offset(0) = Kdb.dot(xCoM_state_end - xfP1);
   offset(1) = Kdb.dot(yCoM_state_end - xfP2);

  // Imposing limits on minimum foot separation
  // Y direction
  if(support_leg == 0) { // LEFT
    if( offset(1) >= -uP2_L/2.0  +base_position(1) - 0.15) {
      offset(1) = -uP2_L/2.0  +base_position(1) - 0.15;
    }
    // Maximum Y separation
    if( offset(1) <= -uP2_L/2.0  +base_position(1) - Ymax) {
      offset(1) = -uP2_L/2.0  +base_position(1) - Ymax;
    }
    
  }else if(support_leg == 1) { // RIGHT
    if( offset(1) <= -uP2_R/2.0  +base_position(1) + 0.15) {
      offset(1) = -uP2_R/2.0  +base_position(1) + 0.15;
    }
    // Maximum Y separation
    if( offset(1) >= -uP2_R/2.0  +base_position(1) + Ymax) {
      offset(1) = -uP2_R/2.0  +base_position(1) + Ymax;
    }    
  }



  // Imposing limits on max foot separation
  // X direction
  if( offset(0) >= -uP1/2.0  +base_position(0) + 0.40) {
    offset(0) = -uP1/2.0  +base_position(0) + 0.40;
  }

  vTheta_int_error += (yaw_reference - Yaw) / 1e-3;
  
  deltaTheta = 0.2 * (yaw_reference - Yaw) + 0.08 * (0 - base_ang_velocity(2))+ vTheta_int_error * 1e-9 ; 
  
  if(support_leg == 0) {
    // Nothing
    //offset(1) *= -1;        
  }else if(support_leg == 1) {
    //offset(1) *= -1;    
    }


  return offset;
}
