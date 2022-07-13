#include "IDQPcpp.hpp"
#include "src/Core/Matrix.h"
//#include "constraint.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <ostream>

IDQPcpp::IDQPcpp(int nStates, int nObj, int nActuators, int nContacts) :
  nStates_(nStates), nObj_(nObj), nActuators_(nActuators), nContacts_(nContacts)
{
  int nContactCstr;
  nContactCstr = nContacts * 6;
  nContactCstr_ = nContactCstr;
  dimX = nStates + nActuators + nContactCstr + 1; // +1: CLF relaxation
  ADyn = MatrixXd::Ones(nStates, dimX);
  BDyn = VectorXd::Zero(nStates);

  ACont = MatrixXd::Ones(nContactCstr, dimX);
  BCont = VectorXd::Zero(nContactCstr);

  ACone = MatrixXd::Ones(5 * nContacts, dimX);
  BCone = VectorXd::Zero(5 * nContacts);

  AZMP = MatrixXd::Ones(4, dimX);
  BZMP = VectorXd::Zero(4);

  Atorque = MatrixXd::Ones(2 * nActuators, dimX);
  Btorque = VectorXd::Zero(2 * nActuators);

  Aacc = MatrixXd::Ones(2 * nStates, dimX);
  Bacc = VectorXd::Zero(2 * nStates);

  yStar = VectorXd::Zero(nObj);

  Hopt = MatrixXd::Zero(dimX,dimX);
  fopt = VectorXd::Zero(dimX);

  MatKp = MatrixXd::Zero(nObj,nObj);
  MatKd = MatrixXd::Zero(nObj,nObj);

  Reg = MatrixXd::Identity(dimX, dimX) * 1e-6;
  Reg(dimX, dimX) = 1; // Delta for CLF relaxation
  ZeroVec = VectorXd::Zero(nStates);


  yaOut = VectorXd::Zero(nObj);
  ydOut = VectorXd::Zero(nObj);
  yOut = VectorXd::Zero(nObj);

  DyaOut = VectorXd::Zero(nObj);
  DydOut = VectorXd::Zero(nObj);
  DyOut = VectorXd::Zero(nObj);

  Xsol = VectorXd::Zero(dimX);
  x0 = VectorXd::Zero(dimX);
  
  nEq = 0;
  nIneq = 0;
  dt = 1e-3;

  // CBF
  double w_cbf = 1/(dt * 10);
  double eps_cbf = 3.0;
  Breact = VectorXd::Zero(dimX);
  
  alpha_cbf = eps_cbf/w_cbf + sqrt(eps_cbf*eps_cbf - 1)/w_cbf;
  k_cbf = eps_cbf * w_cbf + sqrt(eps_cbf*eps_cbf - 1) * w_cbf;

  // CLFs
  eps_clf = 1;
  Flin = MatrixXd::Zero(nObj_ * 2, nObj_*2);
  Glin = MatrixXd::Zero(nObj_ * 2, nObj_);
  Pmat = MatrixXd::Zero(nObj_ * 2, nObj_*2);
  Pemat = Pmat;
  Rmat = MatrixXd::Identity(nObj_, nObj_);
  Qmat = MatrixXd::Identity(nObj_ * 2, nObj_ * 2);
  Ie = MatrixXd::Identity(nObj_ * 2, nObj_ * 2);
  Ie.block(0,0,nObj_, nObj_) *= 1/eps_clf;
  eta = VectorXd::Zero(nObj_ * 2);
  Aclf = MatrixXd::Zero(1, nObj_);
  Bclf = VectorXd::Zero(1);

  Flin.block(0,nObj_, nObj_,nObj_) = MatrixXd::Identity(nObj_,nObj_);
  Glin.block(nObj_, 0, nObj_, nObj_) = MatrixXd::Identity(nObj_,nObj_); 

  MotorjointTrans = VectorXd::Ones(nActuators);

  Jac_cont = MatrixXd::Zero(nContactCstr_, nStates_);
  Jac_cont_dot = MatrixXd::Zero(nContactCstr_, nStates_);  
  Bmat = MatrixXd::Zero(nStates_, nActuators_);
  Bmat.block(0,0,6,nActuators_) = MatrixXd::Zero(6,nActuators_);
  Bmat.block(6,0,nActuators_,nActuators_) = MatrixXd::Identity(nActuators_,nActuators_);    

  activeConstraints.dynamics = false;
  activeConstraints.contacts = false;
  activeConstraints.acceleration = false;
  activeConstraints.position = false;
  activeConstraints.friction = false;
  activeConstraints.torque = false;
  activeConstraints.zmp = false;
  activeConstraints.customEq = false;
  activeConstraints.customIneq = false;
}

int IDQPcpp::addDynamicsCMat(std::function<MatrixXd (VectorXd &)> Mmat_, std::function<MatrixXd (VectorXd &, VectorXd &)> Cmat_, std::function<VectorXd (VectorXd &)> Gvec_){

  Mmat = Mmat_;
  Cmat = Cmat_;
  Gvec = Gvec_;

  Cvec = [this] (VectorXd& q, VectorXd& dq) { return Cmat(q,dq) * dq; };

  nEq += nStates_;
  activeConstraints.dynamics = true;
  return 0;
}

int IDQPcpp::addDynamicsCVec(std::function<MatrixXd (VectorXd &)> Mmat_, std::function<VectorXd (VectorXd &, VectorXd &)> Cvec_, std::function<VectorXd (VectorXd &)> Gvec_){

  Mmat = Mmat_;
  Cvec = Cvec_;
  Gvec = Gvec_;

  nEq += nStates_;  
  activeConstraints.dynamics = true;  
  return 0;
}

int IDQPcpp::addDynamicsHvec(std::function< MatrixXd (VectorXd&)> Mmat_,
                          std::function< VectorXd (VectorXd&, VectorXd&)> Hvec_) {
  Mmat = Mmat_;
  Hvec = Hvec_;
  Cvec = [this] (VectorXd& q, VectorXd& dq) { return Hvec(q,dq); };
  Gvec = [this] (VectorXd& q) { return VectorXd::Zero(nStates_); };

  nEq += nStates_;
  activeConstraints.dynamics = true;  
  return 0;
}

int IDQPcpp::addDynamicsHvecNum(MatrixXd Mmat_, VectorXd Hvec_)
{
  MmatNum = Mmat_;
  HvecNum = Hvec_;

  Mmat = [this] (VectorXd &q) {return MmatNum; };
  Cvec = [this] (VectorXd& q, VectorXd& dq) { return HvecNum; };
  Gvec = [this] (VectorXd& q) { return VectorXd::Zero(nStates_); };

  nEq += nStates_;
  activeConstraints.dynamics = true;  
  
  return 0;
}

int IDQPcpp::updateDynamicsHvecNum(MatrixXd Mmat_, VectorXd Hvec_)
{
  MmatNum = Mmat_;
  HvecNum = Hvec_;

  return 0;
}

int IDQPcpp::addContactJacobians(IDQP_type::List_Jac list_Jac, IDQP_type::List_Jac_dot list_Jac_dot){

  List_Jac_cont = list_Jac;
  List_Jac_cont_dot = list_Jac_dot;

  nEq += nContactCstr_;
  activeConstraints.contacts = true;
  return 0;
}

int IDQPcpp::addFrictionCone(double mu_friction)
{
  // Acone X \leq Bcone
  // Friction Cone
  ACone.setZero();
  BCone.setZero();

  for(int i = 0; i < nContacts_; i++) {
    ACone(5 * i + 0, 6 * i + nStates_+nActuators_)   = 1;
    ACone(5 * i + 0, 6 * i + nStates_+nActuators_+2) = -mu_friction/sqrt(2);  
    ACone(5 * i + 1, 6 * i + nStates_+nActuators_)   = -1;
    ACone(5 * i + 1, 6 * i + nStates_+nActuators_+2) = -mu_friction/sqrt(2);  
    ACone(5 * i + 2, 6 * i + nStates_+nActuators_+1) = 1;
    ACone(5 * i + 2, 6 * i + nStates_+nActuators_+2) = -mu_friction/sqrt(2);  
    ACone(5 * i + 3, 6 * i + nStates_+nActuators_+1) = -1;
    ACone(5 * i + 3, 6 * i + nStates_+nActuators_+2) = -mu_friction/sqrt(2);  
    ACone(5 * i + 4, 6 * i + nStates_+nActuators_+2) = -1;
  }

  BCone.setZero();

  nIneq += nContacts_ * 5;
  activeConstraints.friction = true;  
  return 0;
}

int IDQPcpp::addControlLimits(VectorXd &uMin, VectorXd &uMax)
{
  Atorque.setZero();
  Atorque.block(0, nStates_, nActuators_,nActuators_) = MatrixXd::Identity(nActuators_, nActuators_);
  Atorque.block(nActuators_, nStates_, nActuators_, nActuators_) = -MatrixXd::Identity(nActuators_, nActuators_);

  VectorXd uMin_trans(nActuators_);
  VectorXd uMax_trans(nActuators_);

  for(int i = 0; i < nActuators_; i++) {
    uMin_trans(i) = uMin(i) / MotorjointTrans(i);
    uMax_trans(i) = uMax(i) / MotorjointTrans(i);
  }

  Btorque << uMax_trans, -uMin_trans;

  nIneq += 2 * nActuators_;

  activeConstraints.torque = true;
  return 0;
}

int IDQPcpp::addPositionLimits(VectorXd &qMin, VectorXd &qMax, VectorXd &IndexMap)
{
  // Aposition X <= [bmin;bmax]
  PositionIndexMap = IndexMap;
  qmin = qMin;
  qmax = qMax;
  
  Aposition = MatrixXd::Zero(2 * IndexMap.size(), dimX);
  for(int i = 0; i < IndexMap.size(); i++) {
    Aposition(i, IndexMap(i)) = -alpha_cbf;
    Aposition(i + IndexMap.size(), IndexMap(i)) = alpha_cbf;
  }

  Bposition = VectorXd::Zero(2 * IndexMap.size());  

  nIneq += 2 * IndexMap.size();
  activeConstraints.position = true;
  return 0;
}

int IDQPcpp::setAccelerationLimits(VectorXd &qddMin, VectorXd &qddMax)
{
  Aacc.setZero();
  Aacc.block(0, 0,nStates_,nStates_) = MatrixXd::Identity(nStates_, nStates_);
  Aacc.block(nStates_, 0, nStates_, nStates_) = -MatrixXd::Identity(nStates_, nStates_);

  Bacc << qddMax, -qddMin;
  nIneq += 2 * nStates_;

  activeConstraints.acceleration = true;
  return 0;
}

int IDQPcpp::addZMP(double LengthXmin, double LengthXmax, double LengthYmin, double LengthYmax, List_Ad Ad_foot, List_Position Pos_foot)
{
  // Populate List of Adjoint Matrices
  List_Ad_contacts = Ad_foot;
  List_Pos_contacts = Pos_foot;
  // We will Transform into unique Reference frame Lambda_ref = \sum Ad_c * Lambda_c
  // during the solve phase. In the reference frame, we will solve the ZMP with AZMPr
  // and then project back into AZMP
  AZMPr = MatrixXd::Zero(4, 6);

  AZMPr(0, 2) = -LengthXmax;
  AZMPr(0, 4) = -1;

  AZMPr(1, 2) = LengthXmin;
  AZMPr(1, 4) = 1;

  AZMPr(2, 2) = -LengthYmax;
  AZMPr(2, 3) = 1;

  AZMPr(3, 2) = LengthYmin;
  AZMPr(3, 3) = -1;
  
  BZMP.setZero();
  
  nIneq += 4; // It will have its own qp constraint
  activeConstraints.zmp = true;    
  return 0;
}

int IDQPcpp::addCustomEq(MatrixXd &Acustom, VectorXd &Bcustom)
{
  Acustom_eq = Acustom;
  Bcustom_eq = Bcustom;

  activeConstraints.customEq = true;

  nEq += Acustom.rows();
  
  return 0;
}

int IDQPcpp::addCustomIneq(MatrixXd &Acustom, VectorXd &Bcustom)
{
  Acustom_ineq = Acustom;
  Bcustom_ineq = Bcustom;

  nIneq += Acustom.rows();
  activeConstraints.customIneq = true;
  return 0;
}


int IDQPcpp::addTypeIconstraint(List_Jac M_I, int nConst_I)
{
  List_MI = M_I;
  int nRelations = M_I.size();
  
  M_I_mat = MatrixXd::Zero(nConst_I * nRelations, dimX);
  b_I_vec = VectorXd::Zero(nConst_I * nRelations);

  this->nConstI = nConst_I;
  nEq += nConst_I * nRelations;
  activeConstraints.typeIEq = true;
  return 0;
}


int IDQPcpp::addTypeIIconstraint(List_Jac M_II, List_Jac_dot M_II_dot, int nConst_II)
{
  List_MII = M_II;
  List_MII_dot = M_II_dot;
  
  int nRelations = M_II.size();
  
  M_II_mat = MatrixXd::Zero(nConst_II * nRelations, dimX);
  b_II_vec = VectorXd::Zero(nConst_II * nRelations);

  this->nConstII = nConst_II;
  nEq += nConst_II * nRelations;
  activeConstraints.typeIIEq = true;
  return 0;
}

int IDQPcpp::specifyMotorTransmissions(VectorXd &transmissions)
{
  MotorjointTrans = transmissions;
  return 0;
}

int IDQPcpp::addObjectives(std::function<VectorXd (VectorXd &)> y_a,
                        std::function<MatrixXd (VectorXd &)> Jac_y,
                        std::function<MatrixXd (VectorXd &, VectorXd&)> Jac_y_dot,
                        std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> y_d,
                        std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> Dy_d,                        
                        VectorXd &Kp, VectorXd &Kd)
{
  ya = y_a;
  yd = y_d;
  Dyd = Dy_d;
  DDyd = [this] (double t, VectorXd&, VectorXd&, VectorXd&, VectorXd&) { return VectorXd::Zero(nObj_); };
  this->Jac_y = Jac_y;
  this->Jac_y_dot = Jac_y_dot;

  MatKp = Kp.array().matrix().asDiagonal();
  MatKd = Kd.array().matrix().asDiagonal();

  return 0;
}

int IDQPcpp::addObjectives(std::function<VectorXd (VectorXd &)> y_a,
                        std::function<MatrixXd (VectorXd &)> Jac_y,
                        std::function<MatrixXd (VectorXd &, VectorXd&)> Jac_y_dot,
                        std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> y_d,
                        std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> Dy_d,
                        std::function<VectorXd (double, VectorXd&, VectorXd&, VectorXd&, VectorXd&)> DDy_d,                        
                        VectorXd &Kp, VectorXd &Kd)
{
  ya = y_a;
  yd = y_d;
  Dyd = Dy_d;
  DDyd = DDy_d;
  this->Jac_y = Jac_y;
  this->Jac_y_dot = Jac_y_dot;

  MatKp = Kp.array().matrix().asDiagonal();
  MatKd = Kd.array().matrix().asDiagonal();

  return 0;
}

int IDQPcpp::setRegularizationTerm(double q_ddot, double u_control, double lambda)
{
  Reg.block(0,0,nStates_,nStates_) =  MatrixXd::Identity(nStates_, nStates_) * q_ddot;
  Reg.block(nStates_,nStates_,nActuators_,nActuators_) = MatrixXd::Identity(nActuators_, nActuators_) * u_control;
  Reg.block(nStates_ + nActuators_, nStates_ + nActuators_,
            nContactCstr_, nContactCstr_) = MatrixXd::Identity(nContactCstr_, nContactCstr_) * lambda;
  return 0;
}

MatrixXd IDQPcpp::solveRiccatiArimotoPotter() {

  MatrixXd Px;
  const uint dim_x = Flin.rows();
  const uint dim_u = Glin.cols();

  // set Hamilton matrix
  Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
  Ham << Flin, -Glin * Rmat.inverse() * Glin.transpose(), -Qmat, -Flin.transpose();
  std::cout<<"Ham: "<<Ham<<std::endl;
  // calc eigenvalues and eigenvectors
  Eigs = new Eigen::EigenSolver<Eigen::MatrixXd>(Ham);
  std::cout<<"Eigs defined"<<std::endl;
  // check eigen values
  std::cout << "eigen values：\n" << Eigs->eigenvalues() << std::endl;
  std::cout << "eigen vectors：\n" << Eigs->eigenvectors() << std::endl;

  // extract stable eigenvectors into 'eigvec'
  Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);

  int j = 0;
  for (int i = 0; i < 2 * dim_x; ++i) {
    if (Eigs->eigenvalues()[i].real() < 0.) {
      eigvec.col(j) = Eigs->eigenvectors().block(0, i, 2 * dim_x, 1);
      ++j;
    }
  }
  std::cout<<"eigvec: "<<eigvec<<std::endl;  

  // calc P with stable eigen vector matrix
  Eigen::MatrixXcd Vs_1, Vs_2;
  Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
  Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
  Px = (Vs_2 * Vs_1.inverse()).real();

  return Px;

}

int IDQPcpp::buildOptimization()
{
  std::cout<<"Building Optimization"<<std::endl;
  std::cout<<"nObj: "<<nObj_<<std::endl;
  std::cout<<"Flin: "<<Flin<<std::endl;
  std::cout<<"Glin: "<<Glin<<std::endl;
  std::cout<<"Qmat: "<<Qmat<<std::endl;
  std::cout<<"Rmat: "<<Rmat<<std::endl;
  // Compute the RES-CLF constraint!
  //Pmat = solveRiccatiArimotoPotter();
  //Pemat = Ie * Pmat * Ie;
  // Hard Coding (TODO: Debug Eigen Error on Eigensolver call)
  Pmat.block<6,6>(0,0) = MatrixXd::Identity(nObj_ * 2, nObj_ * 2) * 10.0995;
  Pmat.block<6,6>(6,6) = MatrixXd::Identity(nObj_ * 2, nObj_ * 2) * 1.00995;  
  Pmat.block<6,6>(0,6) = MatrixXd::Identity(6,6) * 0.001;
  Pmat.block<6,6>(6,0) = MatrixXd::Identity(6,6) * 0.001;

  Pemat = Ie * Pmat * Ie;
  
  // VectorXcd eivals = Pemat.eigenvalues();
  // std::cout<<"computing largest"<<std::endl;      
  // VectorXd largest = eivals.real();
  // std::cout<<"c_clf"<<std::endl;    
  c_clf = 1/(eps_clf  *  183.205 );


  std::cout<<"Optimization Built"<<std::endl;
  
  IneqCounter = 0;
  nIneq +=1; // Include CLF condition
  
  // It will start to populate in order:
  // Friction Cone Constraints
  // Torque Limit Constraints
  // ZMP Constraints
  // Acceleration Constraints
  // If they are available
  
  Aeq = MatrixXd::Zero(nEq, dimX);

  Beq = VectorXd::Zero(nEq);

  Aineq = MatrixXd::Zero(nIneq, dimX);
  Bineq = VectorXd::Zero(nIneq);

  if( activeConstraints.friction ) {
    Aineq.block(0,0, nContacts_ * 5, dimX) = ACone;
    Bineq.segment(0, nContacts_ * 5) = BCone;
    IneqCounter += nContacts_ * 5;
  }

  if( activeConstraints.torque ) {
    Aineq.block(IneqCounter, 0, nActuators_ * 2, dimX) = Atorque;
    Bineq.segment(IneqCounter, nActuators_ * 2) = Btorque;
    IneqCounter += nActuators_ * 2;
  }

  if( activeConstraints.acceleration ) {
    Aineq.block(IneqCounter, 0, nStates_ * 2, dimX) = Aacc;
    Bineq.segment(IneqCounter, nStates_ * 2) = Bacc;
    IneqCounter += nStates_ * 2;
  }

  if( activeConstraints.customIneq ) {
    //qp.addConstraint( lessThan(par(Acustom_ineq) * X, par(Bcustom_ineq) ) );
    Aineq.block(IneqCounter, 0, Acustom_ineq.rows(), dimX) = Acustom_ineq;
    Bineq.segment(IneqCounter, Acustom_ineq.rows()) = Bcustom_ineq;
    IneqCounter += Acustom_ineq.rows();
  }

  
  qp.reset(dimX, nEq, nIneq);
  //solver->setAlpha(1.2); //1.8
  //solver->setPolish(true);    
  //solver->setEpsAbs(1e-2); //1e-10
  //solver->setEpsRel(1e-2); //1e-10
  //solver->setMaxIter(15);
  return 0;
}

VectorXd IDQPcpp::solve(double t, VectorXd &q, VectorXd &dq, VectorXd &qdes, VectorXd &dqdes)
{
  // It will start to populate in order:
  // Dynamic Constraints
  // Contact Constraints
  // If they are available
  EqCounter = 0;
  dyn_IneqCounter = IneqCounter;


  if( activeConstraints.dynamics ) {
    // Dynamics
    Aeq.block(0,0,nStates_,nStates_) = Mmat(q);
    Aeq.block(0, nStates_, nStates_, nActuators_) = -Bmat;
    if( !activeConstraints.contacts ){
      Aeq.block(0, nStates_ + nActuators_, nStates_, nContactCstr_) = -Jac_cont.transpose();
    }else{
      for(int nc = 0; nc < nContacts_; nc++) {
        Jac_cont.block(6*nc,0,6, nStates_) = List_Jac_cont[nc](q);
        Jac_cont_dot.block(6*nc,0,6,nStates_) = List_Jac_cont_dot[nc](q,dq);
      }
      Aeq.block(0, nStates_ + nActuators_, nStates_, nContactCstr_) = -Jac_cont.transpose();
    }  
    Beq.segment(0,nStates_) = -Cvec(q,dq) - Gvec(q);
    EqCounter += nStates_;
  }

  //  Contact Constraints
  if( activeConstraints.contacts ) {
    Aeq.block(EqCounter, 0, nContactCstr_, nStates_) = Jac_cont;
    Aeq.block(EqCounter, nStates_, nContactCstr_, dimX - nStates_) = MatrixXd::Zero(nContactCstr_, dimX - nStates_);

    //Beq.segment(0,nStates_) = -Cmat(q,dq) * dq - Gvec(q);
    Beq.segment(EqCounter, nContactCstr_) = - Jac_cont_dot * dq;
    EqCounter += nContactCstr_;
  }

  if( activeConstraints.customEq ) {
    Aeq.block(EqCounter, 0, Acustom_eq.rows(), dimX) = Acustom_eq;
    Beq.segment(EqCounter, Acustom_eq.rows()) = Bcustom_eq;
    EqCounter += Acustom_eq.rows();
  }


  //Aeq =  (Aeq.array().abs()<1e-1).select(0,Aeq);

  //Aeq_sp = Aeq.sparseView();

  // Position Constraints
  if(activeConstraints.position) {
    for(int i = 0; i < PositionIndexMap.size(); i++) {
      Bposition(i) = k_cbf * (q(PositionIndexMap(i)) - qmin(i)) + (1 + k_cbf * alpha_cbf) * dq(PositionIndexMap(i));
      Bposition(PositionIndexMap.size() + i) = k_cbf * (qmax(i) - q(PositionIndexMap(i))) - (1 + k_cbf * alpha_cbf) * dq(PositionIndexMap(i));
      Breact(PositionIndexMap(i)) = 1e-8 * (-alpha_cbf) / ((q(PositionIndexMap(i)) - qmin(i)) + alpha_cbf * dq(PositionIndexMap(i)) + 0.01);
    }

    Aineq.block(dyn_IneqCounter, 0, Aposition.rows(), dimX) = Aposition;
    Bineq.segment(dyn_IneqCounter, Aposition.rows()) = Bposition;
    dyn_IneqCounter += Aposition.rows();
  }

  // ZMP Constraints
  if( activeConstraints.zmp ) {
    MatrixXd AZMP_delta = MatrixXd::Zero(4,6);
    VectorXd Avg_Pos = VectorXd::Zero(3);

    for(int i = 0; i < nContacts_; i++) {
      Avg_Pos += List_Pos_contacts[i](q)/nContacts_;        
    }

    // Making sure the adjoint matrices are wrt foot contact plane
    // MatrixXd AdPlane = MatrixXd::Identity(6,6);
    // AdPlane(3, 1) = -(-q(2)) ;
    // AdPlane(5, 0) = (-q(2));

    AZMP_delta(0, 2) = -Avg_Pos[0];
    AZMP_delta(1, 2) = Avg_Pos[0];
    AZMP_delta(2, 2) = -Avg_Pos[1];
    AZMP_delta(3, 2) = Avg_Pos[1];
      
    for(int i = 0; i < nContacts_; i++) {
      AZMP.block(0, nStates_ + nActuators_ + 6 * i, 4, 6) = (AZMPr + AZMP_delta) * List_Ad_contacts[i](q);
    }

    Aineq.block(dyn_IneqCounter, 0, 4, dimX) = AZMP;
    Bineq.segment(dyn_IneqCounter, 4) = BZMP;
    dyn_IneqCounter += 4;
  }
  
  // Type I and II constraints
  if( activeConstraints.typeIEq ) {
    for(int i = 0; i < List_MI.size(); i++) {
      M_I_mat.block(nConstI*i, 0, this->nConstI, dimX ) = List_MI[i](q);
    }

    Aeq.block(EqCounter, 0, nConstI * List_MI.size(), dimX) = M_I_mat;
    Beq.segment(EqCounter, nConstI * List_MI.size()) = VectorXd::Zero(nConstI * List_MI.size());

    EqCounter += nConstI * List_MI.size();
  }

  if( activeConstraints.typeIIEq ) {
    for(int i = 0; i < List_MII.size(); i++) {
      M_II_mat.block(nConstII*i, 0, this->nConstII, nStates_) = List_MII[i](q);
      b_II_vec.segment(nConstII*i, this->nConstII) = -List_MII_dot[i](q, dq) * dq;
    }

    Aeq.block(EqCounter, 0, nConstII * List_MII.size(), dimX) = M_II_mat;
    Beq.segment(EqCounter, nConstII * List_MII.size()) = b_II_vec;
    EqCounter += nConstII * List_MII.size();
  }
  
  // Cost Function
  yaOut = ya(q);
  ydOut = yd(t, q, dq, qdes, dqdes);

  Jacy = Jac_y(q);
  DJacy = Jac_y_dot(q,dq);
  
  DyaOut = Jacy * dq;
  DydOut = Dyd(t,q,dq, qdes, dqdes);
  DDydOut = DDyd(t, q, dq,qdes, dqdes);  
  
  yOut = yaOut - ydOut;
  DyOut = DyaOut - DydOut;
  
  yStar = DDydOut - MatKp * yOut - MatKd * (DyOut);

  // CLF Constraints
  eta << yOut, DyOut;
  Aclf = 2 * eta.transpose() * Pemat * Glin * Jacy;
  Bclf = -eta.transpose() * (Flin.transpose() * Pemat + Pemat * Flin + c_clf * Pemat) * eta
    - 2 * eta.transpose() * Pemat * Glin * DJacy * dq;

  Aineq.block(dyn_IneqCounter, 0, 1, nStates_) = Aclf;
  Aineq(dyn_IneqCounter, dimX) = 1; // Consider delta for CLF relaxation
  Bineq.segment(dyn_IneqCounter, 1) = Bclf;

  Hopt = Reg; // Regularization term  
  //Hopt.block(0,0, nStates_, nStates_) += Jacy.transpose() * Jacy;
  //Hopt = 0.5*(Hopt + Hopt.transpose());

  //fopt.segment(0, nStates_) = 2 * dq.transpose() * DJacy.transpose() * Jacy
  //  - 2 * yStar.transpose() * Jacy;		 
  //fopt += Breact;  

  VectorXd x0_guess = x0;
  solutionStatus = qp.solve_quadprog(2 * Hopt, fopt, -Aeq, Beq, -Aineq, Bineq, x0_guess);

  if(solutionStatus != 0) {
    Xsol = x0;
  }else{
    x0 = x0_guess;
    Xsol = x0_guess;
  }

  // if( solver->getInfo().status_val != 1) {
  //   std::cout<<"Not solved, error status: "<<solver->getInfo().status_val<< std::endl;
  //   solutionStatus = 0;
    // Check if some constraints failed
    Eigen::ArrayXd test=  (Aineq * Xsol - Bineq);
    if( ( test > 1e-3).any()) {
      std::cout<< "Inequality violation" <<std::endl;
    }

    Eigen::ArrayXd test2 =  (Aeq * Xsol - Beq);
    if( ( abs(test2) > 1e-3).any()) {
      std::cout<< "Equality violation" <<std::endl;
    }

    Eigen::ArrayXd test3 =  (Aineq.block(nIneq,0,1,dimX) * Xsol - Bclf);
    if( ( test > 1e-3).any()) {
      std::cout<< "CLF violation" <<std::endl;
    }

    Eigen::ArrayXd test4 =  (M_II_mat * Xsol - b_II_vec);
    if( ( test > 1e-3).any()) {
      std::cout<< "CLF violation" <<std::endl;
    }    
  //  if(activeConstraints.zmp) {
  //     Eigen::ArrayXd test3 =  (AZMP * Xsol - BZMP);
  //     if( ( test3 > 1e-3).any()) {
  //       std::cout<< "ZMP Violation" <<std::endl;
  //     }      
  //   }
  // }else{
  //   solutionStatus = 1;
  // }


  // // Check jacobians Constraints and outputs are full rank
  // MatrixXd JacTotal(18,18);
  // JacTotal << List_Jac_cont[0](q), List_Jac_cont[1](q), Jac_y(q);
  // Eigen::ColPivHouseholderQR< MatrixXd > RankMatrix(JacTotal);
  // std::cout<<"Rank : "<< RankMatrix.rank() << std::endl;

  // if( RankMatrix.rank() != 18) {
  //   std::cout<<" ---------------------- Singularity ------------------ " <<std::endl;
  // }
  return Xsol;
}

int IDQPcpp::testEvaluation(VectorXd &q, VectorXd &dq) {
  std::cout<<"Mmat: "<<Mmat(q)<<std::endl;
  std::cout<<"Cmat: "<<Cvec(q,dq)<<std::endl;
  std::cout<<"Gvec: "<<Gvec(q)<<std::endl;
  std::cout<<"Contact Jacobian 1: " <<  List_Jac_cont[0](q) << std::endl;
  std::cout<<"Contact Jacobian 2: " <<  List_Jac_cont[1](q) << std::endl;
  std::cout<<"Contact Jacobian dot 1: " <<  List_Jac_cont_dot[0](q, dq) << std::endl;
  std::cout<<"Contact Jacobian dot 2: " <<  List_Jac_cont_dot[1](q, dq) << std::endl;

  std::cout<<"yOutput: " <<  ya(q) << std::endl;
  std::cout<<"JyOutput: " <<  Jac_y(q)<< std::endl;
  std::cout<<"DJyOutput: " <<  Jac_y_dot(q,dq) << std::endl;

  std::cout<<"Aeq: "<<Aeq<<std::endl;
  // Get Rank of Aeq
  FullPivLU<MatrixXd> lu_decomp(Aeq);
  auto rank = lu_decomp.rank();
  std::cout<<"Aeq rank: "<<rank<<std::endl;
  std::cout<<"Aeq rows: "<<Aeq.rows()<<std::endl;
  std::cout<<"nEQ: "<<nEq<<std::endl;
  std::cout<<"Aineq: "<<Aineq<<std::endl;
  return 0;
}


VectorXd IDQPcpp::getZMPprediction(VectorXd &q) {
  VectorXd act_ZMP(2);
  VectorXd lambda_prj(6);
  VectorXd avgPos(3);  
  
  lambda_prj.setZero();
  avgPos.setZero();

  if(activeConstraints.zmp) {
    for(int i = 0; i < nContacts_; i++) {
      lambda_prj += List_Ad_contacts[i](q)*Xsol.segment<6>(nStates_ + nActuators_ + i * 6);
      avgPos += (List_Pos_contacts[i](q))/nContacts_;
    
    }

    if(lambda_prj(2) != 0) {
      act_ZMP[0] = -lambda_prj(4)/lambda_prj(2) - avgPos(0);
      act_ZMP[1] = lambda_prj(3)/lambda_prj(2) - avgPos(1);
    }
  }
  
  return act_ZMP;
}

VectorXd IDQPcpp::get_qdd() { return Xsol.segment(0, nStates_); }

VectorXd IDQPcpp::get_torque()
{
  VectorXd torque_idqp(nActuators_);

  // Apply transmission
  for(int i = 0; i < nActuators_; i++)
    {
      torque_idqp(i) = Xsol(nStates_ + i) * MotorjointTrans(i);
    }
  
  return torque_idqp;
}
