# README

This repo contains information of the QP based inverse dynamic controllers. For now there are two types of controlles:

- Task Space PD control

- CLFQP

The structure of the equation of motion is given by:

$$
M(q)\ddot{q}+C(q,\dot{q})\dot{q}+G(q)=Bu + \sum J_c^T\lambda
$$

where,

- $q$ is the robot coordinates, including the floating base [pos, euler, $q_b$], its dimension is $n_q$

- $q_b$ is the relative joint angles in the robotic system, its dimension is $n_b$

- $u$ are the joint torques, note that sometimes we include passive joints or fixed joints that can be dealt adding constraints.

- B maps the torques to the joint accelerations, to simplify the problem definition we assume that $n_q = n_b+6$ and:
  
  $$
  B = \begin{bmatrix} 0_{6,6} & 0_{6,n_b} \\ 0_{n_b, 0} & I_{nb,nb} \end{bmatrix}
  $$

- $J_c$ is the contact jacobian for each contact point

- $\lambda$ is the contact wrench for each contact point

And the control objectives (Assumed relative degree two) are specified as:

$$
y=h^a(t)-h^d(t,q,\dot{q})
$$
$$
\dot{y}=J_y \dot{q}-\dot{h}^d
$$
$$
\ddot{y} = \dot{J}_y \dot{q}+ J_y \ddot{q} - \ddot{h}^d
$$

In both QP controllers, we are solving a similar opimization problem:

$$
\min_{\mathcal{X}} C(\mathcal{X}) + W(\mathcal{X}) \\
s.t \hspace{1em} A_{ineq} \mathcal{X} \leq B_{ineq} \mathcal{X} \\
s.t \hspace{1em} A_{eq} \mathcal{X} \leq B_{eq} \mathcal{X} \\
$$

where,

$ \mathcal{X} = [\ddot{q}, u, \lambda]^T$

and,

$W(\mathcal{X}) = \mathcal{X}^T E \mathcal{X}$ is a regularization term with $E=\sigma I$

For **Task space PD control** we are using:

$$
C(\mathcal{X}) = || \dot{J}_y \dot{q}+J_y \ddot{q}-y*||^2\\
y*=-Kpy-K_d\dot{y}
$$

For **CLFQP**

$$
C(\mathcal{X})=0
$$

Note theta we include the CLFQP stability constraint in the inequality constraints. (Link to paper)

$$
L_gL_fy*u \leq -L_f^2y-cV(\eta) \\
V(\eta) = \eta^T P_{\epsilon} \eta
$$

## Controller declaration

The selection of the controller is specified when creating the IDQP class.

```c
  controller = new IDQPcpp(nQ, nObj, nMotors, nContacts);
```

where,

- nQ: is the dimension of $q$, in other words $nQ = n_q$ 

- nObj: is the number of ouput objectives, given by the dimension of the vector $y$ 

- nMotors: is the number of motors or actuated joints.

- nContacts: is the number of contact points

Note that for every contact point we must provide the corresponding jacobians and its time-derivative. 

## Adding Constraints

The constraints can be added sequentially with function calls, IDQP provides function members for typical constraints as well as custom constraints.

### Dynamic Constraints

### + Function Based Dynamics

```c
controller.addDynamicsHvec(MassMatrix, Hvector);
```

The ``` MatrixXd MassMatrix(VectorXd &) ``` is a function pointer that receives $q$ and computes the matrix $M(q)$. 
The ``` VectorXd Hvector(VectorXd &, VectorXd &)` is a function pointer that computes Coriolis and Gravity as a single vector $H(q, \dot{q}) = C(q, \dot{q}) \dot{q} + G(q)$ 

The function pointer takes care of the update step whenever you call the method to solve the QP.

### + Numerical Dynamics

If a symbolic expression for the dynamics isn't available, we can use any numerical solution externally (by using RBDL, Pinocchio or hand calculated for instance)

```c
controller.addDynamicsHvecNum(MassMatrix, Hvector);
```

Here, ``` MassMatrix``` is a variable of type ``` MatrixXd ``` and ``` Hvector ``` is a vector of type ```VectorXd```. When you use this function, you must take responsibility of the update event every control step with the function: 

```c
controller.updateDynamicsHvecNum(MassMatrix, Hvector);
```

### Contact Constraints

Each contact constraint is specified by a contact jacobian ($J_c$) and its derivative ($\dot{J}_c$). If we have several contacts, we specify each pair as part of a list. The List types are provided by the class.

```c
  IDQP_type::List_Jac List_Jac_contact;
  IDQP_type::List_Jac_dot List_Jac_contact_dot;
```

Then, we populate each list with pushback operators:

```c
  List_Jac_contact.push_back(rExpr.JLeft);
  List_Jac_contact_dot.push_back(rExpr.DJLeft);
```

The Jacobians $J(q)$  are function pointers of the form ``` MatrixXd J(VectorXd &)```, while the time-derivative of the jacobians $\dot{J}(q,\dot{q})$ have the form ```MatrixXd DJ(VectorXd &, VectorXd &)```. Finally, we add the lists to specify the robot's contact constraints.

```c
controller.addContactJacobians(List_Jac_contact, List_Jac_contact_dot);
```

### Friction Constraints

Specify the friction coefficient with the variable ```mu``` of type ```double``` 

```c
controller.addFrictionCone(mu);
```

### Control Limits

Specify the limits of the control action, where ```uMin``` and ```uMax``` are ```VectorXd``` constants.

```c
controller.addControlLimits(uMin, uMax);
```

### Custom Constraints

### + Simple Constant Constraints

This can be manually constructed and represents:

$$
A^c_{eq} \mathcal{X} \leq B^c_{eq}
$$

```c
controller.addCustomEq(Aeq_DS, Beq_DS);
```

### + Type I constraints

These are referred to symbolic constraints of the form: ($c(q,\dot{q})=0$)

$$
J_1(q,\dot{q}) \ddot{q}+J_2(q,\dot{q}) \dot{q} = 0
$$

```c
controller.addTypeIconstraint(List_Jac, List_Jac_dot, nRel);
```

where, ```nRel``` is the number of rows of each jacobian. Note that, (for now) all the jacobians on the list must have the same number of rows. 

### + Type II constraints

These are referred to symbolic constraints of the form: ($c(q)=0$)

$$
J(q) \ddot{q} + \dot{J}(q, \dot{q})\dot{q}=0
$$

```c
controller.addTypeIIconstraint(List_JacII, List_JacII_dot, nRel);
```

where, `nRel` is the number of rows of each jacobian. Note that, (for now) all the jacobians on the list must have the same number of rows.

## Adding the Control Objectives

The control objectives specify the control ouputs to be solved. ```Kp``` and ```Kd``` are used only for Task Space PD control. The other parameters are function pointers of the form:

- ```c
  VectorXd y(VectorXd &q, VectorXd &dq)
  ```

- ```c
  MatrixXd Jy(VectorXd &q, VectorXd &dq)
  ```

- ```c
  MatrixXd DJy(VectorXd &q, VectorXd &dq)
  ```

- ```c
  VectorXd yd(VectorXd &q, VectorXd &dq, VectorXd &var, VectorXd &dvar)
  ```

- ```c
  VectorXd dyd(VectorXd &q, VectorXd &dq, VectorXd &var, VectorXd &dvar)
  ```

- ```c
  (optional) VectorXd ddyd(VectorXd &q, VectorXd &dq, VectorXd &var, VectorXd &dvar)
  ```

```c
controller.addObjectives(y, Jy, DJy, yd, dyd, Kp, Kd);
```

## Adding Regularization

$$
W(\mathcal{X}) = \ddot{q}^T I_q r_q \ddot{q}+u^TI_ur_uu + \lambda^TI_{\lambda}r_{\lambda}\lambda
$$

```c
controller.addRegularizationTerm(RegQ, RegU, RegL);
```

## Build Optimization

```c
controller.buildOptimization();
```

## Solve Controller

```c
controller.solve(t, q, dq, var, dvar);
```
