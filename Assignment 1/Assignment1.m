%% Question 1a
clc
% Defining variables
syms theta phi L phi_dot theta_dot m2 g m1 u real
u = sym('u',[3,1],'real');
p1 = sym('p1_',[3,1],'real');
dp1 = sym('dp1_',[3,1],'real');

% Setting up q vector and q_dot
q = [p1; theta; phi];
dq = [dp1; theta_dot; phi_dot];

% Setting up p_2 and p_2_dot
p2 = [p1(1) + cos(theta)*sin(phi)*L;
      p1(2) + sin(theta)*sin(phi)*L;
       p1(3)-L*cos(phi)]; 
dp2 = jacobian(p2,q)*dq;

% Kinetic energy
T_1 = (1/2)*m1*(dp1)'*dp1; 
T_2 = (1/2)*m2*(dp2)'*dp2; 
T = T_1 + T_2;

% Potential energy
V_1 = m1*g*[0 0 1]*p1;
V_2 = m2*g*[0 0 1]*p2;
V = V_1 + V_2;

% Lagrangian
L = T-V;

% External force vector
F = [u; 0; 0];

% Gradient of L wrt dq  ( = W(q)dq )
L_gradient_q_dot = (jacobian(L,dq))';

% Gradient of L wrt q
L_gradient_q  = (jacobian(L,q))';

% Since d/dt(dL/dqdot) - dL/dq = F, we can compute:   Eq. 2.112
L_gradient_q_dot_dt = F+L_gradient_q;

% Calculating d/dt(grad_qdot(L))  Eq. 2.107b last term
W_qdot_jacobian_q_dot = jacobian(L_gradient_q_dot,q)*dq;

% Extracting W(q)*ddq from eq 2.107b
Mv = simplify(L_gradient_q_dot_dt - W_qdot_jacobian_q_dot);



%% Question 1b 
clc
% Defining variables
syms theta phi L phi_dot theta_dot m2 g m1 u z real

u = sym('u',[3,1],'real');

p1 = sym('p1_',[3,1],'real');
p2 = sym('p2_',[3,1],'real');

dp1 = sym('dp1_',[3,1],'real');
dp2 = sym('dp2_',[3,1],'real');

ddp1 = sym('ddp1_',[3,1],'real');
ddp2 = sym('ddp2_',[3,1],'real');

% q and its derivates vectors
q = [p1;p2];
dq = [dp1; dp2];
ddq = [ddp1;ddp2];

% Calculating the constraint
e = p1 - p2;
c = (1/2)*(e'*e-L^2);

% Kinetic energy
T_1 = (1/2)*m1*(dp1)'*dp1; 
T_2 = (1/2)*m2*(dp2)'*dp2; 
T_tot = T_1 + T_2;

% Potential energy
V_1 = m1*g*[0 0 1]*p1;
V_2 = m2*g*[0 0 1]*p2;
V_tot = V_1+V_2;

% Calculating Lagrangian
L = T_tot - V_tot - z'*c;

% External force vector
F = [u; 0; 0; 0];

% dL/dq
grad_q_L = jacobian(L,q)';

% dL/ddq 
grad_dq_L = jacobian(L,dq)';

dq_Wdq_dq = jacobian(grad_dq_L,q)*dq;

grad_dq_dt_L = F + grad_q_L;

% Finding Mv with eq. 2.167a. 
W_ddq = grad_dq_dt_L - dq_Wdq_dq;


%% Question 2a

% Q = F here. 
grad_dq_dt_L = jacobian(p1,q)'*u + grad_q_L;

% Finding Mv with eq. 2.167a. 
W_ddq = grad_dq_dt_L - dq_Wdq_dq;

% Computing the W matrix
W = jacobian(grad_dq_L,dq);

a = jacobian(c,q)';

M = [W, a;
    a' 0];

c_qdqu =  [F-dq_Wdq_dq+jacobian((T_tot),q)'-jacobian((V_tot),q)';
           -jacobian(jacobian(c,q)*dq,q)*dq];


% rank(M) = 7

sol = simplify(inv(M)*c_qdqu);

% We can see that the explicit form is way more complex than the implicit
% form

