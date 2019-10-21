%% Question 1a

syms theta phi L phi_dot theta_dot m2 g m1 z real

u = sym('u',[3,1],'real');

p1 = sym('p1',[3,1],'real');
dp1 = sym('dp1',[3,1],'real');

q = [p1; theta; phi];
dq = [dp1; theta_dot; phi_dot];

p2 = [p1(1) + cos(theta)*sin(phi)*L;
      p1(2) + sin(theta)*sin(phi)*L;
       p1(3)-L*cos(phi)]; 
dp2 = jacobian(p2,q)*dq;

T_1 = (1/2)*m1*(dp1)'*dp1; 
V_1 = m1*g*[0 0 1]*p1;

T_2 = (1/2)*m2*(dp2)'*dp2; 
V_2 = m2*g*[0 0 1]*p2;

L = (T_1+T_2)-(V_1+V_2);

F = [u; 0; 0];

L_gradient_q_dot = (jacobian(L,dq))';
W_qdot_jacobian = jacobian(L_gradient_q_dot,q);

L_gradient_q  = (jacobian(L,q))';
L_gradient_q_dot_dt = F+L_gradient_q;

Mv = L_gradient_q_dot_dt - W_qdot_jacobian*dq;

simplify(Mv)


%% Question 1b 
syms theta phi L phi_dot theta_dot m2 g m1 z real

u = sym('u',[3,1],'real');

p1 = sym('p1',[3,1],'real');
p2 = sym('p2',[3,1],'real');

dp1 = sym('dp1',[3,1],'real');
dp2 = sym('dp2',[3,1],'real');

ddp1 = sym('ddp1',[3,1],'real');
ddp2 = sym('ddp2',[3,1],'real');

q = [p1;p2];
dq = [dp1; dp2];
ddq = [ddp1; ddp2];

e = p1 - p2;
c = (1/2)*(e'*e-L^2);

T_1 = (1/2)*m1*(dp1)'*dp1; 
V_1 = m1*g*[0 0 1]*p1;

T_2 = (1/2)*m2*(dp2)'*dp2; 
V_2 = m2*g*[0 0 1]*p2;

L = (T_1+T_2)-(V_1+V_2) - z'*c;

F = [u; 0; 0; 0];
gradient_L_qdot = (jacobian(L,dq))';
jacobian_W_qdot = jacobian(gradient_L_qdot,q)*dq;

L_gradient_q  = (jacobian(L,q))';
L_gradient_q_dot_dt = F+L_gradient_q;

Mv = L_gradient_q_dot_dt - jacobian_W_qdot;
