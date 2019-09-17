%% Question 1a

syms theta phi L phi_dot theta_dot m2 g m1 u real

u = sym('u',[3,1],'real');

p_1 = sym('p_1',[3,1],'real');
p_1_dot = sym('p_1_dot',[3,1],'real');

q = [p_1; theta; phi];
q_dot = [p_1_dot; theta_dot; phi_dot];

p_2 = [p_1(1) + cos(theta)*sin(phi)*L;
      p_1(2) + sin(theta)*sin(phi)*L;
       p_1(3)-L*cos(phi)]; 
p_2_dot = jacobian(p_2,q)*q_dot;

T_1 = (1/2)*m1*(p_1_dot)'*p_1_dot; 
V_1 = m1*g*[0 0 1]*p_1;

T_2 = (1/2)*m2*(p_2_dot)'*p_2_dot; 
V_2 = m2*g*[0 0 1]*p_2;

L = (T_1+T_2)-(V_1+V_2);

L_gradient_q_dot = (jacobian(L,q_dot))';
W_jacobian = jacobian(L_gradient_q_dot,q);

F = [u; 0; 0];

L_gradient_q  = (jacobian(L,q))';
L_gradient_q_dot_dt = F+L_gradient_q;

Mv = L_gradient_q_dot_dt - W_jacobian*q_dot;

simplify(Mv)


%% Question 1b 

q_new = [p_1; p_2];
q_new_dot = [p_1_dot; p_2_dot];
e = p_1 - p_2;

C = (1/2)*(e'*e-L^2);
C = simplify(C);

