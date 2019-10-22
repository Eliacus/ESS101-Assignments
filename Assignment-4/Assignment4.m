%% 1a

% Butcher 
a = [1/4 (1/4)-sqrt(3)/6; 
     (1/4)+sqrt(3)/6 1/4];
b = [.5 .5];
c = [.5-sqrt(3)/6 .5+sqrt(3)/6];
butcher = struct('a',a,'b',b,'c',c);

% Initial values
x0 = 1;

% Sizes
[s, ~] = size(butcher.a);
n = length(x0);

% Tolerance 
tol = 1e-2;

% Constants
lambda = -2;
dt = 0.1;
tf = 2;
N = tf/dt;
alpha = 0.1;

% Function 
f = @(x) lambda * x;

% K

K = sym('k',[s,n],'real');
syms x
% r
r = get_r(butcher,x,dt,f,K);
r_function=matlabFunction(r,'Vars',{'k1','k2','x'});

t = 0:dt:tf;

% Initialize matrices
xk = zeros(N,n);
xk(1,:) = x0;
k = zeros(s,n);

for i=2:N+1
    i
    for n=1:length(x0)
        k(:,n) = [0; 0];
        nor = 1;
        while nor > tol
            dr = jacobian(r(:,n),K(:,n));
            dk = -dr\r(:,n);
            dk = subs(dk,K(:,n),k(:,n));
            dk = subs(dk,x,xk(i-1,n));
            
            k(:,n) = k(:,n) + alpha*dk;
            
            nor = norm(double(r_function(k(1,n),k(2,n),xk(i-1,n))));
        end
    end
        xk(i,:) = xk(i-1,:) + dt* butcher.b*k;
end

plot(t,xk)
hold on
%% 2c


% Butcher 
a = [1/4 (1/4)-sqrt(3)/6; 
     (1/4)+sqrt(3)/6 1/4];
b = [.5 .5];
c = [.5-sqrt(3)/6 .5+sqrt(3)/6];
butcher = struct('a',a,'b',b,'c',c);


% Initial values
x0 = [1 0];

% Sizes
[s, ~] = size(butcher.a);
n = length(x0);

% Function
f = @(x) [x(2);
          5*(1-x(1)^2)*x(2)-x(1)];

% Tolerance 
tol = 1e-4;

% Constants
dt = 1e-2;
tf = 25;
N = tf/dt;
alpha = 0.01;

% K

K = sym('k',[s,n],'real');
x = sym('x',[1,2],'real');

% r
r = get_r(butcher,x,dt,f,K);
r_function=matlabFunction(r,'Vars',{'k1_1','k2_1','k1_2','k2_2','x1','x2'});
r_r = reshape(r,s*n,1);
K_r = reshape(K,1,s*n);

dr = jacobian(r_r,K_r);

% Initialize matrices 
xk = zeros(N,n);
xk(1,:) = x0;
k = zeros(s,n);
t = 0:dt:tf;

dk = -dr\r_r;
dk_fun = matlabFunction(dk,'Vars',{'k1_1','k2_1','k1_2','k2_2','x1','x2'});

for i=2:N+1
    k = zeros(4,1);
    nor = 1;
    while nor > tol
        dkk = dk_fun(k(1),k(2),k(3),k(4),xk(i-1,1),xk(i-1,2));

        k = k + alpha*dkk;

        nor = norm(double(r_function(k(1),k(2),k(3),k(4),xk(i-1,1),xk(i-1,2))));
    end
    xk(i,:) = xk(i-1,:) + dt* butcher.b*reshape(k,2,2);
end

plot(t,xk(:,1))
hold on
plot(t,xk(:,2))

%% 3a


% Butcher 
a = [1/4 (1/4)-sqrt(3)/6; 
     (1/4)+sqrt(3)/6 1/4];
b = [.5 .5];
c = [.5-sqrt(3)/6 .5+sqrt(3)/6];
butcher = struct('a',a,'b',b,'c',c);

% Create variables
p = sym('p',[3,1],'real');
v = sym('v',[3,1],'real');
%dp = sym('dp',[3,1],'real');
%dv = sym('dv',[3,1],'real');
syms z real 

x = [p; v];
x0 = [1; 0; 0; 0; 1; 0];

% Dimension variables
n = length(x0);
s = 2;

% Initialize constants
m = 1;
g = 9.82;
L = 1;

% Simulation parameters
tol = 1e-4;
tf = 2;
dt = 0.1;
N = tf/dt;
t = 0:dt:tf;

% Equations 
dv = -g*[0;0;1]-((z*p)/m);
gf = (p'*dv)+(v'*v);
f = [v; dv];
fc=0.5*(((p')*p)-L^2);

% Function generations
f_func = matlabFunction(f,'Vars',{x,'z'});
fc_func = matlabFunction(gf,'Vars',{x,'z'});

K = sym('k',[s,n],'real');
Z = sym('z',[s,1],'real');

r_f = get_r_from_w(butcher,x,dt,f_func,K,Z);

r_g = get_r2_from_w_a(butcher,x,dt,fc_func,K,Z);

r = [r_f;r_g];

K_r = reshape(K,[2*n,1]);

w = [K_r;Z];

dr = jacobian(r,w);

r_func = matlabFunction(r,'Vars',{w,x});

dr_func = matlabFunction(dr,'Vars',{w,x});

xk = zeros(n, N);
xk(:,1) = x0';
alpha=0.1;
Cq=zeros(1,N);
Cq(1)=0.5*((([1;0;0]')*[1;0;0])-L^2);

for i=2:N+1
    W = zeros(14,1);
    nor = 1;
    while nor > tol
        dW = -dr_func(W,xk(:,i-1))\r_func(W,xk(:,i-1)); 

        W = W + alpha*dW;
    
        nor = norm(double(r_func(W,xk(:,i-1))));
    end
    xk(:,i) = xk(:,i-1) + dt* (butcher.b*reshape(W(1:12),2,6))';
    Cq(i)=0.5*(((xk(1:3,i)')*xk(1:3,i))-L.^2);
end

%Plott them here: 
figure(1)
plot(t,Cq); title('C(q)');

figure(2)
plot(t,xk); title('x');
legend('p1','p2', 'p3', 'v1','v2','v3');


%% 3b

% Butcher 
a = [1/4 (1/4)-sqrt(3)/6; 
     (1/4)+sqrt(3)/6 1/4];
b = [.5 .5];
c = [.5-sqrt(3)/6 .5+sqrt(3)/6];
butcher = struct('a',a,'b',b,'c',c);

% Create variables
p = sym('p',[3,1],'real');
v = sym('v',[3,1],'real');
%dp = sym('dp',[3,1],'real');
%dv = sym('dv',[3,1],'real');
syms z real 

x = [p; v];
x0 = [1; 0; 0; 0; 1; 0];

% Dimension variables
n = length(x0);
s = 2;

% Initialize constants
m = 1;
g = 9.82;
L = 1;

% Simulation parameters
tol = 1e-4;
tf = 2;
dt = 0.1;
N = tf/dt;
t = 0:dt:tf;

% Equations 
dv = -g*[0;0;1]-((z*p)/m);
gf = (p'*dv)+(v'*v);
f = [v; dv];
fc=0.5*(((p')*p)-L^2);

% Function generations
f_func = matlabFunction(f,'Vars',{x,'z'});
fc_func = matlabFunction(fc,'Vars',{x});

K = sym('k',[s,n],'real');
Z = sym('z',[s,1],'real');

r_f = get_r_from_w(butcher,x,dt,f_func,K,Z);

r_g = get_r2_from_w(butcher,x,dt,fc_func,K,Z);

r = [r_f;r_g];

K_r = reshape(K,[2*n,1]);

w = [K_r;Z];

dr = jacobian(r,w);

r_func = matlabFunction(r,'Vars',{w,x});

dr_func = matlabFunction(dr,'Vars',{w,x});

xk = zeros(n, N);
xk(:,1) = x0';
alpha=0.1;
Cq=zeros(1,N);
Cq(1)=0.5*((([1;0;0]')*[1;0;0])-L^2);

for i=2:N+1
    W = zeros(14,1);
    nor = 1;
    while nor > tol
        dW = -dr_func(W,xk(:,i-1))\r_func(W,xk(:,i-1));

        W = W + alpha*dW;
    
        nor = norm(double(r_func(W,xk(:,i-1))));
    end
    xk(:,i) = xk(:,i-1) + dt* (butcher.b*reshape(W(1:12),2,6))';
    Cq(i)=0.5*(((xk(1:3,i)')*xk(1:3,i))-L.^2);
end

%Plots
figure(1)
plot(t,Cq,'k'); title('C(q)');

figure(2)
plot(t,xk); title('x');
legend('p1','p2', 'p3', 'v1','v2','v3');


