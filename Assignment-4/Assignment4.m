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
tol = 1e-2;

% Constants
dt = 1e-2;
tf = 25;
N = tf/dt;
alpha = 0.1;

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
