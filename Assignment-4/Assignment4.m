%% 1a

% Butcher 
a = [1/4 (1/4)-sqrt(3)/6; 
     (1/4)+sqrt(3)/6 1/4];
b = [.5 .5];
c = [.5-sqrt(3)/6 .5+sqrt(3)/6];
butcher = struct('a',a,'b',b,'c',c);

% Sizes
[s, ~] = size(butcher.a);
n = length(x0);

% Tolerance 
tol = 1e-3;

% Initial values
x0 = [1 1];
xk = zeros(s,n);

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

% r
r = get_r(butcher,x0,dt,f,K);

k = zeros(s,n);
for n=1:length(x0)
    k(:,n) = [1; 1];
    for i=1:length(K)
        nor = 1;
        while nor > tol
            dr = jacobian(r(:,n),K(:,n));
            dk = -dr\r(:,n);
            dk = subs(dk,K(:,n),k(:,n));
            
            k(:,n) = k(:,n) + alpha*dk;
            
            nor = norm(double(get_r(butcher,xk(i,n),dt,f,k(:,n))));
        end
    end
end


