close all
lambda = -2;

dts = [1 0.5 0.1 0.01 0.001];
global_error = zeros(4,3);

for tt=1:length(dts)
    
dt = dts(tt);

tf = 2;

x = zeros(tf/dt,1);

x(1) = 1;

RK_orders = [1 2 4];


for k = 1:length(RK_orders)
    
RK_order = RK_orders(k);

switch RK_order
    case 1
        a = [0];
        b = [1];
        c = [0];
        butcher = struct('a',a,'b',b,'c',c);
    case 2
        a = [0 0; .5 0];
        b = [0 1];
        c = [0 .5];
        butcher = struct('a',a,'b',b,'c',c);
    case 4
        a = [0 0 0 0; .5 0 0 0; 0 .5 0 0; 0 0 1 0];
        b = [1/6 1/3 1/3 1/6];
        c = [0 .5 .5 1];
        butcher = struct('a',a,'b',b,'c',c);
end
f = @(x) lambda * x;

u = @(t) 0;

for i=2:tf/dt+1
   x(i) = generic_RK(butcher,x(i-1),dt,f,i*dt,u);
end
global_error(tt,k) = abs(x(end)-xsol(tf));

end
end

% b)

% Plotting

loglog(dts,global_error(:,1),dts,global_error(:,2),dts,global_error(:,3))
xlabel('dt')
ylabel('global error')
legend('RK1','RK2','RK4')
