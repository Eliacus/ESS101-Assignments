function [r, dr] = get_r(butcher,xk,dt,f, k)

[s, ~] = size(butcher.a);

r = sym('r',[s,length(xk)],'real');

for i=1:s
    asum = zeros(1,length(xk));
    for j=1:s
        asum = asum + butcher.a(i,j).*k(j,:);
    end
    r(i,:) = (f(xk + dt.*asum) - (k(i,:)));
end
end