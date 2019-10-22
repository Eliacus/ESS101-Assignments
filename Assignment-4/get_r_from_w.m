function r = get_r_from_w(butcher,x,dt,f, k, z)

[s, ~] = size(butcher.a);
m = length(z);
n = length(x);

r = sym('r',[s,n],'real');

for i=1:s
    asum = zeros(1,length(x));
    for j=1:s
        asum = asum + butcher.a(i,j).*k(j,:);
    end
    r(i,:) = f(x + dt.*asum',z(i))'-k(i,:);
end
    r = reshape(r,12,1);
 
end