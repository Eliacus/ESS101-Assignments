function r = get_r2_from_w(butcher,x,dt,f, k, z)

[s, ~] = size(butcher.a);
m = length(z);
n = length(x);

r = sym('r',[m,1],'real');

for i=1:s
    asum = zeros(1,length(x));
    for j=1:s
        asum = asum + butcher.a(i,j).*k(j,:);
    end
    r(i,:) = f(x + dt.*asum')';
end
    %r = reshape(r,12,1);
 
end