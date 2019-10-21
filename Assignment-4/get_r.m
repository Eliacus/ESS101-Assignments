function r = get_r(butcher,x,dt,f, k)

[s, ~] = size(butcher.a);

r = sym('r',[s,length(x)],'real');

for i=1:s
    asum = zeros(1,length(x));
    for j=1:s
        asum = asum + butcher.a(i,j).*k(j,:);
    end
    r(i,:) = (f(x + dt.*asum)' - (k(i,:)));
end
end