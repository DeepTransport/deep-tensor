function [f1, f2, g1, g2] = joint_banana(z, sigma)

y = z(1,:);
u = z(2:3,:);
F = log( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
G = [F-3; F-5];
%
mllkd = sum((G-y).^2,1)/(2*sigma^2);
mlp = 0.5*sum(u.^2,1);

f1 = mllkd + mlp - 0.5*sum(z.^2,1);
f2 = 0.5*sum(z.^2,1);

if nargout > 2
    tmp = [2*(u(1,:)-1)+400*(u(1,:).^2-u(2,:)).*u(1,:); 200*(u(2,:)-u(1,:).^2)]...
        ./( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
    g1 = [-(2*F-2*y-8); (2*F-2*y-8)*tmp]/sigma^2;
    g2 = z;
end

end