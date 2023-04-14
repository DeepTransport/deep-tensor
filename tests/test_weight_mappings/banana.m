function [f,g] = banana(X, sigma, a)

x = X(1,:);
y = X(2,:);

tmp = (y-a*x.^2+1);
f = tmp.^2*10;
f = f/sigma + 0.5*x.^2;

g = x + (20/sigma)*tmp.*[-2*a*x ; ones(1,size(x,2))];

end
