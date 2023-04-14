function [f1, f2] = banana_dirt(X, sigma, a)

x = X(1,:);
y = X(2,:);

f1 = (y-a*x.^2+2).^2*10;
f1 = f1/sigma;

f2 = (y.^2+x.^2)/2;

end
