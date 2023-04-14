function f = banana(X, sigma, a)

x = X(1,:);
y = X(2,:);

f = (y-a*x.^2+1).^2*10 + x.^2;
f = f/sigma;

end
