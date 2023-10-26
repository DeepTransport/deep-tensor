function r = ttrankest(bases, N)
% Solve the quadratic to estimate the TT rank matching N evaluations
% (n(1)+n(d))*r + sum(n(2:d-1))*r^2 = N
% ar^2 + br - N = 0          a = sum(n(2:d-1)), b = n(1)+n(d)
% D = b^2 + 4*a*N
% r = (-b + sqrt(D))/(2*a)

b = numel(bases.oneds{1}.nodes) + numel(bases.oneds{end}.nodes);
a = sum(cellfun(@(b)numel(b.nodes), bases.oneds(2:end-1)));
if a==0
    % Linear equation in 2D
    r = N/b;
else
    r = (-b + sqrt(b^2 + 4*a*N))/(2*a);
end
r = round(r);
end
