function f = eval_potential(obj, x, k)
% Evaluate the deep Rosenblatt transport Z = T(X), where Z is the reference
% random variable and X is the target random variable.
%   logf = LOG_PDF(dirt, x, k)
%
%   x    - random variable drawn form the pdf defined by DIRT
%   k    - number of layers in the evaluations
%   logf - log of the DIRT density

if nargin <= 2
    k = num_layers(obj);
else
    k = min(k, num_layers(obj));
end
[~,f] = eval_rt(obj, x, k);

end
