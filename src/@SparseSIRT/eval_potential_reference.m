function f = eval_potential_reference(obj, x)
% Evaluate the normalised (marginal) pdf represented by squared FTT.
%   f = EVAL_PDF(irt, x)
%
%   x - input variables
%   f - marginal density at x

d = ndims(obj.approx);
nk = cardinal(obj.approx);
[dx,n] = size(x);
brl = repmat(obj.approx.data.coeff(:)', n, 1); % n x nk
for k = 1:dx
    bk = eval_basis(obj.approx.base.oneds{obj.order(k)}, reshape(x(obj.order(k),:),[],1));
    brl = brl.*bk(:,obj.approx.data.I.array(:,obj.order(k))+1);
end
if dx < d
    f = sum((reshape(brl,n,nk)*obj.cum_P{obj.order(dx)}).^2,2);
else
    f = sum(reshape(brl,n,nk),2).^2;
end
mlogw = eval_measure_potential_reference(obj.approx.base, x, 1:dx);
f = log(obj.z) - log(f(:)'+obj.tau) + mlogw;

end
