
function [mllkd,mlp] = log_target_pullback_pcn(irt,func,ry,z)

mlf = pullback(irt, func, [repmat(ry,1,size(z,2));z]);
%
mlp = 0.5*sum(z.^2,1);
mllkd = mlf - mlp;

end