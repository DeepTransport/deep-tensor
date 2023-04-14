function [mlf,gmlf] = log_target_pullback_nuts(irt,func,ry,z)

[mlf,gmlf] = pullback(irt, func, [repmat(ry,1,size(z,2));z]);
%
ind = length(ry) + (1:size(z,1));
gmlf = gmlf(ind,:);

end