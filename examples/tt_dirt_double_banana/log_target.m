function [mlf,gmlf] = log_target(func,y,x)

[mllkd,mlp,gmllkd,gmlp] = func([repmat(y,1,size(x,2));x]);
%
mlf = mllkd + mlp;
ind = length(y) + (1:size(x,1));
gmlf = gmllkd(ind,:)+gmlp(ind,:);

end