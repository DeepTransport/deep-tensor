function x = quadratic_root(a,b,c,left,right)

x = zeros(size(c));
tmp1 = zeros(size(c));
tmp2 = zeros(size(c));
%
indp = (b>=0);
indm = ~indp;
%
tmp = b.^2-4.*a.*c;
if norm(tmp(tmp<0)) > 1E-6
    warning('root finding does not have a real root')
end
tmp(tmp<0) = 0;
tmp = sqrt(tmp);
%
tmp1(indp) = (-b(indp)-tmp(indp))./(2*a(indp));
tmp2(indp) = (2*c(indp))./(-b(indp)-tmp(indp));
%
tmp1(indm) = (2*c(indm))./(-b(indm)-tmp(indm));
tmp2(indm) = (-b(indm)-tmp(indm))./(2*a(indm));
%
% sort
x1 = min(tmp1, tmp2);
x2 = max(tmp1, tmp2);
%
b_ind_x1l = x1<=left;
b_ind_x2r = x2>=right;
%
ind_lot_rot = b_ind_x1l  &  b_ind_x2r; % both roots are not in the interval
ind_lot_rin = b_ind_x1l  & ~b_ind_x2r; % left root out, right root in
ind_lin_rot = ~b_ind_x1l &  b_ind_x2r; % left root in,  right root out
ind_lin_rin = ~b_ind_x1l & ~b_ind_x2r; % both roots are in the interval
%
if sum(ind_lin_rin) > 0
    x(ind_lin_rin) = x1(ind_lin_rin);
end
%
if sum(ind_lot_rot) > 0
    d1 = abs(x1(ind_lot_rot)-left(ind_lot_rot));
    d2 = abs(x2(ind_lot_rot)-right(ind_lot_rot));
    if sum(d1<d2) > 0
        x(ind_lot_rot(d1<d2)) = x1(ind_lot_rot(d1<d2));
    end
    if sum(d1>=d2) > 0
        x(ind_lot_rot(d1>=d2)) = x2(ind_lot_rot(d1>=d2));
    end
end
%
if sum(ind_lot_rin) > 0 
    x(ind_lot_rin) = x2(ind_lot_rin);
end
%
if sum(ind_lin_rot) > 0 
    x(ind_lin_rot) = x1(ind_lin_rot);
end

end