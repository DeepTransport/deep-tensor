function fx = eval_potential_reference(obj, x)
% Evaluate the normalised (marginal) pdf represented by squared FTT.
%   f = EVAL_POTENTIAL_REFERENCE(irt, x)
%
%   x - input variables in the reference coordinates
%   f - potential function at x

% map to the reference coordinate

d = ndims(obj.approx);
dz = size(x,1);
if obj.int_dir > 0
    % marginalised from the right
    fxl = eval_block(obj.approx, x, obj.int_dir);
    if dz < d
        fx = sum((fxl*obj.ms{dz+1}).^2, 2)';
    else
        fx = fxl'.^2;
    end
    mlogw = eval_measure_potential_reference(obj.approx.base, x, 1:dz);
    fx = log(obj.z) - log(fx+obj.tau) + mlogw;
else
    % marginalised from the left
    fxg = eval_block(obj.approx, x, obj.int_dir);
    if dz < d
        ie = (d-dz)+1;
        fx = sum((obj.ms{ie-1}*fxg).^2, 1);
    else
        fx = fxg.^2;
    end
    ie  = (d-dz)+1;
    mlogw = eval_measure_potential_reference(obj.approx.base, x, d:-1:ie);
    fx = log(obj.z) - log(fx+obj.tau) + mlogw;
end

end

