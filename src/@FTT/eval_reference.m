function fx = eval_reference(obj, x)
% Evaluate the FTT for either the first or last k variables.
%   f = EVAL_REFERENCE(tt, x)
%
%   x   - input variables, k x n, in the reference domain
%   dir - direction evlauation, >0: from left to right
%                               <0: from right to left
%   f   - function values at x, m x n

fx  = eval_block(obj, x, obj.direction);
if obj.direction > 0
    fx  = fx';
end

end