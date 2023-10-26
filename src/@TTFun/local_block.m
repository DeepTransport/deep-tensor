function [f, f_evals] = local_block(oned, xleft, xright, func)
% Evaluate the user function at cross indices and qudrature or
% intepolation points weighed by the interplation matrix
%
nnodes = length(oned.nodes);
if isempty(xleft)
    % left boudnary
    nileft  = 1;
    niright = size(xright, 2);
    % f       = zeros(niright, oned_def.nquad);
    % define parameters, each xright binding with a block of xquad
    param   = [ repmat(oned.nodes(:)', 1, niright); ...
        reshape(repmat(xright, nnodes, 1), size(xright,1), niright*nnodes) ];
    % return of func is a vector, reshape to a matrix
    % 1st index: xquad, 2nd: xright
    %f = feval_reference(obj,func,param);
    f = func(param);
    f = reshape(f', 1, nnodes, niright, []);
elseif isempty(xright)
    % right boundary
    nileft  = size(xleft, 2);
    niright = 1;
    % f       = zeros(nileft, oned_def.nquad);
    % define parameters, each xquad binding with a block of xleft
    param   = [ repmat(xleft, 1, nnodes); ...
        reshape( repmat(oned.nodes(:)', nileft, 1), 1, nileft*nnodes ) ];
    % return of func is a vector, reshape to a matrix
    % 1st index: xleft, 2nd: xquad
    %f = feval_reference(obj,func,param);
    f = func(param);
    f = reshape(f', nileft, nnodes, 1, []);
else
    %
    nileft  = size(xleft, 2);
    niright = size(xright, 2);
    % f       = zeros(nileft, niright, oned_def.nquad);
    
    tmp_lq  = [ repmat(xleft, 1, nnodes); ...
        reshape( repmat(oned.nodes(:)', nileft, 1), 1, nileft*nnodes ) ];
    
    param   = [ repmat(tmp_lq, 1, niright); ...
        reshape(repmat(xright, nnodes*nileft, 1), size(xright,1), nileft*nnodes*niright) ];
    
    % return of func is a vector, reshape to tensor
    %f = feval_reference(obj,func,param);
    f = func(param);
    f = reshape(f', nileft, nnodes, niright, []);
end

%mesh(reshape(f, nileft*nnodes, []));

if isa(oned, 'Spectral')
    f = oned.node2basis*reshape(permute(f, [2,1,3,4]), nnodes, []);
    f = permute(reshape(f, cardinal(oned), nileft, niright, []), [2 1 3 4]);
end
f_evals = size(param, 2);

end