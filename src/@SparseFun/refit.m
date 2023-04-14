function obj = refit(obj,y)

[coeff,err] = SparseFun.ls_solve(obj.A,y',obj.w);
        
fprintf('Refit CV err: %4.4e |\n',norm(err));

obj.data = coeff;
obj.err = err;
obj.y = y;

end
