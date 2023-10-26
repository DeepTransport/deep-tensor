classdef SparseData
    % SparseData class: data of a functional TT approximation
    %
    % SparseData Properties:
    %   I     - A multi-index.
    %   coeff - Coefficients of the approximation.
    %   x     - Sampling points. 
    %   A     - Vandermonde/design matrix.
    %   w     - Weights.
    %   y     - Response variable/funcion outputs.
    %   qA,rA - Factors of QR factorization of A.
    %   err   - Estimated L2 error from cross validation
    %   opt_err - 
    %           Estimated approximated projector norm
    %
    % SparseData Methods:
    %   clean - Cleans intermediate data used for builing SparseFun and  
    %           keeps coeff and I for evlauations.
    % 
    % See also SparseFun
    
    properties
        I MultiIndices
        coeff % coefficients
        x % sample points
        A % Vandermore matrix
        w % weights
        y % function evaluations
        qA % QR of A
        rA % QR of A
        err % estimated error
        opt_err % estimated operator error
    end
    
    methods
        function obj = clean(obj)
            obj.x = [];
            obj.A = [];
            obj.w = [];
            obj.y = [];
            obj.qA = [];
            obj.rA = [];
            obj.err = Inf;
            obj.opt_err = Inf;
        end
        
        function obj = SparseData()
            obj.I = MultiIndices();
            obj.coeff = [];
            obj.x = [];
            obj.A = [];
            obj.w = [];
            obj.y = [];
            obj.qA = [];
            obj.rA = [];
            obj.err = Inf;
            obj.opt_err = Inf;
        end
    end 
end