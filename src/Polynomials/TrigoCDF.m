classdef TrigoCDF < SpectralCDF
    % SpectralCDF class
    %
    % For Fourier basis, FourierCDF is used. For other spectral polynomials
    % in bounded domains, we first transform the polynomial to the 2nd
    % Chebyshev basis, and then apply the inversion. See Chebyshev2ndCDF.
    %
    % Before applying root findings, a grid search based on sampling_nodes
    % is applied to locate the left and right boundary of root finding.
    %
    % See also ChebyshevCDF and FourierCDF.
    
    
    methods
        function obj = TrigoCDF(varargin)
            obj@SpectralCDF(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function x = grid_measure(obj, n)
            x = linspace(-pi, 0, n);
            x(1)   = -pi+eps;
            x(end) = 0-eps;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_int_deri(obj, pk, r)
            % r in [-1,1] to theta in [-pi,0]
            z = eval_int_deri@SpectralCDF(obj, pk, -real(acos(r)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf(obj, pk, r)
            % r in [-1,1] to theta in [-pi,0]
            z = eval_cdf@SpectralCDF(obj, pk, -real(acos(r)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf(obj, pk, xi)
            r = invert_cdf@SpectralCDF(obj, pk, xi);
            r = cos(-r); % theta in [-pi,0] to r in [-1,1]
        end
        
    end
    
end