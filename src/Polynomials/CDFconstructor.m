function cdf = CDFconstructor(poly, varargin)
% Select the one dimensional CDF function for a given collocation basis. 
%   cdf = CDFconstructor(poly)

if isa(poly, 'Lagrange1')
    cdf = Lagrange1CDF(poly, varargin{:});
elseif isa(poly, 'Lagrangep')
    cdf = LagrangepCDF(poly, varargin{:});
elseif isa(poly, 'Legendre')
    cdf = BoundedPolyCDF(poly, varargin{:});
elseif isa(poly, 'Chebyshev1st')
    cdf = Chebyshev1stTrigoCDF(poly, varargin{:});
    %cdf = Chebyshev1stCDF(poly, varargin{:});
elseif isa(poly, 'Chebyshev2nd')
    cdf = Chebyshev2ndTrigoCDF(poly, varargin{:});
    %cdf = Chebyshev1stCDF(poly, varargin{:});
elseif isa(poly, 'Jacobi11')
    error('not implemented for general Jacobi in [-1,1]');
elseif isa(poly, 'Fourier')
    cdf = FourierCDF(poly, varargin{:});
elseif isa(poly, 'Hermite')
    cdf = HermiteCDF(poly, varargin{:});
elseif isa(poly, 'Laguerre') 
    cdf = LaguerreCDF(poly, varargin{:});
end

end