function vGLD = NCR_Open(x, x1, x2, y, y1, y2, z, z1, H, n, str)
% NCR_Open  Evaluates a spherical kernel integral over a tesseroid using
%           open Newton–Cotes quadrature.
%
%   This function computes the same volume integral as GLQ_L, but
%   replaces Gauss–Legendre quadrature with an *open* Newton–Cotes rule.
%   The integration nodes lie strictly inside the interval (endpoints
%   excluded), which can be advantageous when the integrand is singular
%   or nearly singular at the boundaries. The number of integration
%   points n per dimension may be 2, 3, 4, 5, 6, or 7 (the default is
%   5).
%
%   *** ALL ANGULAR INPUTS ARE IN DEGREES ***
%   (observation coordinates x, y and the tesseroid boundaries x1, x2,
%    y1, y2).  They are internally converted to radians.
%
%   Inputs:
%       x   - Observation longitude (degrees)
%       x1  - West longitude boundary of the tesseroid (degrees)
%       x2  - East longitude boundary of the tesseroid (degrees)
%       y   - Observation latitude (degrees)
%       y1  - South latitude boundary of the tesseroid (degrees)
%       y2  - North latitude boundary of the tesseroid (degrees)
%       z   - Observation radial coordinate (e.g., radius from Earth's
%             centre, metres)
%       z1  - Bottom radial boundary of the tesseroid (metres)
%       H   - Top radial boundary grid. For 'rect'/'constant' mode the
%             mean of H is used. For 'non'/'variable' mode H must be an
%             n‑by‑n matrix of top radii evaluated at the quadrature
%             nodes (the calling code is responsible for providing the
%             correctly formatted H).
%       n   - Number of quadrature nodes per dimension (2…7; default 5)
%       str - (Optional) Surface type:
%               'rect' | 'constant'   : constant top, mean(H(:)) (default)
%               'non'  | 'variable'   : variable top, directly uses H
%
%   Output:
%       vGLD - Scalar value of the integrated kernel. Multiply by
%              G * density to obtain physical gravity quantities.
%
%   Open Newton–Cotes nodes and weights:
%       The nodes are constructed as
%           t_k = (2k - 1) / (2n),  k = 1, 2, ..., n
%       which are the midpoints of the n equal sub‑intervals of [0, 1].
%       The corresponding weights are provided by the local function
%       nodes() and are normalised such that sum(w) = 1.
%       For n = 2…7 these weights are exact for polynomials up to degree
%       n−1 (open Newton–Cotes formulas).
%
%   Example:
%       % Evaluate the kernel with n = 5 (open Newton–Cotes)
%       H = 6371e3 * ones(5, 5);   % constant top at 6371 km
%       val = NCR_Open(0, -1, 1, 45, 44, 46, 6371e3, 6370e3, H, 5);
%
%   See also: GLQ_L, NCR_Close, NODES.

% ---------------------------------------------------------------------
% Input parsing and defaults
% ---------------------------------------------------------------------
narginchk(9, 11);
if nargin == 9
    n = 5;
    str = 'rect';
elseif nargin == 10
    str = 'rect';
end
str = lower(str);

% Validate surface mode
if ~any(strcmp(str, {'rect', 'constant', 'non', 'variable'}))
    error('NCR_Open:InvalidSurfaceMode', ...
          'Surface mode must be ''rect''/''constant'' or ''non''/''variable''.');
end

% ---------------------------------------------------------------------
% Open Newton–Cotes weights and node parameters
% ---------------------------------------------------------------------
w = nodes(n);               % column vector of length n, sum(w) = 1

% Open-node index vector: t = 1, 3, 5, ..., 2n-1
t_idx = 1 : 2 : (2 * n - 1);

% Equally spaced open nodes in each dimension (in degrees)
lon_nodes = x1 + t_idx' * (x2 - x1) / (2 * n);
lat_nodes = y1 + t_idx' * (y2 - y1) / (2 * n);

% ---------------------------------------------------------------------
% Top surface z2
% ---------------------------------------------------------------------
if strcmp(str, 'non') || strcmp(str, 'variable')
    % Variable top: H is assumed to be an n-by-n matrix already
    % evaluated at the quadrature nodes.
    z2 = H;
else
    % Constant top: use mean of H
    z2 = mean(H(:)) * ones(n, n);
end

% ---------------------------------------------------------------------
% Radial nodes (open, equally spaced between z1 and z2)
% ---------------------------------------------------------------------
% t_idx reshaped to 1x1xn to enable implicit expansion
t_3d = reshape(t_idx, 1, 1, n);
ZQ   = z1 + t_3d .* (z2 - z1) / (2 * n);    % n x n x n array

% ---------------------------------------------------------------------
% Weights for the three dimensions
% ---------------------------------------------------------------------
% For longitude: weights * interval length (converted to radians)
wx = w * deg2rad(x2 - x1);
% For latitude: weights * interval length (converted to radians)
wy = w * deg2rad(y2 - y1);
% Radial weights (will be multiplied by the radial interval length later)
wz = w;

% ---------------------------------------------------------------------
% 3D grids of angular nodes (in radians)
% ---------------------------------------------------------------------
[X_3d, Y_3d, ~] = ndgrid(lon_nodes, lat_nodes, ones(n, 1));
XQ = deg2rad(X_3d);         % integration longitudes [rad]
YQ = deg2rad(Y_3d);         % integration latitudes  [rad]

% Observation point coordinates (in radians)
Xobs = deg2rad(x);
Yobs = deg2rad(y);

% ---------------------------------------------------------------------
% Spherical kernel
% ---------------------------------------------------------------------
% cos(angular distance) between observation and integration points
cos_psi = sin(Yobs) .* sin(YQ) + ...
          cos(Yobs) .* cos(YQ) .* cos(XQ - Xobs);

% Distance
dist = sqrt(z.^2 + ZQ.^2 - 2 .* z .* ZQ .* cos_psi);

% Kernel: ZQ^2 * cos(latitude) / distance
K = ZQ.^2 .* cos(YQ) ./ dist;

% ---------------------------------------------------------------------
% Assemble 3D weight arrays and compute the integral
% ---------------------------------------------------------------------
[WX, WY, Wz] = ndgrid(wx, wy, wz);
WZ = Wz .* (z2 - z1);       % multiply by radial interval length

% Summation (absolute value as in the original algorithm)
vGLD = abs(sum(K(:) .* WZ(:) .* WY(:) .* WX(:)));

end

% =====================================================================
function w = nodes(n)
% nodes   Open Newton–Cotes quadrature weights.
%
%   w = nodes(n) returns the column vector of weights for an n‑point
%   *open* Newton–Cotes rule on the standard interval [0, 1]. The
%   formula uses the midpoints of the n subintervals and does not
%   evaluate the integrand at the endpoints.
%   The weights are normalised so that sum(w) = 1.
%   Supported values: n = 2, 3, 4, 5, 6, 7. For any other n a single
%   weight w = 1 is returned (fallback, not a standard rule).
%
%   Reference:
%   Abramowitz, M. & Stegun, I. A. (1964), Handbook of Mathematical
%   Functions, Ch. 25.4.

switch n
    case 2
        w = [1/2; 1/2];
    case 3
        w = [3; 2; 3] / 8;
    case 4
        w = [13; 11; 11; 13] / 48;
    case 5
        w = [275; 100; 402; 100; 275] / 1152;
    case 6
        w = [247; 139; 254; 254; 139; 247] / 1280;
    case 7
        w = [4949/27648; 49/7680; 6223/15360; -6257/34560; ...
             6223/15360; 49/7680; 4949/27648];
    otherwise
        w = 1;
end
end