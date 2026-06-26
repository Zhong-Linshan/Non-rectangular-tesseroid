function vGLD = GLQ_L(x, x1, x2, y, y1, y2, z, z1, H, t, w, str)
% GLQ_L  Evaluates the volume integral of a spherical kernel over a 
%        tesseroid (spherical prism) using Gauss-Legendre quadrature.
%
%   This function computes the integral of a kernel of the form
%         K = ZQ^2 * cos(lat) / distance
%   over a tesseroid defined by longitude bounds [x1, x2], latitude bounds
%   [y1, y2], and radial bounds [z1, z2]. The top surface z2 can be either
%   constant (mean of H) or spatially variable (interpolated from H).
%   The result is often used as a building block for gravitational
%   potential, gravity, or gravity gradient modelling in spherical
%   coordinates (e.g., tesseroid forward modelling).
%
%   Input arguments:
%       x   - Longitude of the observation point (decimal degrees)
%       x1  - Western longitude boundary of the tesseroid (degrees)
%       x2  - Eastern longitude boundary of the tesseroid (degrees)
%       y   - Latitude of the observation point (decimal degrees)
%       y1  - Southern latitude boundary of the tesseroid (degrees)
%       y2  - Northern latitude boundary of the tesseroid (degrees)
%       z   - Radial coordinate of the observation point (e.g., radius from
%             Earth's centre, in metres)
%       z1  - Radial coordinate of the bottom surface of the tesseroid (m)
%       H   - 2D matrix of the top surface radial coordinates (m).
%             The matrix dimensions must correspond to a regular grid
%             covering [x1, x2] in longitude and [y1, y2] in latitude.
%             By default, H is expected to be arranged such that rows
%             correspond to latitude and columns to longitude.
%       t   - Vector of Gauss-Legendre nodes (length N, values in [-1, 1])
%       w   - Vector of Gauss-Legendre weights (length N)
%       str - (Optional) character vector specifying how the top surface
%             radius is determined:
%               'constant' | 'rect'  : use the mean value of H (default)
%               'variable' | 'non'   : use spatially variable top surface
%                                       by linear interpolation of H
%
%   Output:
%       vGLD - Scalar value of the integrated kernel. Multiply by the
%              product of the density and the gravitational constant to
%              obtain physical gravity quantities (the exact factor depends
%              on the modelling context).
%
%   Notes:
%       * The integration is performed in spherical coordinates. The
%         trigonometric terms convert angular distances and include the
%         cosine of latitude for proper area element weighting.
%       * Gauss-Legendre quadrature of order N is applied in all three
%         dimensions using the provided nodes and weights.
%       * The code uses implicit array expansion and therefore requires
%         MATLAB R2016b or later.
%
%   Example:
%       N = 4;
%       [t, w] = lgwt(N, -1, 1);          % obtain nodes and weights
%       H = 6371e3 * ones(101, 101);      % constant top radius (in metres)
%       val = GLQ_L(0, -1, 1, 0, -1, 1, 6371e3, 6370e3, H, t, w, 'const');
%
%   References:
%       [1] Heck, B., & Seitz, K. (2007). A comparison of the tesseroid,
%           prism and point-mass approaches for mass reductions in gravity
%           field modelling. Journal of Geodesy, 81, 121-136.
%       [2] Uieda, L., Barbosa, V. C. F., & Braitenberg, C. (2016).
%           Tesseroids: forward modelling in spherical coordinates.
%           Geophysical Journal International, 204(2), 1164-1179.
%
%   See also: MESHGRID, NDGRID, INTERP2, DEG2RAD, LGWT.

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
narginchk(11, 12);

% Set default top-surface mode
if nargin < 12 || isempty(str)
    str = 'constant';
end

% Validate and parse the surface mode (accept both new and legacy keywords)
str = lower(str);
switch str
    case {'constant', 'rect'}
        useVariableTopo = false;
    case {'variable', 'non'}
        useVariableTopo = true;
    otherwise
        error('GLQ_L:InvalidSurfaceMode', ...
              ['Unrecognised surface mode: ''%s''.\n', ...
               'Valid options are ''constant''/''rect'' or ''variable''/''non''.'], str);
end

% -------------------------------------------------------------------------
% Quadrature setup
% -------------------------------------------------------------------------
N = length(t);          % Number of quadrature nodes per dimension

% Map nodes from [-1, 1] to the physical coordinate ranges
lonNodes = t * (x2 - x1) / 2 + (x2 + x1) / 2;   % longitude nodes (Nx1)
latNodes = t * (y2 - y1) / 2 + (y2 + y1) / 2;   % latitude nodes  (Nx1)

% -------------------------------------------------------------------------
% Top surface elevation (z2)
% -------------------------------------------------------------------------
if useVariableTopo
    % Build grid matching the original H matrix for interpolation
    % Assume H has dimensions: (nLatitude, nLongitude) i.e. rows = latitude
    [nLat, nLon] = size(H);
    [lonH, latH] = meshgrid(linspace(x1, x2, nLon), ...
                             linspace(y1, y2, nLat));

    % Quadrature node grids for interpolation
    [lonGrid2d, latGrid2d] = meshgrid(lonNodes, latNodes);

    % Linearly interpolate the top surface at the quadrature nodes
    z2 = interp2(lonH, latH, H, lonGrid2d, latGrid2d, 'linear');
else
    % Constant top surface: use the mean value of H
    z2 = mean(H(:)) * ones(N, N);
end

% -------------------------------------------------------------------------
% Radial nodes (3D)
% -------------------------------------------------------------------------
% Expand Gauss-Legendre nodes to the third dimension
t_3d = reshape(t, 1, 1, N);              % [1 x 1 x N]

% Map radial nodes to [z1, z2]  (implicit expansion: z2 is NxN, t_3d is 1x1xN)
ZQ = t_3d .* (z2 - z1) / 2 + (z2 + z1) / 2;   % [N x N x N]

% -------------------------------------------------------------------------
% Integration weights for each dimension
% -------------------------------------------------------------------------
% Longitude: weight * (interval half-length) * (degree-to-radian factor)
wx = w * deg2rad((x2 - x1) / 2);
% Latitude
wy = w * deg2rad((y2 - y1) / 2);
% Radial (the interval half-length factor is applied later)
wz = w;

% -------------------------------------------------------------------------
% Compute the kernel on the 3D quadrature grid
% -------------------------------------------------------------------------
% Longitude and latitude nodes in 3D
[X_3d, Y_3d, ~] = ndgrid(lonNodes, latNodes, ones(N, 1));

% Convert angular coordinates to radians
XQ = deg2rad(X_3d);    % longitude of integration points [rad]
YQ = deg2rad(Y_3d);    % latitude  of integration points [rad]

% Observation point coordinates in radians
xObsRad = deg2rad(x);
yObsRad = deg2rad(y);

% Spherical distance between observation point and integration points
% cos(psi) = sin(phi_obs)*sin(phi) + cos(phi_obs)*cos(phi)*cos(lambda - lambda_obs)
cosPsi = sin(yObsRad) .* sin(YQ) + ...
         cos(yObsRad) .* cos(YQ) .* cos(XQ - xObsRad);

% L = sqrt(r_obs^2 + r^2 - 2 * r_obs * r * cos(psi))
distance = sqrt(z.^2 + ZQ.^2 - 2 .* z .* ZQ .* cosPsi);

% Kernel: ZQ^2 * cos(latitude) / distance
K = ZQ.^2 .* cos(YQ) ./ distance;

% -------------------------------------------------------------------------
% Apply Gaussian quadrature weights
% -------------------------------------------------------------------------
% Build 3D arrays of weights
[WX, WY, Wz] = ndgrid(wx, wy, wz);

% Radial weights additionally multiplied by the interval half-length
WZ = Wz .* (z2 - z1) / 2;

% Final integration (absolute value as per the original algorithm)
vGLD = abs(sum(K(:) .* WZ(:) .* WY(:) .* WX(:)));

end