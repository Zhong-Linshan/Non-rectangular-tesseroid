function vGLD = NCR_Close(x, x1, x2, y, y1, y2, z, z1, H, n, str)
% NCR_Close  Evaluates a spherical kernel integral over a tesseroid using
%            closed Newton–Cotes quadrature.
%
%   This function computes the same volume integral as GLQ_L, but
%   replaces Gauss–Legendre quadrature with an equally‑spaced
%   Newton–Cotes rule. The number of integration points n in each
%   dimension can be chosen as 2, 3, 4, 5, 6, or 7 (higher values are
%   not implemented in the default weight table).
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
%   Notes:
%       * The integration uses the closed Newton–Cotes weights provided
%         by the internal function nodes().  The weights already sum to
%         unity, so the total weight for each dimension is multiplied by
%         the interval length.
%       * The radial nodes are equally spaced between z1 and z2 (the
%         top of the tesseroid).  This differs from the Gauss–Legendre
%         approach, where the nodes are mapped from [-1,1].
%       * Because Newton–Cotes weights can become negative for n >= 9,
%         only n = 2…7 are supported.
%
%   Example:
%       % Evaluate the kernel with n = 5 points in each direction
%       H = 6371e3 * ones(5, 5);   % constant top at 6371 km
%       val = NCR_Close(0, -1, 1, 45, 44, 46, 6371e3, 6370e3, H, 5);
%
%   See also: GLQ_L, MESHGRID, NDGRID, DEG2RAD.

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
    error('NCR_Close:InvalidSurfaceMode', ...
          'Surface mode must be ''rect''/''constant'' or ''non''/''variable''.');
end

% ---------------------------------------------------------------------
% Newton–Cotes weights and nodes
% ---------------------------------------------------------------------
w = nodes(n);               % column vector of length n, sum(w) = 1

% Equally spaced nodes in each dimension (in original units)
lon_nodes = linspace(x1, x2, n)';   % column vector (degrees)
lat_nodes = linspace(y1, y2, n)';   % column vector (degrees)

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
% Radial nodes (equally spaced from z1 to z2)
% ---------------------------------------------------------------------
% t is a 1x1xn array with values 0, 1, 2, ..., n-1
t_3d = reshape(0:n-1, 1, 1, n);
ZQ   = z1 + t_3d .* (z2 - z1) / (n - 1);    % n x n x n array

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
% cos(angular distance) between observation point and integration points
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
function [w] = nodes(n)
% nodes   Closed Newton–Cotes quadrature weights.
%
%   w = nodes(n) returns the column vector of weights for an n‑point
%   closed Newton–Cotes rule on the standard interval [0, 1].
%   The weights are normalised so that sum(w) = 1.
%   Supported values: n = 2 (trapezoidal), 3 (Simpson's 1/3),
%   4 (Simpson's 3/8), 5, 6, 7. For any other n a single weight
%   w = 1 is returned (midpoint rule, not strictly Newton–Cotes).
%
%   Reference:
%   Abramowitz, M. & Stegun, I. A. (1964), Handbook of Mathematical
%   Functions, Ch. 25.4.

switch n
    case 2
        w = [1/2; 1/2];
    case 3
        w = [1; 4; 1] / 6;
    case 4
        w = [1; 3; 3; 1] / 8;
    case 5
        w = [7; 32; 12; 32; 7] / 90;
    case 6
        w = [19/288; 25/96; 25/144; 25/144; 25/96; 19/288];
    case 7
        w = [41/840; 9/35; 9/280; 34/105; 9/280; 9/35; 41/840];
    otherwise
        % fallback (not a standard Newton–Cotes rule)
        w = 1;
end
end