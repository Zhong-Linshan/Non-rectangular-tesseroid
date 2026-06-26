function result = Vnorth_pole(r, r1, r2, phi1, phi2, lam1, lam2)
% Vnorth_pole  Vertical gravitational attraction of a tesseroid
%              (spherical prism) – analytical integration over
%              longitude and numerical integration over latitude.
%
%   This function computes the vertical component of the gravitational
%   attraction produced by a tesseroid at an observation point located
%   on the polar axis (or, more generally, at any point, because the
%   formula is valid for arbitrary locations after a coordinate
%   rotation).  The expression results from analytically integrating
%   the kernel with respect to longitude and then evaluating the
%   remaining integral with respect to co‑latitude using the auxiliary
%   function Vvr.
%
%   *** ALL ANGULAR INPUTS ARE IN RADIANS ***
%       phi1, phi2 : geocentric latitude  (not co‑latitude)
%       lam1, lam2 : longitude
%
%   Inputs:
%       r     - Radial coordinate of the observation point (metres)
%       r1    - Inner (bottom) radius of the tesseroid (metres)
%       r2    - Outer (top) radius of the tesseroid (metres)
%       phi1  - Southern latitude boundary (radians, -pi/2 … pi/2)
%       phi2  - Northern latitude boundary (radians, -pi/2 … pi/2)
%       lam1  - Western longitude boundary (radians, 0 … 2*pi)
%       lam2  - Eastern longitude boundary (radians, 0 … 2*pi)
%
%   Output:
%       result - Vertical gravitational attraction produced by the
%                tesseroid (units: m/s^2 if r, r1, r2 in metres and
%                multiplied by G*rho; the present function returns
%                only the geometric integral, without G and density).
%
%   Mathematical formulation:
%       result = dlam * [ Vvr(r, r1, r2, pi/2 - phi1) -
%                         Vvr(r, r1, r2, pi/2 - phi2) ]
%       where dlam = lam2 - lam1 and Vvr is the kernel integrated
%       with respect to longitude.
%
%   Example:
%       r_obs = 6371e3;
%       r_bot = 6370e3;
%       r_top = 6371e3;
%       val = Vnorth_pole(r_obs, r_bot, r_top, ...
%                           deg2rad(44), deg2rad(46), ...
%                           deg2rad(-1), deg2rad(1));
%       % val is the geometric factor; multiply by G*rho for gravity.
%
%   Reference:
%       Heck, B., & Seitz, K. (2007). A comparison of the tesseroid,
%       prism and point‑mass approaches … Journal of Geodesy, 81, 121‑136.
%       (see Eq. … for the analytical longitude integration)
%
%   See also: Vvr, GLQ_L, NCR_Close, NCR_Open.

% Convert latitude boundaries to co‑latitude
the2 = pi/2 - phi1;   % co‑latitude of northern boundary
the1 = pi/2 - phi2;   % co‑latitude of southern boundary
dlam = lam2 - lam1;   % longitude extent

% Evaluate the integrated kernel at the two co‑latitude limits
result = dlam .* ( Vvr(r, r1, r2, the2) - Vvr(r, r1, r2, the1) );

end

% =====================================================================
function result = Vvr(r, r1, r2, psi)
% Vvr   Kernel of the analytical longitude integration for a
%       tesseroid's vertical gravity.
%
%   result = Vvr(r, r1, r2, psi) returns the value of the integral
%       ∫ cos(psi) / sqrt(r^2 + r'^2 - 2 r r' cos(psi)) dlambda
%   evaluated analytically after integration over longitude.  The
%   formula is valid for an arbitrary co‑latitude psi (in radians).
%
%   Inputs:
%       r   - Radial coordinate of observation point (metres)
%       r1  - Bottom radius of the tesseroid (metres)
%       r2  - Top radius of the tesseroid (metres).
%             **Important**: r2 must be an array of the same size as
%             psi because the code uses logical indexing on r2.
%             If a scalar is passed, the user must expand it
%             beforehand (e.g., r2 * ones(size(psi))).
%       psi - Co‑latitude (angle from the polar axis) in radians.
%             May be a scalar or an array.
%
%   Output:
%       result - Same size as psi, containing the kernel value.
%
%   Remarks:
%       The function treats the special case psi == 0 (observation
%       point on the polar axis) separately, using a simplified
%       expression that avoids division by zero.
%
%   Reference:
%       Derived from the analytical integration of the tesseroid
%       potential; see, e.g., Heck & Seitz (2007) or Grombein et al.
%       (2013, Journal of Geodesy).
%
%   See also: Vnorth_pole.

% Logical masks for the general case (psi ~= 0) and the polar case
idx = (psi == 0);          % points on the polar axis
idx0 = ~idx;               % all other points

result = zeros(size(psi));

% -----------------------------------------------------------------
% General case (psi ~= 0)
% -----------------------------------------------------------------
% Distances to bottom (r1) and top (r2) surfaces
l1 = sqrt( (r - r1).^2 + 2 * r .* r1 .* (1 - cos(psi(idx0))) );
l2 = sqrt( (r - r2(idx0)).^2 + 2 * r .* r2(idx0) .* (1 - cos(psi(idx0))) );

% Radial projections
ra1 = r1 - r .* cos(psi(idx0));
ra2 = r2(idx0) - r .* cos(psi(idx0));

% Analytic formula for the vertical kernel
result(idx0) = (l2 - l1) .* (l2.^2 + l1.*l2 + l1.^2) / (3 * r) ...
             + (l2 .* ra2 - l1 .* ra1) / 2 .* cos(psi(idx0)) ...
             + (r^2 / 2) .* sin(psi(idx0)).^2 .* cos(psi(idx0)) ...
                .* log( (l2 + ra2) ./ (l1 + ra1) );

% -----------------------------------------------------------------
% Polar case (psi == 0) – simplified expression
% -----------------------------------------------------------------
result(idx) = abs( r2(idx) - r ) .* ( (r2(idx) - r).^2 / (3 * r) ...
                                    + (r2(idx) - r) / 2 ) ...
            - abs( r1 - r ) .* ( (r1 - r).^2 / (3 * r) ...
                               + (r1 - r) / 2 );

end