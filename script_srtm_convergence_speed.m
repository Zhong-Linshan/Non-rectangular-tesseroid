%% Computation of vertical gravity over a 1°×1° SRTM block:
%  accuracy and speed comparison of different integration schemes.
%
%  This script loads a 1201×1201 elevation matrix from SRTM tile
%  srtm_57_06.tif (Qinghai–Xizang Plateau region), computes a highly
%  accurate reference solution using Vnorth_pole on a 3″ grid, and
%  then compares the performance (relative error and computation time)
%  of three quadrature families:
%    - Gauss–Legendre (GLQ) with variable or constant top
%    - Closed Newton–Cotes (NCR_Close)
%    - Open Newton–Cotes (NCR_Open)
%
%  For the GLQ schemes the number of quadrature nodes N is varied from
%  2 to 600 to study convergence.  For the Newton–Cotes schemes a
%  fixed 6‑point rule is used.  Each timing is repeated M = 100 times
%  and the mean of the last 90 runs (after 10 warm‑up runs) is reported.

datetime('now')   % display start time

% -------------------------------------------------------------------------
% 1.  Load topography
% -------------------------------------------------------------------------
h_qingzang = imread('srtm_57_06.tif');     % read GeoTIFF
H = double(h_qingzang(1:1201, 1:1201));    % extract 1201×1201 window (3″ grid)

% -------------------------------------------------------------------------
% 2.  Physical and geometrical constants
% -------------------------------------------------------------------------
G  = 6.67349e-11 * 2670;    % gravitational constant × rock density (SI)
r1 = 6371000;               % bottom radius (Earth's mean radius) [m]
r  = 6471000;               % observation radius (100 km above surface) [m]
% IMPORTANT: GLQ_L, NCR_Close, NCR_Open (version B) expect angles in DEGREES,
% so phi and lam are given in degrees.  They are converted to radians
% internally by those functions.
phi = 90;                   % observation latitude  = north pole [deg]
lam = 0;                    % observation longitude = arbitrary    [deg]

% Geographic boundaries of the 1°×1° test area (degrees)
phi1 = 34;                  % southern latitude
phi2 = 35;                  % northern latitude
lam1 = 100;                 % western longitude
lam2 = 101;                 % eastern longitude

% -------------------------------------------------------------------------
% 3.  Reference solution (Vnorth_pole on a 3″ grid)
% -------------------------------------------------------------------------
% Build the tesseroid bottom radius: flip H so that rows run south→north
H_tess = r1 + flip(H);

% Average the four corners of each 1″ sub‑cell for the analytical reference
H_tess1 = (H_tess(1:end-1, 1:end-1) + ...
           H_tess(2:end,   1:end-1) + ...
           H_tess(1:end-1, 2:end)   + ...
           H_tess(2:end,   2:end)) / 4;

% Sub‑grid for the reference: spacing 1/1200° (3 arc‑seconds)
% Vnorth_pole requires radians, so we convert the grid here.
[lamg, phig] = meshgrid(deg2rad(lam1 : 1/1200 : lam2), ...
                        deg2rad(phi1 : 1/1200 : phi2));

% Analytical integration for each sub‑tesseroid
aa = Vnorth_pole(r, r1, H_tess1, ...
                 phig(1:end-1, 1:end-1), phig(2:end, 2:end), ...
                 lamg(1:end-1, 1:end-1), lamg(2:end, 2:end));

% Total vertical gravity (multiplied by G×ρ) – reference value
u = G * sum(aa(:));

% -------------------------------------------------------------------------
% 4.  Experiment setup for GLQ with varying node numbers
% -------------------------------------------------------------------------
N = [2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 25, 30, 40, 48, ...
     50, 60, 75, 80, 100, 120, 150, 200, 240, 300, 400, 600];
M = 1;                    % number of timing repeats

% Pre‑allocate storage
du1 = zeros(length(N), 1);  % relative error  – GLQ variable top
du2 = zeros(length(N), 1);  % relative error  – GLQ constant top
ta1 = zeros(length(N), 1);  % mean time       – GLQ variable top
ta2 = zeros(length(N), 1);  % mean time       – GLQ constant top

% -------------------------------------------------------------------------
% 5.  Loop over Gauss–Legendre node counts
% -------------------------------------------------------------------------
for i = 1 : length(N)
    disp(i)                 % show current index

    [t, w] = nodes(N(i));   % get nodes and weights for the current N

    t1 = zeros(M, 1);       % timing storage for variable top
    t2 = zeros(M, 1);       % timing storage for constant top

    for m = 1 : M
        % Variable top
        tic
        u1 = G * glq(r, r1, H, phi, lam, phi2, lam1, t, w, 'non');
        t1(m) = toc;

        % Constant top
        tic
        u2 = G * glq(r, r1, H, phi, lam, phi2, lam1, t, w, 'rect');
        t2(m) = toc;
    end

    % Relative errors (values are scalars, same for all runs)
    du1(i) = abs(u1 - u) / u;
    du2(i) = abs(u2 - u) / u;

    % Average time (skip the first 10 warm‑up runs)
    ta1(i) = mean(t1(10:M));
    ta2(i) = mean(t2(10:M));
end

% -------------------------------------------------------------------------
% 6.  Timing and error for Newton–Cotes schemes (fixed n = 6)
%     These functions use their own internal node construction, so we
%     only need to pass the external parameters.  Each is run M times.
% -------------------------------------------------------------------------
t3 = zeros(1, M);           % time – CNCR variable top
t4 = zeros(1, M);           % time – CNCR constant top
t5 = zeros(1, M);           % time – ONCR variable top
t6 = zeros(1, M);           % time – ONCR constant top

for m = 1 : M
    % Closed Newton–Cotes, variable top
    tic
    u3 = G * cncr(r, r1, H, phi, lam, phi2, lam1, 'non');
    t3(m) = toc;

    % Closed Newton–Cotes, constant top
    tic
    u4 = G * cncr(r, r1, H, phi, lam, phi2, lam1, 'rect');
    t4(m) = toc;

    % Open Newton–Cotes, variable top
    tic
    u5 = G * oncr(r, r1, H, phi, lam, phi2, lam1, 'non');
    t5(m) = toc;

    % Open Newton–Cotes, constant top
    tic
    u6 = G * oncr(r, r1, H, phi, lam, phi2, lam1, 'rect');
    t6(m) = toc;
end

% Relative errors (scalar, same for all M runs)
du3 = abs(u3 - u) / u;
du4 = abs(u4 - u) / u;
du5 = abs(u5 - u) / u;
du6 = abs(u6 - u) / u;

% Mean times (excluding warm‑up)
ta3 = mean(t3(10:M));
ta4 = mean(t4(10:M));
ta5 = mean(t5(10:M));
ta6 = mean(t6(10:M));

datetime('now')   % display end time


% =========================================================================
%  Local subfunctions: wrappers that partition a 1°×1° block into
%  sub‑tesseroids of size inc×inc and apply a quadrature method.
%  The top surface may be variable (the actual SRTM heights) or
%  constant (the mean height).
% =========================================================================

function u = glq(r, r1, H, phi, lam, phiq, lamq, t, w, str)
% glq   Wrapper for Gauss–Legendre quadrature over a 1°×1° block.
%       The block is divided into sub‑tesseroids whose size is
%       inc = n/1200 degrees, where n = length(t).
%       All angles (phi, lam, phiq, lamq, and the computed boundaries)
%       are in DEGREES, as required by GLQ_L.
    n = length(t);
    inc = n / 1200;
    v = zeros(1200, 1200);
    for i = 1 : 1200/n
        phi1 = phiq - i * inc;      % southern boundary [deg]
        phi2 = phi1 + inc;          % northern boundary [deg]
        for j = 1 : 1200/n
            lam1 = lamq + (j-1) * inc;   % western boundary [deg]
            lam2 = lam1 + inc;           % eastern boundary [deg]
            % Extract the corresponding sub‑block and flip to N→S order
            H_tess = r1 + flip(H((i-1)*n+1 : i*n+1, ...
                                  (j-1)*n+1 : j*n+1));
            v(i, j) = GLQ_L(lam, lam1, lam2, phi, phi1, phi2, ...
                            r, r1, H_tess, t, w, str);
        end
    end
    u = sum(v(:));
end

function u = cncr(r, r1, H, phi, lam, phiq, lamq, str)
% cncr  Wrapper for closed Newton–Cotes quadrature (fixed n = 6).
%       Uses NCR_Close with n+1 points (7 points per dimension).
%       All angles are in DEGREES, consistent with NCR_Close.
    n = 6;
    inc = n / 1200;
    v = zeros(1200, 1200);
    for i = 1 : 1200/n
        phi1 = phiq - i * inc;
        phi2 = phi1 + inc;
        for j = 1 : 1200/n
            lam1 = lamq + (j-1) * inc;
            lam2 = lam1 + inc;
            H_tess = r1 + flip(H((i-1)*n+1 : i*n+1, ...
                                  (j-1)*n+1 : j*n+1));
            v(i, j) = NCR_Close(lam, lam1, lam2, phi, phi1, phi2, ...
                                r, r1, H_tess, n+1, str);
        end
    end
    u = sum(v(:));
end

function u = oncr(r, r1, H, phi, lam, phiq, lamq, str)
% oncr  Wrapper for open Newton–Cotes quadrature (fixed n = 6).
%       Uses NCR_Open with n points.  For the variable top mode the
%       averaged corner heights H_tess1 are used, matching the
%       convention of the reference solution.
%       All angles are in DEGREES, consistent with NCR_Open.
    n = 6;
    inc = n / 1200;
    v = zeros(1200, 1200);
    for i = 1 : 1200/n
        phi1 = phiq - i * inc;
        phi2 = phi1 + inc;
        for j = 1 : 1200/n
            lam1 = lamq + (j-1) * inc;
            lam2 = lam1 + inc;
            H_tess = r1 + flip(H((i-1)*n+1 : i*n+1, ...
                                  (j-1)*n+1 : j*n+1));
            % Average the four corners for the open Newton–Cotes version
            H_tess1 = (H_tess(1:end-1,1:end-1) + H_tess(2:end,1:end-1) + ...
                       H_tess(1:end-1,2:end)   + H_tess(2:end,2:end)) / 4;
            v(i, j) = NCR_Open(lam, lam1, lam2, phi, phi1, phi2, ...
                               r, r1, H_tess1, n, str);
        end
    end
    u = sum(v(:));
end