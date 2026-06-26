%% Real SRTM topography processing and gravity forward modelling
%
%  This script loads 3‑arcsecond SRTM elevation tiles for four regions
%  (Grand Canyon, Qinghai–Xizang Plateau, Mount Qomolangma, and Sichuan
%  Basin), computes basic statistics, and sets up parameters for a
%  tesseroid gravity forward modelling experiment.  The actual gravity
%  computation loop is commented out but left for reference; it uses
%  the Vnorth_pole analytical reference and several quadrature methods.
%
%  Requirements:
%    - SRTM GeoTIFF files: srtm_15_05.tif, srtm_57_06.tif,
%      srtm_54_07.tif, srtm_58_06.tif, srtm_57_07.tif, srtm_58_07.tif
%    - Functions: Vnorth_pole, GLQ_L, NCR_Close, NCR_Open, nodes
%      (and possibly GLQ_ONCR, if uncommented)

% =========================================================================
% 1.  Load SRTM data for the four regions
% =========================================================================
% Pre‑allocate a 4801×4801×4 array to hold the elevation (metres)
H = zeros(4801, 4801, 4);

% (a) Grand Canyon region (tile srtm_15_05)
h_Colorado = imread('srtm_15_05.tif');
% Extract a 4801×4801 window: rows 1:4801, columns 1200:6000
H(:, :, 1) = double(h_Colorado(1:4801, 1200:6000));

% (b) Qinghai–Xizang (Tibetan) Plateau (tile srtm_57_06)
h_qingzang = imread('srtm_57_06.tif');
% Full 4801×4801 tile
H(:, :, 2) = double(h_qingzang(1:4801, 1:4801));

% (c) Mount Qomolangma (Everest) region (tile srtm_54_07)
h_zhufeng = imread('srtm_54_07.tif');
H(:, :, 3) = double(h_zhufeng(1:4801, 1:4801));

% (d) Sichuan Basin region – constructed from four adjacent tiles
h_sichuan = [imread('srtm_57_06.tif') imread('srtm_58_06.tif'); ...
             imread('srtm_57_07.tif') imread('srtm_58_07.tif')];
% Extract a 4801×4801 window: rows 3600:8400, columns 4800:9600
H(:, :, 4) = double(h_sichuan(3600:8400, 4800:9600));

% =========================================================================
% 2.  Compute global statistics over all regions
% =========================================================================
% Mean and standard deviation of all pixels (four regions combined)
H_age = mean(mean(H));   % equivalent to mean(H(:))
H_std = std(std(H));     % equivalent to std(H(:))

% =========================================================================
% 3.  Statistics for the Mount Qomolangma region (index 3)
%     Separate elevations above and below 500 m
% =========================================================================
h1 = 0;   % cumulative elevation for points with elevation > 500 m
h2 = 0;   % cumulative elevation for points with elevation ≤ 500 m
n1 = 0;   % count of points with elevation > 500 m
n2 = 0;   % count of points with elevation ≤ 500 m

for i = 1 : 4800
    for j = 1 : 4800
        aa = H(i, j, 3);
        if aa > 500
            n1 = n1 + 1;
            h1 = h1 + aa;
        else
            n2 = n2 + 1;
            h2 = h2 + aa;
        end
    end
end

% Mean elevation for the two categories
h1_age = h1 / n1;   % mean above 500 m
h2_age = h2 / n2;   % mean below or equal to 500 m

% =========================================================================
% 4.  Set up parameters for the gravity forward modelling
% =========================================================================
n = 6;                      % number of quadrature nodes per dimension
[t, w] = nodes(n);          % Gauss–Legendre nodes and weights on [-1,1]
inc = n / 1200;             % grid spacing of the observation grid (degrees)
inc0 = 1 / 1200;            % sub‑grid spacing for the analytical reference (degrees)

% Physical constants
G = 6.67349e-11 * 2670;     % gravitational constant × density (SI units)

% Geometry (observation point 100 km above the reference sphere)
r1 = 6371000;               % Earth's mean radius (bottom of tesseroids)
r  = 6471000;               % observation radius

% Observation point coordinates (north pole, in DEGREES for the quadrature
% routines, which internally convert to radians)
phi = 90;                   % latitude  (degrees)
lam = 0;                    % longitude (degrees)

% Approximate central coordinates of the four regions (degrees)
PHI = [40; 35; 30; 32];     % latitude  (degrees)
LAM = [251; 100; 85; 104];  % longitude (degrees)

% =========================================================================
% 5.  Topography visualisation (Figure 4)
%     These commands plot each region as a mesh surface, viewed from
%     south‑west, with latitude and longitude tick labels.
% =========================================================================
str1 = {'\bf(a) \rmGrand Canyon'; ...
        '\bf(b) \rmQinghai–Xizang Plateau'; ...
        '\bf(c) \rmMount Qomolangma'; ...
        '\bf(d) \rmSichuan Basin'};

for i = 1 : 4
    figure
    [x, y] = meshgrid(PHI(i)-4 : 1/1200 : PHI(i), ...
                      LAM(i)    : 1/1200 : LAM(i)+4);
    meshc(x, y, H(:, :, i)');
    view(45, 50)
    set(gca, 'XTick', PHI(i)-4:1:PHI(i), ...
             'XTickLabel', [PHI(i) PHI(i)-1 PHI(i)-2 PHI(i)-3 PHI(i)-4])
    xlabel('Latitude unit: \circ', 'FontSize', 10)
    ylabel('Longitude unit: \circ', 'FontSize', 10)
    zlabel('Elevation (m)', 'FontSize', 10)
    title(str1(i), 'FontSize', 10)
    % Print with a correct file name (e.g., figure4a.png)
    print(gcf, sprintf('figure4%c.png', 'a'+i-1), '-r300', '-dpng');
end

% =========================================================================
% 6.  Gravity forward modelling loop
%     For each region, the vertical gravity of 800×800 tesseroids is
%     computed using the analytical reference (Vnorth_pole) and several
%     quadrature methods (GLQ, NCR_Close, NCR_Open, and a constant‑top
%     GLQ).  The total gravity and relative errors are stored.
% =========================================================================
u0 = zeros(4, 1);   % analytical reference total
u  = zeros(4, 4);   % total gravity for 4 approximate methods
du = zeros(4, 4);   % relative error

for k = 1 : 4
    tic
    k   % display current region index
    H1 = H(:, :, k);
    v0 = zeros(800, 800);   % reference (Vnorth_pole)
    v1 = v0;                 % GLQ, variable top
    v2 = v0;                 % NCR_Close, variable top
    v3 = v0;                 % NCR_Open, variable top
    v4 = v0;                 % GLQ, constant top

    for i = 1 : 800
        % Latitude boundaries of the current tesseroid
        phi1 = PHI(k) - i * inc;      % southern boundary (degrees)
        phi2 = phi1 + inc;            % northern boundary

        for j = 1 : 800
            % Longitude boundaries
            lam1 = LAM(k) + (j-1) * inc;   % western boundary (degrees)
            lam2 = lam1 + inc;             % eastern boundary

            % Topography block for this tesseroid (flipped to north‑south)
            H_tess = r1 + flip(H1((i-1)*n+1 : i*n+1, ...
                                   (j-1)*n+1 : j*n+1));
            % Averaged corner heights for the reference sub‑division
            H_tess1 = (H_tess(1:end-1,1:end-1) + H_tess(2:end,1:end-1) + ...
                       H_tess(1:end-1,2:end)   + H_tess(2:end,2:end)) / 4;

            % Constant‑top radius (centre of the 6×6 block)
            r2 = r1 + H1((i-1)*n+4, (j-1)*n+4);

            % Sub‑grid for the analytical reference
            [lamg, phig] = meshgrid(deg2rad(lam1 : inc0 : lam2), ...
                                    deg2rad(phi1 : inc0 : phi2));

            % Analytical reference (Vnorth_pole expects radians)
            aa = Vnorth_pole(r, r1, H_tess1, ...
                             phig(1:n,1:n), phig(2:n+1,2:n+1), ...
                             lamg(1:n,1:n), lamg(2:n+1,2:n+1));
            v0(i, j) = sum(aa(:));

            % Approximate methods (all expect angles in DEGREES)
            v1(i,j) = GLQ_L(lam,lam1,lam2,phi,phi1,phi2,r,r1,H_tess,t,w,'non');
            v2(i,j) = NCR_Close(lam,lam1,lam2,phi,phi1,phi2,r,r1,H_tess,n+1,'non');
            v3(i,j) = NCR_Open(lam,lam1,lam2,phi,phi1,phi2,r,r1,H_tess1,n,'non');
            v4(i,j) = GLQ_L(lam,lam1,lam2,phi,phi1,phi2,r,r1,r2,t,w);
        end
    end

    % Total gravity (multiplied by G * density)
    u0(k) = G * sum(v0(:));
    u(k,:) = [G*sum(v1(:)), G*sum(v2(:)), G*sum(v3(:)), ...
              G*sum(v4(:))];
    du(k,:) = abs(u(k,:) - u0(k)) ./ u0(k);
    toc
end