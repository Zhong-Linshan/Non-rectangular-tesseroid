%% Numerical experiment: accuracy comparison of six integration methods
%  for the vertical gravity of a tesseroid observed at the north pole.
%
%  This script performs M = 300 independent random realisations of a
%  synthetic topography and, for each realisation, computes the total
%  vertical gravity over a 1°×1° area using one analytical reference
%  and six different quadrature schemes.  The relative errors are
%  saved and plotted to compare the accuracy of the methods.
%
%  Requirements:
%    - Functions: GLQ_L, NCR_Close, NCR_Open, Vnorth_pole, nodes
%    - External utility: tight_subplot (File Exchange)

% =========================================================================
% 1.  General settings
% =========================================================================
n = 6;                  % number of quadrature nodes per dimension
inc = 1/60;             % grid spacing of the observation grid (degrees)
inc0 = 1/60 / n;        % sub‑grid spacing used for the analytical reference

% Obtain Gauss–Legendre nodes and weights for the GLQ method
[t, w] = nodes(n);      % t, w on [-1, 1]

% Geometry constants
r1 = 6371000;           % bottom radius of the tesseroid (Earth's surface)
r  = 6471000;           % observation radius (100 km above the surface)
phi = 90;               % observation latitude  (north pole, degrees)
lam = 0;                % observation longitude (arbitrary, degrees)

% Number of independent random realisations
M = 300;

% Arrays to store the total gravity sums and relative errors
u  = zeros(M, 6);       % total sum for each method
du = zeros(M, 6);       % relative error w.r.t. analytical reference

% =========================================================================
% 2.  Main loop over random realisations
% =========================================================================
for m = 1 : M
    m   % display the current realisation index

    % -----------------------------------------------------------------
    % 2.1  Generate a smooth synthetic topography (without noise)
    % -----------------------------------------------------------------
    [x, y] = meshgrid((1:60*n+101)/(60*n+101), ...
                       (1:60*n+101)/(60*n+101));
    h1  = 800  * exp(2  * (-10*(x-0.8).^2 - 8 *(y-0.9).^2));
    h2  = 900  * exp(2  * (-7 *(x-0.9).^2 - 10*(y-0.2).^2));
    h3  = 1200 * exp(4  * (-6 *(x-0.5).^2 - 7 *(y-0.5).^2));
    h4  = 700  * exp(2  * (-9 *(x-0.15).^2 - 6 *(y-0.15).^2));
    h5  = 800  * exp(2  * (-8 *(x-0.2).^2 - 9 *(y-0.85).^2));
    h6  = 600  * exp(10 * (-(x-0.5).^2 - (y-0).^2));
    h7  = 600  * exp(10 * (-(x-0.5).^2 - (y-1).^2));
    h8  = 600  * exp(10 * (-(x-0).^2   - (y-0.5).^2));
    h9  = 600  * exp(10 * (-(x-1).^2   - (y-0.5).^2));
    h10 = 600  * exp(10 * (-(x-0).^2   - (y-0).^2));
    h11 = 600  * exp(10 * (-(x-0).^2   - (y-1).^2));
    h12 = 600  * exp(10 * (-(x-1).^2   - (y-0).^2));
    h13 = 600  * exp(10 * (-(x-1).^2   - (y-1).^2));
    % The random noise and Gaussian filtering are disabled here to keep
    % the topography perfectly smooth; uncomment them for a rougher case.
    % h14 = 100 * randn(60*n+101, 60*n+101);
    % gaussianKernel = fspecial('gaussian', [5,5], 2);
    % H1 = h1+...+h13+h14; H1 = imfilter(H1,gaussianKernel); ...

    H1 = h1 + h2 + h3 + h4 + h5 + h6 + h7 + h8 + h9 + h10 + h11 + h12 + h13;

    % Remove the 50‑point buffer to avoid filtering edge effects
    H_earth = H1(51 : 60*n+51, 51 : 60*n+51);

    % -----------------------------------------------------------------
    % 2.2  Allocate matrices for the observation grid (60×60 points)
    % -----------------------------------------------------------------
    v0 = zeros(60, 60);   % analytical reference (Vnorth_pole)
    v1 = zeros(60, 60);   % GLQ, variable top
    v2 = zeros(60, 60);   % GLQ, constant top (mean)
    v3 = zeros(60, 60);   % Closed Newton–Cotes, variable top
    v4 = zeros(60, 60);   % Closed Newton–Cotes, constant top
    v5 = zeros(60, 60);   % Open Newton–Cotes, variable top
    v6 = zeros(60, 60);   % Open Newton–Cotes, constant top

    % -----------------------------------------------------------------
    % 2.3  Loop over 60×60 observation cells (each 1° × 1°)
    % -----------------------------------------------------------------
    for i = 1 : 60
        % Southern boundary with a random shift (integer 0..179).
        % Note: the shift is added directly in degrees; the original
        % code does NOT divide by 3600, which may lead to large offsets.
        phi1 = 90 - i*inc - randi([0, 179], 1, 1);
        phi2 = phi1 + inc;   % northern boundary

        for j = 1 : 60
            % Western boundary with a random shift (integer 0..358).
            % Same remark as above: the shift is in degrees, not arc-seconds.
            lam1 = randi([0, 358], 1, 1) + (j-1)*inc;
            lam2 = lam1 + inc;   % eastern boundary

            % ---------------------------------------------------------
            % Extract the topography block for this cell
            % ---------------------------------------------------------
            % H_earth is (60*n+1)×(60*n+1); we select the sub‑block
            % that covers the current 1°×1° cell.  Because the
            % original topography has rows running south‑to‑north,
            % we flip it to obtain north‑to‑south ordering.
            H_tess = r1 + flip(H_earth((i-1)*n+1 : i*n+1, ...
                                       (j-1)*n+1 : j*n+1));

            % Average the four corner heights for the reference
            % (each block is subdivided into n×n sub‑cells)
            H_tess1 = (H_tess(1:end-1, 1:end-1) + ...
                       H_tess(2:end,   1:end-1) + ...
                       H_tess(1:end-1, 2:end)   + ...
                       H_tess(2:end,   2:end)  ) / 4;

            % Constant‑top elevation: centre of the 6×6 block (index 4)
            r2 = r1 + H_earth((i-1)*n+4, (j-1)*n+4);

            % Sub‑grid for the analytical reference
            [lamg, phig] = meshgrid(deg2rad(lam1 : inc0 : lam2), ...
                                    deg2rad(phi1 : inc0 : phi2));

            % ---------------------------------------------------------
            % (a) Analytical reference using Vnorth_pole
            % ---------------------------------------------------------
            aa = Vnorth_pole(r, r1, H_tess1, ...
                             phig(1:n,   1:n),   phig(2:n+1, 2:n+1), ...
                             lamg(1:n,   1:n),   lamg(2:n+1, 2:n+1));
            v0(i, j) = sum(aa(:));

            % ---------------------------------------------------------
            % (b) Six quadrature‑based methods
            % ---------------------------------------------------------
            %  1: GLQ, variable top
            v1(i, j) = GLQ_L(lam, lam1, lam2, phi, phi1, phi2, ...
                             r, r1, H_tess, t, w, 'non');

            %  2: GLQ, constant top
            v2(i, j) = GLQ_L(lam, lam1, lam2, phi, phi1, phi2, ...
                             r, r1, r2, t, w);

            %  3: Closed Newton–Cotes, variable top  (n+1 points)
            v3(i, j) = NCR_Close(lam, lam1, lam2, phi, phi1, phi2, ...
                                 r, r1, H_tess, n+1, 'non');

            %  4: Closed Newton–Cotes, constant top
            v4(i, j) = NCR_Close(lam, lam1, lam2, phi, phi1, phi2, ...
                                 r, r1, r2, n+1);

            %  5: Open Newton–Cotes, variable top
            v5(i, j) = NCR_Open(lam, lam1, lam2, phi, phi1, phi2, ...
                                r, r1, H_tess1, n, 'non');

            %  6: Open Newton–Cotes, constant top
            v6(i, j) = NCR_Open(lam, lam1, lam2, phi, phi1, phi2, ...
                                r, r1, r2, n);
        end
    end

    % -----------------------------------------------------------------
    % 2.4  Compute the total sum over all 60×60 cells
    % -----------------------------------------------------------------
    u0 = sum(v0(:));   % analytical reference total

    u(m, :)  = [sum(v1(:)), sum(v2(:)), sum(v3(:)), ...
                sum(v4(:)), sum(v5(:)), sum(v6(:))];

    % Relative error with respect to the reference
    du(m, :) = abs(u(m, :) - u0) / u0;
end

% =========================================================================
% 3.  Statistical summary of the relative errors
% =========================================================================
mu = mean(du)';         % mean relative error
ma = max(du)';          % maximum relative error
mi = min(du)';          % minimum relative error
st = std(du)';          % standard deviation

% =========================================================================
% 4.  Plot the logarithmic relative errors
% =========================================================================
figure
set(0, 'DefaultFigureRenderer', 'zbuffer');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10 10 18 10]);

% tight_subplot creates a subplot with minimal margins
ha = tight_subplot(1, 1, [.0 .0], [.11 .02], [.07 .01]);

% Plot each method with a distinct symbol and colour
plot(log10(du(:,1)), 'p', 'Color', [0   255 0  ]/255)  % non GLQ
hold on
plot(log10(du(:,2)), 'o', 'Color', [153 102 51 ]/255)  % rect GLQ
plot(log10(du(:,3)), 's', 'Color', [255 165 0  ]/255)  % non CNCR
plot(log10(du(:,4)), '+', 'Color', [0   0   0  ]/255)  % rect CNCR
plot(log10(du(:,5)), 'd', 'Color', [255 0   50 ]/255)  % non ONCR
plot(log10(du(:,6)), 'h', 'Color', [0   100 255]/255)  % rect ONCR

% Axis limits and ticks
set(gca, 'XTick', 0:30:300, 'XTickLabel', 0:30:300, ...
         'YTick', -10:1:-1, 'YTickLabel', -10:1:-1, ...
         'XLim', [-5 305], 'YLim', [-9 -2.8])

% Legend with six columns, placed slightly to the left
hLegend = legend('non GLQ', 'rect GLQ', 'non CNCR', 'rect CNCR', ...
                 'non ONCR', 'rect ONCR');
hLegend.Position = hLegend.Position + [-0.4 0.1 0 0];
hLegend.NumColumns = 6;
hLegend.ItemTokenSize = [15, 9];
hLegend.Box = "off";

hXLabel = xlabel('Number of experimental replicates');
hYLabel = ylabel('log_{10}(Relative error)');

set(gca, 'FontName', 'Arial', 'FontSize', 10)
set([hLegend, hXLabel, hYLabel], 'FontSize', 10, 'FontName', 'Arial')

% Save the figure at 300 dpi
print(gcf, 'figure3.png', '-r300', '-dpng');