%% Synthetic topography generation and filtering for tesseroid examples
%  This script creates a high‑resolution synthetic elevation model,
%  applies successive Gaussian filters, and plots the original and
%  filtered surfaces.  The figures are designed to illustrate the
%  effect of smoothing on the topography used in the tesseroid
%  forward modelling examples.

% Number of quadrature nodes per dimension (used to size the grid)
n = 6;

% =========================================================================
% 1.  Build a smooth, multi‑peak synthetic topography
% =========================================================================
% High‑resolution normalised coordinate grid (0 to 1)
[x, y] = meshgrid((1 : 60*n+101) / (60*n+101), ...
                   (1 : 60*n+101) / (60*n+101));

% Eleven Gaussian hills placed at different (x, y) locations
h1  = 800  * exp(2  * (-10*(x-0.8 ).^2 - 8 *(y-0.9 ).^2));
h2  = 900  * exp(2  * (-7 *(x-0.9 ).^2 - 10*(y-0.2 ).^2));
h3  = 1200 * exp(4  * (-6 *(x-0.5 ).^2 - 7 *(y-0.5 ).^2));
h4  = 700  * exp(2  * (-9 *(x-0.15).^2 - 6 *(y-0.15).^2));
h5  = 800  * exp(2  * (-8 *(x-0.2 ).^2 - 9 *(y-0.85).^2));
h6  = 600  * exp(10 * (- (x-0.5 ).^2 -   (y-0   ).^2));
h7  = 600  * exp(10 * (- (x-0.5 ).^2 -   (y-1   ).^2));
h8  = 600  * exp(10 * (- (x-0   ).^2 -   (y-0.5 ).^2));
h9  = 600  * exp(10 * (- (x-1   ).^2 -   (y-0.5 ).^2));
h10 = 600  * exp(10 * (- (x-0   ).^2 -   (y-0   ).^2));
h11 = 600  * exp(10 * (- (x-0   ).^2 -   (y-1   ).^2));
h12 = 600  * exp(10 * (- (x-1   ).^2 -   (y-0   ).^2));
h13 = 600  * exp(10 * (- (x-1   ).^2 -   (y-1   ).^2));

% Add random Gaussian noise to simulate realistic roughness
h14 = 100 * randn(60*n+101, 60*n+101);

% Total synthetic elevation (metres)
H = h1 + h2 + h3 + h4 + h5 + h6 + h7 + h8 + h9 + h10 + h11 + h12 + h13 + h14;

% =========================================================================
% 2.  Progressive Gaussian low‑pass filtering
% =========================================================================
gaussianKernel = fspecial('gaussian', [5, 5], 2);
H1 = imfilter(H,  gaussianKernel);   % filtered once
H2 = imfilter(H1, gaussianKernel);   % filtered twice
H3 = imfilter(H2, gaussianKernel);   % filtered three times

% =========================================================================
% 3.  Extract the central region of interest
% =========================================================================
% The original grid has a buffer of 50 points on each side, which is
% discarded to avoid boundary effects introduced by the Gaussian filter.
H00 = H (51 : 60*n+51, 51 : 60*n+51);   % unfiltered
H01 = H1(51 : 60*n+51, 51 : 60*n+51);   % 1x filtered
H02 = H2(51 : 60*n+51, 51 : 60*n+51);   % 2x filtered
H03 = H3(51 : 60*n+51, 51 : 60*n+51);   % 3x filtered

% Flat reference surface at 700 m (mean base level)
H0 = 700 * ones(60*n+1, 60*n+1);

% =========================================================================
% 4.  Coordinate grid for the final plots
% =========================================================================
[x, y] = meshgrid(linspace(0, 1, 60*n+1), linspace(0, 1, 60*n+1));

% =========================================================================
% 5.  Figure (a): Unfiltered topography (larger view with LaTeX labels)
% =========================================================================
figure
surf(x, y, H00, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on
surf(x, y, H0, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
colormap(turbo)          % perceptually uniform, better than jet
colorbar
lighting gouraud         % smooth shading
camlight headlight       % light from camera position
material dull            % non‑shiny surface

% Axes and tick formatting
set(gca, 'FontName', 'Arial', 'XTick', 0:0.25:1, 'XTickLabel', [], ...
               'YTick', 0:0.25:1, 'YTickLabel', [], 'TickDir', 'out')

% LaTeX axis labels for latitude (φ) and longitude (λ) with prime marks
text('String', '$$\varphi_2$$',       'Interpreter', 'latex', ...
     'Position', [0.18 1.27], 'FontSize', 12);
text('String', '$$\varphi_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 1.12], 'FontSize', 12);
text('String', '$$\varphi_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.88], 'FontSize', 12);
text('String', '$$\varphi_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.63], 'FontSize', 12);
text('String', '$$\varphi_1$$',       'Interpreter', 'latex', ...
     'Position', [0.19 0.28], 'FontSize', 12);

text('String', '$$\lambda_1$$',       'Interpreter', 'latex', ...
     'Position', [0.25 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.5  0.23], 'FontSize', 12);
text('String', '$$\lambda_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.75 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [1.0  0.23], 'FontSize', 12);
text('String', '$$\lambda_2$$',       'Interpreter', 'latex', ...
     'Position', [1.25 0.23], 'FontSize', 12);

view(315, 50)
xlabel('Longitude unit: \circ', 'FontSize', 12)
ylabel('Latitude unit: \circ',   'FontSize', 12)
zlabel('Elevation (m)',          'FontSize', 12)
title('\bf(a) \rmUnfiltered',    'FontSize', 12)
grid on
box on

% =========================================================================
% 6.  Figure (a) repeated: Unfiltered topography (tighter view)
%     This variant is cropped more closely to the data and is used as
%     the baseline panel in the multi‑panel figure 2 of the paper.
% =========================================================================
figure
surf(x, y, H00, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on
surf(x, y, H0, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
colormap(turbo)
colorbar
lighting gouraud
camlight headlight
material dull

set(gca, 'FontName', 'Arial', 'XTick', 0:0.25:1, 'XTickLabel', [], ...
               'YTick', 0:0.25:1, 'YTickLabel', [], 'TickDir', 'out')

% Adjusted label positions for the tighter axes
text('String', '$$\varphi_2$$',       'Interpreter', 'latex', ...
     'Position', [0.073 1.155], 'FontSize', 12);
text('String', '$$\varphi_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [-0.067 1.005], 'FontSize', 12);
text('String', '$$\varphi_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [-0.067 0.765], 'FontSize', 12);
text('String', '$$\varphi_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [-0.067 0.515], 'FontSize', 12);
text('String', '$$\varphi_1$$',       'Interpreter', 'latex', ...
     'Position', [0.083 0.165], 'FontSize', 12);

text('String', '$$\lambda_1$$',       'Interpreter', 'latex', ...
     'Position', [0.14 0.12], 'FontSize', 12);
text('String', '$$\lambda_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.39 0.12], 'FontSize', 12);
text('String', '$$\lambda_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.64 0.12], 'FontSize', 12);
text('String', '$$\lambda_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [0.89 0.12], 'FontSize', 12);
text('String', '$$\lambda_2$$',       'Interpreter', 'latex', ...
     'Position', [1.14 0.12], 'FontSize', 12);

view(315, 50)
xlabel('Longitude unit: \circ', 'FontSize', 12)
ylabel('Latitude unit: \circ',   'FontSize', 12)
zlabel('Elevation (m)',          'FontSize', 12)
title('\bf(a) \rmUnfiltered',    'FontSize', 12)
grid on
box on
print(gcf, 'figure2a.png', '-r300', '-dpng');

% =========================================================================
% 7.  Figure (b): After one Gaussian filter pass
% =========================================================================
figure
surf(x, y, H01, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on
surf(x, y, H0, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
colormap(turbo)
colorbar
lighting gouraud
camlight headlight
material dull

set(gca, 'FontName', 'Arial', 'XTick', 0:0.25:1, 'XTickLabel', [], ...
               'YTick', 0:0.25:1, 'YTickLabel', [], 'TickDir', 'out')

% Labels (same positions as the first view)
text('String', '$$\varphi_2$$',       'Interpreter', 'latex', ...
     'Position', [0.18 1.27], 'FontSize', 12);
text('String', '$$\varphi_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 1.12], 'FontSize', 12);
text('String', '$$\varphi_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.88], 'FontSize', 12);
text('String', '$$\varphi_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.63], 'FontSize', 12);
text('String', '$$\varphi_1$$',       'Interpreter', 'latex', ...
     'Position', [0.19 0.28], 'FontSize', 12);

text('String', '$$\lambda_1$$',       'Interpreter', 'latex', ...
     'Position', [0.25 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.5  0.23], 'FontSize', 12);
text('String', '$$\lambda_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.75 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [1.0  0.23], 'FontSize', 12);
text('String', '$$\lambda_2$$',       'Interpreter', 'latex', ...
     'Position', [1.25 0.23], 'FontSize', 12);

view(315, 50)
xlabel('Longitude unit: \circ', 'FontSize', 12)
ylabel('Latitude unit: \circ',   'FontSize', 12)
zlabel('Elevation (m)',          'FontSize', 12)
title('\bf(b) \rmApply filtering once', 'FontSize', 12)
grid on
box on
print(gcf, 'figure2b.png', '-r300', '-dpng');

% =========================================================================
% 8.  Figure (c): After two Gaussian filter passes
% =========================================================================
figure
surf(x, y, H02, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on
surf(x, y, H0, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
colormap(turbo)
colorbar
lighting gouraud
camlight headlight
material dull

set(gca, 'FontName', 'Arial', 'XTick', 0:0.25:1, 'XTickLabel', [], ...
               'YTick', 0:0.25:1, 'YTickLabel', [], 'TickDir', 'out')

text('String', '$$\varphi_2$$',       'Interpreter', 'latex', ...
     'Position', [0.18 1.27], 'FontSize', 12);
text('String', '$$\varphi_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 1.12], 'FontSize', 12);
text('String', '$$\varphi_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.88], 'FontSize', 12);
text('String', '$$\varphi_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.63], 'FontSize', 12);
text('String', '$$\varphi_1$$',       'Interpreter', 'latex', ...
     'Position', [0.19 0.28], 'FontSize', 12);

text('String', '$$\lambda_1$$',       'Interpreter', 'latex', ...
     'Position', [0.25 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.5  0.23], 'FontSize', 12);
text('String', '$$\lambda_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.75 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [1.0  0.23], 'FontSize', 12);
text('String', '$$\lambda_2$$',       'Interpreter', 'latex', ...
     'Position', [1.25 0.23], 'FontSize', 12);

view(315, 50)
xlabel('Longitude unit: \circ', 'FontSize', 12)
ylabel('Latitude unit: \circ',   'FontSize', 12)
zlabel('Elevation (m)',          'FontSize', 12)
title('\bf(c) \rmApply filtering twice', 'FontSize', 12)
grid on
box on
print(gcf, 'figure2c.png', '-r300', '-dpng');

% =========================================================================
% 9.  Figure (d): After three Gaussian filter passes
% =========================================================================
figure
surf(x, y, H03, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on
surf(x, y, H0, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
colormap(turbo)
colorbar
lighting gouraud
camlight headlight
material dull

set(gca, 'FontName', 'Arial', 'XTick', 0:0.25:1, 'XTickLabel', [], ...
               'YTick', 0:0.25:1, 'YTickLabel', [], 'TickDir', 'out')

text('String', '$$\varphi_2$$',       'Interpreter', 'latex', ...
     'Position', [0.18 1.27], 'FontSize', 12);
text('String', '$$\varphi_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 1.12], 'FontSize', 12);
text('String', '$$\varphi_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.88], 'FontSize', 12);
text('String', '$$\varphi_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.04 0.63], 'FontSize', 12);
text('String', '$$\varphi_1$$',       'Interpreter', 'latex', ...
     'Position', [0.19 0.28], 'FontSize', 12);

text('String', '$$\lambda_1$$',       'Interpreter', 'latex', ...
     'Position', [0.25 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+15''$$', 'Interpreter', 'latex', ...
     'Position', [0.5  0.23], 'FontSize', 12);
text('String', '$$\lambda_1+30''$$', 'Interpreter', 'latex', ...
     'Position', [0.75 0.23], 'FontSize', 12);
text('String', '$$\lambda_1+45''$$', 'Interpreter', 'latex', ...
     'Position', [1.0  0.23], 'FontSize', 12);
text('String', '$$\lambda_2$$',       'Interpreter', 'latex', ...
     'Position', [1.25 0.23], 'FontSize', 12);

view(315, 50)
xlabel('Longitude unit: \circ', 'FontSize', 12)
ylabel('Latitude unit: \circ',   'FontSize', 12)
zlabel('Elevation (m)',          'FontSize', 12)
title('\bf(d) \rmApply filtering three times', 'FontSize', 12)
grid on
box on
print(gcf, 'figure2d.png', '-r300', '-dpng');