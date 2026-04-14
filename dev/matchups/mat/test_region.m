%% Test Overlap Region Finding (Approach B)
%
% Tests the intersection-centered line-fitting algorithm on one bow-tie
% boundary from the development latlon.mat data.
%
% Mirrors the Python notebook Test_region.ipynb.

clear; close all;

% Get function handles
fr = find_regions();

% Output directory for figures
fig_dir = '/home/xavier/Projects/overleaf/modis-correct';

%% Load data and reshape into bow-ties
data_file = fullfile('..', 'data', 'latlon.mat');
[lat, lon] = fr.load_latlon(data_file);
[lat_bt, lon_bt, n_skip] = fr.reshape_to_bowties(lat, lon);
n_bt = size(lat_bt, 1);
fprintf('%d complete bow-ties, skipped %d lines\n', n_bt, n_skip);
fprintf('Bow-tie array size: %s\n', mat2str(size(lat_bt)));

%% Pick a bow-tie boundary and demo detector pair
k = 901;       % 1-indexed bow-tie (Python k=900 is MATLAB k=901)
det_i = 10;    % trailing D10 from B_k
det_j = 1;     % leading D1 from B_{k+1}
w = 30;        % window width
q_max = 0.05;  % 5% of pixel spacing

% Extract scan lines (MATLAB indexing: lat_bt(bowtie, detector, pixel))
lat_meas = squeeze(lat_bt(k, det_i, :))';     % D10 of B_k
lon_meas = squeeze(lon_bt(k, det_i, :))';
lat_ref  = squeeze(lat_bt(k+1, det_j, :))';   % D1 of B_{k+1}
lon_ref  = squeeze(lon_bt(k+1, det_j, :))';

fprintf('Boundary k=%d: D%d(B_%d) vs D%d(B_%d)\n', k, det_i, k, det_j, k+1);
fprintf('Lat range: %.3f to %.3f\n', min(lat_ref), max(lat_ref));
fprintf('Lon range: %.3f to %.3f\n', min(lon_ref), max(lon_ref));

%% Step 1: Find the Scan-Line Intersection
NADIR = fr.NADIR_COL;
N_PIX = fr.N_PIXELS;
left_cols  = 1:NADIR;
right_cols = (NADIR+1):N_PIX;

delta_left  = lat_meas(left_cols) - lat_ref(left_cols);
delta_right = lat_meas(right_cols) - lat_ref(right_cols);

cross_left  = fr.find_intersection(lat_ref, lon_ref, lat_meas, lon_meas, left_cols);
cross_right = fr.find_intersection(lat_ref, lon_ref, lat_meas, lon_meas, right_cols);

if ~isnan(cross_left)
    fprintf('Left  crossing: P = %.1f\n', cross_left);
else
    fprintf('Left  crossing: None\n');
end
if ~isnan(cross_right)
    fprintf('Right crossing: P = %.1f\n', cross_right);
else
    fprintf('Right crossing: None\n');
end

% Figure: Delta_lat vs column
figure('Position', [100 100 1200 400]);
subplot(1,2,1);
plot(left_cols, delta_left * 111000, '.', 'MarkerSize', 2, 'Color', [0.27 0.51 0.71]);
hold on;
yline(0, 'r--', 'LineWidth', 0.8);
if ~isnan(cross_left)
    xline(cross_left, 'g-', sprintf('P=%.1f', cross_left), 'LineWidth', 1.5);
end
xlabel('Column P'); ylabel('\Delta lat (approx m)');
title(sprintf('D%d-D%d left side', det_i, det_j));

subplot(1,2,2);
plot(right_cols, delta_right * 111000, '.', 'MarkerSize', 2, 'Color', [0.27 0.51 0.71]);
hold on;
yline(0, 'r--', 'LineWidth', 0.8);
if ~isnan(cross_right)
    xline(cross_right, 'g-', sprintf('P=%.1f', cross_right), 'LineWidth', 1.5);
end
xlabel('Column P'); ylabel('\Delta lat (approx m)');
title(sprintf('D%d-D%d right side', det_i, det_j));
sgtitle(sprintf('Scan-line latitude difference at boundary k=%d', k));
exportgraphics(gcf, fullfile(fig_dir, 'region_intersection_mat.png'), 'Resolution', 150);

%% Figure: Scan-line crossing in lat/lon space (zoomed)
if ~isnan(cross_left)
    cx = round(cross_left);
    zoom_hw = 40;
    c0 = max(1, cx - zoom_hw);
    c1 = min(N_PIX, cx + zoom_hw);
    cols_zoom = c0:c1;

    figure('Position', [100 100 800 600]);
    plot(lon_ref(cols_zoom), lat_ref(cols_zoom), '.-', ...
         'MarkerSize', 6, 'Color', 'r', 'DisplayName', ...
         sprintf('D%d B_{k+1} (reference)', det_j));
    hold on;
    plot(lon_meas(cols_zoom), lat_meas(cols_zoom), '.-', ...
         'MarkerSize', 6, 'Color', 'b', 'DisplayName', ...
         sprintf('D%d B_k (measured)', det_i));

    % Mark crossing (interpolated)
    lat_cx = interp1(1:N_PIX, double(lat_ref), cross_left);
    lon_cx = interp1(1:N_PIX, double(lon_ref), cross_left);
    plot(lon_cx, lat_cx, 'p', 'MarkerSize', 15, 'Color', 'g', ...
         'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', ...
         'DisplayName', sprintf('Crossing P=%.1f', cross_left));

    xlabel('Longitude'); ylabel('Latitude');
    title(sprintf('D%d-D%d scan-line crossing (left side, k=%d)', det_i, det_j, k));
    legend('Location', 'best', 'FontSize', 9);
    exportgraphics(gcf, fullfile(fig_dir, 'region_scanline_crossing_mat.png'), 'Resolution', 150);
end

%% Step 2: Fit Reference Line and Compute Quality
result_left = fr.find_region_one_side(lat_ref, lon_ref, lat_meas, lon_meas, ...
                                       left_cols, w, q_max);
if ~isempty(result_left)
    fprintf('\nOverlap region: P = %d to %d\n', result_left.cols(1), result_left.cols(end));
    fprintf('  Width: %d columns\n', length(result_left.cols));
    fprintf('  Crossing at P = %.1f\n', result_left.cross_col);
    fprintf('  q range: [%.5f, %.5f]\n', min(result_left.q), max(result_left.q));
    fprintf('  |d_perp| range: [%.1f, %.1f] m\n', ...
            min(abs(result_left.d_perp)), max(abs(result_left.d_perp)));
    fprintf('  pixel spacing range: [%.0f, %.0f] m\n', ...
            min(result_left.pixel_spacing), max(result_left.pixel_spacing));
else
    fprintf('No overlap region found on left side.\n');
end

%% Figure: Fitted line on projected scan lines (zoomed)
if ~isempty(result_left)
    line = result_left.line;
    win_cols = result_left.win_cols;

    [x_ref, y_ref] = fr.equirect_project(double(lat_ref(win_cols)), ...
                                          double(lon_ref(win_cols)), ...
                                          line.lat0, line.lon0);
    [x_meas, y_meas] = fr.equirect_project(double(lat_meas(win_cols)), ...
                                            double(lon_meas(win_cols)), ...
                                            line.lat0, line.lon0);
    x_fit = linspace(min(x_ref), max(x_ref), 100);
    y_fit = line.slope * x_fit + line.intercept;

    figure('Position', [100 100 1000 500]);
    plot(x_ref, y_ref, 'ro', 'MarkerSize', 5, 'DisplayName', ...
         sprintf('D%d B_{k+1} (reference)', det_j));
    hold on;
    plot(x_meas, y_meas, 'bs', 'MarkerSize', 5, 'DisplayName', ...
         sprintf('D%d B_k (measured)', det_i));
    plot(x_fit, y_fit, '-', 'Color', [0.6 0 0], 'LineWidth', 2, ...
         'DisplayName', 'Fitted line');

    % Mark overlap columns
    overlap_mask = ismember(win_cols, result_left.cols);
    plot(x_meas(overlap_mask), y_meas(overlap_mask), 'gs', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.3 0.8 0.3], ...
         'DisplayName', sprintf('Overlap (%d cols)', length(result_left.cols)));

    xlabel('x (east, m)'); ylabel('y (north, m)');
    title(sprintf('D%d-D%d projected scan lines with fitted line (left side, k=%d)', ...
                  det_i, det_j, k));
    legend('Location', 'best', 'FontSize', 8);
    axis equal;
    exportgraphics(gcf, fullfile(fig_dir, 'region_fitted_line_mat.png'), 'Resolution', 150);
end

%% Figure: Quality profile
if ~isempty(result_left)
    figure('Position', [100 100 1200 450]);

    % Left: quality metric q(P)
    subplot(1,2,1);
    plot(result_left.win_cols, result_left.q_full, '.-', ...
         'MarkerSize', 6, 'Color', [0.27 0.51 0.71]);
    hold on;
    yline(q_max, 'r--', sprintf('q_{max} = %.2f', q_max), 'LineWidth', 1);
    % Shade overlap
    x_shade = [result_left.cols(1), result_left.cols(end), ...
               result_left.cols(end), result_left.cols(1)];
    yl = ylim;
    y_shade = [yl(1), yl(1), yl(2), yl(2)];
    patch(x_shade, y_shade, 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
          'DisplayName', sprintf('Overlap (%d cols)', length(result_left.cols)));
    xline(result_left.cross_col, ':', 'Color', [1 0.5 0], 'LineWidth', 1.5, ...
          'DisplayName', 'Crossing');
    xlabel('Column P'); ylabel('q = |d_{perp}| / pixel spacing');
    title(sprintf('D%d-D%d quality profile (left side)', det_i, det_j));
    legend('Location', 'best', 'FontSize', 8);
    ylim([0, min(0.15, max(result_left.q_full) * 1.2)]);

    % Right: perpendicular distance
    subplot(1,2,2);
    plot(result_left.win_cols, result_left.d_perp_full, '.-', ...
         'MarkerSize', 6, 'Color', [0.27 0.51 0.71]);
    hold on;
    yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    threshold = q_max * result_left.ps_full;
    fill([result_left.win_cols, fliplr(result_left.win_cols)], ...
         [threshold, fliplr(-threshold)], ...
         'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    xline(result_left.cross_col, ':', 'Color', [1 0.5 0], 'LineWidth', 1.5);
    xlabel('Column P'); ylabel('Perpendicular distance (m)');
    title(sprintf('D%d-D%d perpendicular distance to fitted line', det_i, det_j));
    exportgraphics(gcf, fullfile(fig_dir, 'region_quality_profile_mat.png'), 'Resolution', 150);
end

%% Step 3: Run All Detector Pairs
% lat_bt/lon_bt are (n_bt, 10, 1354); extract (10, 1354) for each bow-tie
lat_k  = squeeze(lat_bt(k, :, :));
lon_k  = squeeze(lon_bt(k, :, :));
lat_k1 = squeeze(lat_bt(k+1, :, :));
lon_k1 = squeeze(lon_bt(k+1, :, :));

regions = fr.find_regions_one_boundary(lat_k, lon_k, lat_k1, lon_k1, w, q_max);

fprintf('\nFound %d overlap regions at boundary k=%d\n', length(regions), k);
fprintf('%8s %6s %6s %6s %6s %8s %8s %8s\n', ...
        'Pair', 'Side', 'P_min', 'P_max', 'Width', 'Cross_P', 'q_min', 'q_max');
fprintf('%s\n', repmat('-', 1, 64));
for i = 1:length(regions)
    r = regions(i);
    fprintf('D%d-D%d  %6s %6d %6d %6d %8.1f %8.5f %8.5f\n', ...
            r.det_i, r.det_j, r.side, ...
            r.cols(1), r.cols(end), length(r.cols), ...
            r.cross_col, min(r.q), max(r.q));
end

%% Step 4: Summary Map of All Overlap Regions
figure('Position', [100 100 1200 600]);

% Background: trailing detectors (D7-10) of B_k in blue
for d = 7:10
    plot(squeeze(lon_bt(k, d, :)), squeeze(lat_bt(k, d, :)), ...
         '.', 'MarkerSize', 0.5, 'Color', [0 0 1 0.2]);
    hold on;
end
% Leading detectors (D1-3) of B_{k+1} in red
for d = 1:3
    plot(squeeze(lon_bt(k+1, d, :)), squeeze(lat_bt(k+1, d, :)), ...
         '.', 'MarkerSize', 0.5, 'Color', [1 0 0 0.2]);
end

% Color map for detector pairs
cmap = lines(12);
pairs = fr.DETECTOR_PAIRS;

for i = 1:length(regions)
    r = regions(i);
    % Find color index for this pair
    pair_idx = find(pairs(:,1) == r.det_i & pairs(:,2) == r.det_j, 1);
    c = cmap(pair_idx, :);
    cols = r.cols;

    % Plot overlap pixels
    plot(squeeze(lon_bt(k+1, r.det_j, cols)), ...
         squeeze(lat_bt(k+1, r.det_j, cols)), ...
         'o', 'MarkerSize', 3, 'Color', c, 'MarkerFaceColor', c);
    plot(squeeze(lon_bt(k, r.det_i, cols)), ...
         squeeze(lat_bt(k, r.det_i, cols)), ...
         's', 'MarkerSize', 3, 'Color', c);

    % Mark crossing point
    cx = round(r.cross_col);
    plot(lon_bt(k+1, r.det_j, cx), lat_bt(k+1, r.det_j, cx), ...
         'p', 'MarkerSize', 10, 'Color', c, ...
         'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
end

xlabel('Longitude'); ylabel('Latitude');
title(sprintf('Overlap regions at boundary k=%d (q_{max}=%.2f, w=%d) - %d regions', ...
              k, q_max, w, length(regions)));
exportgraphics(gcf, fullfile(fig_dir, 'region_map_all_mat.png'), 'Resolution', 150);

%% Validation: Haversine distances at overlap regions
figure('Position', [100 100 1200 800]);
n_show = min(4, length(regions));
for idx = 1:n_show
    subplot(2, 2, idx);
    r = regions(idx);
    cols = r.cols;
    wcols = r.win_cols;

    % Haversine over full window
    hav_full = fr.haversine(lat_bt(k, r.det_i, wcols), lon_bt(k, r.det_i, wcols), ...
                            lat_bt(k+1, r.det_j, wcols), lon_bt(k+1, r.det_j, wcols));
    hav_full = squeeze(hav_full);

    % Haversine at overlap columns
    hav_overlap = fr.haversine(lat_bt(k, r.det_i, cols), lon_bt(k, r.det_i, cols), ...
                               lat_bt(k+1, r.det_j, cols), lon_bt(k+1, r.det_j, cols));
    hav_overlap = squeeze(hav_overlap);

    plot(wcols, hav_full, '.-', 'MarkerSize', 4, 'Color', [0.7 0.7 0.7]);
    hold on;
    plot(cols, hav_overlap, '.-', 'MarkerSize', 6, 'Color', [0.27 0.51 0.71]);
    xline(r.cross_col, ':', 'Color', [1 0.5 0], 'LineWidth', 1.5);
    xlabel('Column P'); ylabel('Haversine distance (m)');
    title(sprintf('D%d-D%d %s (%d cols)', r.det_i, r.det_j, r.side, length(cols)));
    legend('Full window', 'Overlap region', 'Crossing', 'Location', 'best', 'FontSize', 8);
end
sgtitle(sprintf('Haversine validation at boundary k=%d', k));
exportgraphics(gcf, fullfile(fig_dir, 'region_haversine_validation_mat.png'), 'Resolution', 150);

fprintf('\nDone. Figures saved to %s\n', fig_dir);
