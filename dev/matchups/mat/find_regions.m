function fr = find_regions()
% FIND_REGIONS  Return a struct of function handles for overlap region finding.
%
%   Implements Approach B from find_regions_plan.tex:
%     Intersection-centered line-fitting with perpendicular distance.
%
%   Algorithm for each detector pair on each side of nadir:
%     1. Find the column where two scan lines cross (delta_lat sign change)
%     2. Center a w=30 pixel window on the crossing point
%     3. Fit a reference line through the leading scan (OLS, equirectangular)
%     4. Compute perpendicular distance from trailing scan to fitted line
%     5. Normalize: q = |d_perp| / pixel_spacing
%     6. Expand from best-overlap column to define region where q < q_max
%
%   Usage:
%     fr = find_regions();
%     cross_col = fr.find_intersection(lat_ref, lon_ref, lat_meas, lon_meas, cols);
%     regions   = fr.find_regions_one_boundary(lat_k, lon_k, lat_k1, lon_k1, w, q_max);

    fr.haversine               = @haversine;
    fr.load_latlon             = @load_latlon;
    fr.reshape_to_bowties      = @reshape_to_bowties;
    fr.equirect_project        = @equirect_project;
    fr.find_intersection       = @find_intersection;
    fr.fit_line                = @fit_line;
    fr.perp_distance           = @perp_distance;
    fr.compute_pixel_spacing   = @compute_pixel_spacing;
    fr.find_region_one_side    = @find_region_one_side;
    fr.find_regions_one_boundary = @find_regions_one_boundary;

    % Constants
    fr.N_DETECTORS = 10;
    fr.N_PIXELS    = 1354;
    fr.FIRST_DETECTOR = 6;
    fr.R_EARTH     = 6371000.0;
    fr.NADIR_COL   = 677;  % 1-indexed
    % Detector pairs: trailing D from B_k, leading D from B_{k+1} (1-indexed)
    pairs = [];
    for di = 7:10
        for dj = 1:3
            pairs = [pairs; di, dj]; %#ok<AGROW>
        end
    end
    fr.DETECTOR_PAIRS = pairs;
end


% ============================================================
% Core: Haversine great-circle distance
% ============================================================
function d = haversine(lat1, lon1, lat2, lon2)
% HAVERSINE  Vectorized great-circle distance in meters.
%   All inputs in degrees; arrays must be same size or scalar.
    R = 6371000.0;
    lat1_r = deg2rad(lat1);
    lat2_r = deg2rad(lat2);
    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    a = sin(dlat/2).^2 + cos(lat1_r).*cos(lat2_r).*sin(dlon/2).^2;
    d = 2 * R * asin(sqrt(a));
end


% ============================================================
% I/O: Load lat/lon from .mat file
% ============================================================
function [lat, lon] = load_latlon(filepath)
% LOAD_LATLON  Load latitude and longitude from a .mat file.
%   Returns arrays of size (N_lines, 1354).
    data = load(filepath, 'lat', 'lon');
    % .mat file stores (1354, N_lines); transpose to (N_lines, 1354)
    lat = single(data.lat');
    lon = single(data.lon');
end


% ============================================================
% Reshape scan lines into complete bow-ties
% ============================================================
function [lat_bt, lon_bt, n_skip] = reshape_to_bowties(lat, lon, first_detector)
% RESHAPE_TO_BOWTIES  Group scan lines into complete 10-detector bow-ties.
%   Discards partial first and trailing bow-ties.
%
%   Returns lat_bt, lon_bt of size (N_bt, 10, 1354) and n_skip.
    if nargin < 3, first_detector = 6; end
    N_DET = 10;
    N_PIX = 1354;
    n_lines = size(lat, 1);

    % Skip partial first bow-tie: D=first_detector through D=10
    n_skip = N_DET - (first_detector - 1);
    n_remaining = n_lines - n_skip;
    n_bt = floor(n_remaining / N_DET);
    last = n_skip + n_bt * N_DET;

    lat_bt = reshape(lat(n_skip+1:last, :), [N_DET, n_bt, N_PIX]);
    lat_bt = permute(lat_bt, [2 1 3]);  % (n_bt, 10, 1354)
    lon_bt = reshape(lon(n_skip+1:last, :), [N_DET, n_bt, N_PIX]);
    lon_bt = permute(lon_bt, [2 1 3]);
end


% ============================================================
% Projection
% ============================================================
function [x, y] = equirect_project(lat, lon, lat0, lon0)
% EQUIRECT_PROJECT  Local equirectangular projection (x_east, y_north) in meters.
    R = 6371000.0;
    x = deg2rad(lon - lon0) .* cos(deg2rad(lat0)) * R;
    y = deg2rad(lat - lat0) * R;
end


% ============================================================
% Step 1: Find intersection
% ============================================================
function cross_col = find_intersection(lat_ref, lon_ref, lat_meas, lon_meas, cols)
% FIND_INTERSECTION  Find the column where two scan lines cross.
%   Uses latitude difference as proxy for along-track displacement.
%
%   cols: 1-indexed column indices to search.
%   Returns interpolated 1-indexed crossing column, or NaN if none found.
    delta = lat_meas(cols) - lat_ref(cols);
    signs = sign(delta);

    % Find sign changes
    d_signs = diff(signs);
    idx = find(d_signs ~= 0, 1, 'first');
    if isempty(idx)
        cross_col = NaN;
        return;
    end

    c0 = cols(idx);
    c1 = cols(idx + 1);
    d0 = double(delta(idx));
    d1 = double(delta(idx + 1));

    % Linear interpolation
    if d1 == d0
        cross_col = c0;
    else
        frac = -d0 / (d1 - d0);
        cross_col = c0 + frac * (c1 - c0);
    end
end


% ============================================================
% Step 2: Fit reference line
% ============================================================
function line = fit_line(lat_ref, lon_ref, center_col, w)
% FIT_LINE  Fit a straight line through w pixels of the reference scan.
%   center_col: 1-indexed center column.
%   Returns struct with slope, intercept, lat0, lon0, cols, x_ref, y_ref.
    if nargin < 4, w = 30; end
    N_PIX = 1354;
    half_w = floor(w / 2);
    col_start = max(1, center_col - half_w);
    col_end   = min(N_PIX, center_col + half_w - 1);
    cols = col_start:col_end;

    % Projection center at window midpoint
    mid = floor(length(cols) / 2) + 1;
    lat0 = double(lat_ref(cols(mid)));
    lon0 = double(lon_ref(cols(mid)));
    [x, y] = equirect_project(double(lat_ref(cols)), double(lon_ref(cols)), lat0, lon0);

    % OLS line fit: y = slope * x + intercept
    p = polyfit(x, y, 1);
    line.slope     = p(1);
    line.intercept = p(2);
    line.lat0      = lat0;
    line.lon0      = lon0;
    line.cols      = cols;
    line.x_ref     = x;
    line.y_ref     = y;
end


% ============================================================
% Step 3: Perpendicular distance
% ============================================================
function d_perp = perp_distance(lat, lon, line)
% PERP_DISTANCE  Signed perpendicular distance (m) from points to fitted line.
    cols = line.cols;
    [x, y] = equirect_project(double(lat(cols)), double(lon(cols)), line.lat0, line.lon0);
    m = line.slope;
    b = line.intercept;
    % Line: m*x - y + b = 0
    d_perp = (m * x - y + b) / sqrt(m^2 + 1);
end


% ============================================================
% Step 4: Pixel spacing
% ============================================================
function ps = compute_pixel_spacing(lat, lon, cols)
% COMPUTE_PIXEL_SPACING  Along-scan pixel spacing (m) via Haversine.
%   Centered differences in interior; one-sided at edges.
    sp = haversine(lat(cols(1:end-1)), lon(cols(1:end-1)), ...
                   lat(cols(2:end)),   lon(cols(2:end)));
    n = length(cols);
    ps = zeros(1, n);
    ps(1)       = sp(1);
    ps(end)     = sp(end);
    ps(2:end-1) = 0.5 * (sp(1:end-1) + sp(2:end));
end


% ============================================================
% Steps 5-6: Quality metric and overlap region
% ============================================================
function result = find_region_one_side(lat_ref, lon_ref, lat_meas, lon_meas, ...
                                       cols, w, q_max)
% FIND_REGION_ONE_SIDE  Find overlap region on one side of nadir.
%   cols: 1-indexed columns for one side.
%   Returns struct with region info, or empty [] if no valid overlap.
    if nargin < 6, w = 30; end
    if nargin < 7, q_max = 0.05; end

    % Step 1: intersection
    cross_col = find_intersection(lat_ref, lon_ref, lat_meas, lon_meas, cols);
    if isnan(cross_col)
        result = [];
        return;
    end
    center = round(cross_col);

    % Step 2: fit reference line
    line = fit_line(lat_ref, lon_ref, center, w);

    % Step 3: perpendicular distances
    d_perp = perp_distance(lat_meas, lon_meas, line);

    % Step 4: pixel spacing and quality metric
    win_cols = line.cols;
    ps = compute_pixel_spacing(lat_ref, lon_ref, win_cols);
    q = abs(d_perp) ./ ps;

    % Step 5: contiguous region expanding from best-overlap column
    [~, best_idx] = min(q);
    mask = q < q_max;

    if ~mask(best_idx)
        result = [];
        return;
    end

    % Expand outward from best-overlap column
    left_idx = best_idx;
    while left_idx > 1 && mask(left_idx - 1)
        left_idx = left_idx - 1;
    end
    right_idx = best_idx;
    while right_idx < length(mask) && mask(right_idx + 1)
        right_idx = right_idx + 1;
    end

    overlap_cols = win_cols(left_idx:right_idx);

    result.cols          = overlap_cols;        % 1-indexed columns in overlap
    result.q             = q(left_idx:right_idx);
    result.d_perp        = d_perp(left_idx:right_idx);
    result.pixel_spacing = ps(left_idx:right_idx);
    result.cross_col     = cross_col;           % interpolated crossing (1-indexed)
    result.center_col    = center;
    result.line          = line;
    % Full window diagnostics
    result.q_full        = q;
    result.d_perp_full   = d_perp;
    result.ps_full       = ps;
    result.win_cols      = win_cols;
end


% ============================================================
% Main entry point
% ============================================================
function regions = find_regions_one_boundary(lat_k, lon_k, lat_k1, lon_k1, w, q_max)
% FIND_REGIONS_ONE_BOUNDARY  Find all overlap regions at one bow-tie boundary.
%   lat_k, lon_k:   (10, 1354) arrays for bow-tie B_k
%   lat_k1, lon_k1: (10, 1354) arrays for bow-tie B_{k+1}
%
%   Returns a struct array of regions.
    if nargin < 5, w = 30; end
    if nargin < 6, q_max = 0.05; end

    NADIR = 677;
    N_PIX = 1354;
    left_cols  = 1:NADIR;           % 1-indexed
    right_cols = (NADIR+1):N_PIX;

    % Detector pairs
    pairs = [];
    for di = 7:10
        for dj = 1:3
            pairs = [pairs; di, dj]; %#ok<AGROW>
        end
    end

    regions = [];
    count = 0;
    for p = 1:size(pairs, 1)
        det_i = pairs(p, 1);  % trailing D from B_k
        det_j = pairs(p, 2);  % leading D from B_{k+1}

        % Trailing = measured, leading = reference
        lat_meas = squeeze(lat_k(det_i, :));
        lon_meas = squeeze(lon_k(det_i, :));
        lat_ref  = squeeze(lat_k1(det_j, :));
        lon_ref  = squeeze(lon_k1(det_j, :));

        sides = {'left', 'right'};
        col_ranges = {left_cols, right_cols};
        for s = 1:2
            result = find_region_one_side(lat_ref, lon_ref, lat_meas, lon_meas, ...
                                          col_ranges{s}, w, q_max);
            if ~isempty(result)
                count = count + 1;
                result.det_i = det_i;
                result.det_j = det_j;
                result.side  = sides{s};
                if count == 1
                    regions = result;
                else
                    regions(count) = result; %#ok<AGROW>
                end
            end
        end
    end

    if count == 0
        regions = [];
    end
end
