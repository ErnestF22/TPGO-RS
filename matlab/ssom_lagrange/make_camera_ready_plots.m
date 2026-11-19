function make_camera_ready_plots(matFiles, outDir, doMatrix)
%MAKE_CAMERA_READY_PLOTS  Redo the TPGO result plots from saved workspaces.
%
%   make_camera_ready_plots()                      % all _n*_mindeg*.mat in cwd
%   make_camera_ready_plots("_n5_mindeg2.mat")     % one file
%   make_camera_ready_plots({"a.mat","b.mat"}, "figs")
%   make_camera_ready_plots(files, "figs", true/false)   % force matrix plots on/off
%
%   The RS usefulness matrices (plot_rs_useful_matrix) are built once, after
%   all files of the call are done. Default: on for several files, off for one.
%   If you loop over files yourself, call plot_rs_useful_matrix once afterwards.
%
%   Each .mat must contain `results` and `sigmas` (as saved after your
%   plotting script). The file name stem (e.g. "_n5_mindeg2") is appended to
%   the output names, exactly like your old test_str.
%
%   Requires R2020a+ (exportgraphics).

if nargin < 1 || isempty(matFiles)
    d = dir('_n*_mindeg*.mat');
    matFiles = fullfile({d.folder}, {d.name});
end
cfg.connect   = false;    % dots only, no connecting lines
cfg.legendFontSize = 10;  % slightly larger legend text

if nargin < 1 || isempty(matFiles)
    d = dir('_n*_mindeg*.mat');
    matFiles = fullfile({d.folder}, {d.name});
end

if nargin < 2 || isempty(outDir)
    outDir = 'lsom_results';
end

% RS usefulness matrices need the CSVs of all runs, so by default they are
% made only when several files are processed in this one call (i.e. at its end).
% With a single file (e.g. when you call this in your own loop) they are skipped.
if nargin < 3 || isempty(doMatrix), doMatrix = numel(cellstr(matFiles)) > 1; end

if isstring(outDir)
    outDir = char(outDir);
end

if ~(startsWith(outDir, filesep) || (ispc && numel(outDir) >= 2 && outDir(2) == ':'))
    outDir = fullfile(pwd, outDir);
end

if ~isfolder(outDir)
    mkdir(outDir);
end

cfg.connect = false;    % dots only, no connecting lines
cfg.legendFontSize = 10;  % slightly larger legend text


%% ---------------- style configuration (edit here) ----------------------
cfg.widthCm   = 8.8;      % single column ~8.5-9 cm; double column ~17-18 cm
cfg.heightCm  = 6.2;
cfg.fontName  = 'Times';
cfg.fontSize  = 8;
cfg.lineWidth = 0.75;     % data lines
cfg.axWidth   = 0.6;      % axes/box line width
cfg.markerSz  = 5;
cfg.connect   = false;     % thin line connecting the markers
cfg.formats   = {'pdf', 'png'};   % add 'eps' if the venue wants it
cfg.legendLoc = 'best';
cfg.xSpacing  = 'even';   % 'even': sigmas equally spaced (no overlap at 0, .01, .02)
                          % 'linear': true numeric spacing of sigma
cfg.legendFontSize = 10;

% Colour-blind safe (Okabe-Ito) + distinct markers, so it also works in
% greyscale print. Order: ICP, Procrustes, Procrustes-QP, RS.
methods = struct( ...
    'prefix', {'manopt_sep', 'procrustes', 'procrustes_qp', 'ssom'}, ...
    'label',  {'TPGO-ICP', 'TPGO-PROCRUSTES', 'TPGO-PROCRUSTES-QP', 'TPGO-RS'}, ...
    'marker', {'o', 's', 'd', '^'}, ...
    'color',  {[0.835 0.369 0.000], [0.000 0.447 0.698], ...
    [0.800 0.475 0.655], [0.000 0.619 0.451]});

% suffix in `results`, output base name, y label, log-y?
metrics = struct( ...
    'suffix', {'rot_errs', 'transl_errs', 'exec_times', ...
    'scale_errs_max', 'scale_errs_mean'}, ...
    'base',   {'rot_errors', 'transl_errors', 'exec_times', ...
    'max_scale_err', 'mean_scale_err'}, ...
    'ylabel', {'Rotation error [$^\circ$]', 'Translation error [m]', ...
    'Execution time [s]', 'Max scale error', ...
    'Geometric mean scale error'}, ...
    'logy',   {false, false, false, false, false});

% workspace variables to dump as CSV (one file per variable per .mat)
csvVars = {'rs_success_bools', 'rot_dets_ok', ...
    'lambdas_acceptable', 'rs_actually_useful'};

%% ---------------- main loop --------------------------------------------
for k = 1:numel(matFiles)
    S = load(matFiles{k});
    [~, stem] = fileparts(matFiles{k});      % e.g. "_n5_mindeg2"
    sigmas = S.sigmas(:);

    if strcmp(cfg.xSpacing, 'even')
        xv = (1:numel(sigmas))';
    else
        xv = sigmas;
    end
    if k == 1
        sigmasRef = sigmas;
    elseif ~isequal(sigmasRef, sigmas)
        warning('sigmas in %s differ from the first file; matrix plot labels use the first.', ...
        matFiles{k});
    end

    % ---- per-sigma rs_usefulness_structs CSVs ---------------------------
    if isfield(S, 'rs_usefulness_structs_mm')
        sigmaTags = cellstr(compose("%03d", round(sigmas * 100)));
        for nn = 1:numel(sigmas)
            T = struct2table(S.rs_usefulness_structs_mm(nn, :));
            writetable(T, fullfile(outDir, ['rs_usefulness_structs' stem '_sigma' sigmaTags{nn} '.csv']));
        end
    else
        warning('rs_usefulness_structs_mm not found in %s, skipping.', matFiles{k});
    end

    % ---- CSVs (tab-delimited, same names as before + outDir) ------------
    for v = csvVars
        name = v{1};
        if ~isfield(S, name)
            warning('Variable %s not found in %s, skipping CSV.', name, matFiles{k});
            continue
        end
        writematrix(double(S.(name)), fullfile(outDir, [name stem '.csv']), 'Delimiter', 'tab');
    end

    for m = metrics
        fig = figure('Name', m.base, 'Color', 'w', 'Units', 'centimeters', ...
            'Position', [2 2 cfg.widthCm cfg.heightCm]);
        ax = axes(fig); hold(ax, 'on');

        for i = 1:numel(methods)
            fld = [methods(i).prefix '_' m.suffix];
            if ~isfield(S.results, fld)
                warning('Missing field %s in %s, skipping.', fld, matFiles{k});
                continue
            end
            y = S.results.(fld);
            ls = '-'; if ~cfg.connect, ls = 'none'; end
            plot(ax, xv, y(:), ...
                'LineStyle', ls, 'LineWidth', cfg.lineWidth, ...
                'Color', methods(i).color, ...
                'Marker', methods(i).marker, 'MarkerSize', cfg.markerSz, ...
                'MarkerFaceColor', methods(i).color, ...
                'MarkerEdgeColor', methods(i).color, ...
                'DisplayName', methods(i).label);
        end

        % axes cosmetics
        box(ax, 'on'); grid(ax, 'on');
        ax.GridAlpha = 0.15; ax.GridLineStyle = '-';
        ax.FontName = cfg.fontName; ax.FontSize = cfg.fontSize;
        ax.LineWidth = cfg.axWidth;
        ax.TickLabelInterpreter = 'latex';
        ax.TickDir = 'in';
        if m.logy, ax.YScale = 'log'; end

        % pad = 0.05 * (max(sigmas) - min(sigmas));
        % if pad == 0, pad = 0.1; end
        % xlim(ax, [min(sigmas) - pad, max(sigmas) + pad]);
        % if numel(sigmas) <= 10, xticks(ax, sigmas); end
        if strcmp(cfg.xSpacing, 'even')
            % sigmas at equally spaced positions, labelled with their values
            xlim(ax, [0.5, numel(sigmas) + 0.5]);
            lab = compose('%g', sigmas);
            if numel(sigmas) > 10, lab(2:2:end) = ""; end   % thin labels if crowded
            xticks(ax, xv);
            xticklabels(ax, lab);
        else
            pad = 0.05 * (max(sigmas) - min(sigmas));
            if pad == 0, pad = 0.1; end
            xlim(ax, [min(sigmas) - pad, max(sigmas) + pad]);
        end

        xlabel(ax, '$\sigma$', 'Interpreter', 'latex', 'FontSize', cfg.fontSize + 1);
        ylabel(ax, m.ylabel,  'Interpreter', 'latex', 'FontSize', cfg.fontSize + 1);

        lgd = legend(ax, 'Location', cfg.legendLoc, 'Interpreter', 'latex');
        lgd.FontSize = cfg.fontSize - 1;
        lgd.Box = 'off';

        % export (vector PDF + hi-res PNG)
        base = fullfile(outDir, [m.base stem]);
        for f = cfg.formats
            ext = f{1};
            file = [base '.' ext];
            if strcmp(ext, 'png')
                exportgraphics(fig, file, 'Resolution', 600);
            else
                exportgraphics(fig, file, 'ContentType', 'vector');
            end
        end
    close(fig);
    end
    close all;
end

% red/green usefulness matrices (one figure per mindeg), made once, after all
% files have been processed and all CSVs written; they read the "ok" column of
% rs_usefulness_structs_n<N>_mindeg<D>_sigma<S>.csv
if doMatrix
    plot_rs_useful_matrix(outDir, outDir, sigmasRef, cfg);
end

end %file function


function plot_rs_useful_matrix(csvDir, outDir, sigmas, cfg)
%PLOT_RS_USEFUL_MATRIX  Red/green matrix of the RS usefulness "ok" flag, one figure per mindeg.
%
%   plot_rs_useful_matrix(csvDir, outDir, sigmas, cfg)
%
%   Reads csvDir/rs_usefulness_structs_n<N>_mindeg<D>_sigma<S>.csv, where <S> is
%   round(100*sigma) written with 3 digits. Each file has one row per test and
%   the columns at_least_one_step, ok, rec, glob; only the "ok" column is used.
%   For every mindeg found, one figure is produced:
%     x axis : n  (one box per n)
%     y axis : sigma (one row of boxes per noise level)
%     inside each box: one vertical tick per test, green = ok 1, red = ok 0
%
%   sigmas : optional vector with the noise levels. If empty (or if its length
%            does not match the files found), the y labels are taken from the
%            file names (<S>/100).
%   cfg    : optional style struct (same fields as in make_camera_ready_plots).

if nargin < 1 || isempty(csvDir), csvDir = 'figures'; end
if nargin < 2 || isempty(outDir), outDir = csvDir; end
if nargin < 3, sigmas = []; end
if nargin < 4 || isempty(cfg)
    cfg.widthCm = 8.8; cfg.heightCm = 6.2;
    cfg.fontName = 'Times'; cfg.fontSize = 8;
    cfg.axWidth = 0.6; cfg.formats = {'pdf', 'png'};
end
if ~isfield(cfg, 'tickFill'), cfg.tickFill = 0.8; end  % fraction of a tick slot that is coloured
if ~isfolder(outDir), mkdir(outDir); end

% colours: slightly different lightness so it still reads in greyscale
colRed   = [0.75 0.15 0.15];
colGreen = [0.55 0.82 0.45];

%% find and parse the CSV files
d = dir(fullfile(csvDir, 'rs_usefulness_structs_n*_mindeg*_sigma*.csv'));
tok = regexp({d.name}, ...
    '^rs_usefulness_structs_n(\d+)_mindeg(\d+)_sigma(\d+)\.csv$', 'tokens', 'once');
keep = ~cellfun(@isempty, tok);
d = d(keep); tok = tok(keep);
if isempty(d)
    warning('No rs_usefulness_structs_n*_mindeg*_sigma*.csv files found in %s.', csvDir);
    return
end
nVals  = cellfun(@(t) str2double(t{1}), tok);
mdVals = cellfun(@(t) str2double(t{2}), tok);
sgVals = cellfun(@(t) str2double(t{3}), tok);    % round(100*sigma)

%% one figure per mindeg
for md = unique(mdVals)
    inMd = find(mdVals == md);
    ns = unique(nVals(inMd));
    sg = unique(sgVals(inMd));
    nN = numel(ns);
    nS = numel(sg);

    data = cell(nS, nN);                 % data{r, j}: ok flags of sigma r, n j
    for i = inMd
        j = find(ns == nVals(i));
        r = find(sg == sgVals(i));
        ok = readOkColumn(fullfile(d(i).folder, d(i).name));
        data{r, j} = ok;
        fprintf('mindeg %d, n=%d, sigma tag %03d: %d tests, %.1f%% ok, %d NaN\n', ...
            md, ns(j), sg(r), numel(ok), 100 * mean(ok(~isnan(ok)) ~= 0), nnz(isnan(ok)));
    end
    T = max(cellfun(@numel, data(:)));   % tests per box
    if T == 0
        warning('mindeg %d: no usable "ok" data, skipping.', md);
        continue
    end

    if ~isempty(sigmas) && numel(sigmas) == nS
        yLab = compose('%g', sigmas(:)');
    else
        if ~isempty(sigmas)
            warning('mindeg %d: %d sigma levels in files but %d sigmas given; labels taken from file names.', ...
                md, nS, numel(sigmas));
        end
        yLab = compose('%g', sg(:)' / 100);
    end

    % tick centres, split by value
    XC = []; RC = []; V = [];
    for j = 1:nN
        for r = 1:nS
            ok = data{r, j};
            idx = find(~isnan(ok));
            if isempty(idx), continue; end
            XC = [XC; (j - 1) * T + idx(:)];     %#ok<AGROW>
            RC = [RC; r * ones(numel(idx), 1)];  %#ok<AGROW>
            V  = [V;  ok(idx) ~= 0];             %#ok<AGROW>
        end
    end

    h = max(cfg.heightCm, 1.8 + 0.45 * nS);
    fig = figure('Name', sprintf('rs useful, mindeg %d', md), 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 cfg.widthCm h]);
    ax = axes(fig); hold(ax, 'on');

    % real data (not shown in the legend)
    drawCells(ax, XC(V == 0), RC(V == 0), colRed,   cfg.tickFill);
    drawCells(ax, XC(V == 1), RC(V == 1), colGreen, cfg.tickFill);

    % box outlines, one per (n, sigma) that has data
    for j = 1:nN
        for r = 1:nS
            if isempty(data{r, j}), continue; end
            rectangle(ax, 'Position', [(j - 1) * T + 0.5, r - 0.5, T, 1], ...
                'EdgeColor', 'k', 'LineWidth', cfg.axWidth);
        end
    end

    % dummy patches for the legend
    hR = patch(ax, NaN, NaN, colRed,   'EdgeColor', 'none', 'DisplayName', 'ok = 0');
    hG = patch(ax, NaN, NaN, colGreen, 'EdgeColor', 'none', 'DisplayName', 'ok = 1');

    % axes
    xlim(ax, [0.5, nN * T + 0.5]);
    ylim(ax, [0.5, nS + 0.5]);
    xticks(ax, ((0:nN - 1) * T + (T + 1) / 2));
    xticklabels(ax, compose('%d', ns));
    yticks(ax, 1:nS);
    yticklabels(ax, yLab);
    ax.TickLength = [0 0];
    ax.Layer = 'top';
    box(ax, 'on');
    ax.FontName = cfg.fontName; ax.FontSize = cfg.fontSize;
    ax.LineWidth = cfg.axWidth;
    ax.TickLabelInterpreter = 'latex';
    xlabel(ax, '$n$',      'Interpreter', 'latex', 'FontSize', cfg.fontSize + 1);
    ylabel(ax, '$\sigma$', 'Interpreter', 'latex', 'FontSize', cfg.fontSize + 1);

    lgd = legend(ax, [hR hG], 'Location', 'northoutside', ...
        'Orientation', 'horizontal', 'Interpreter', 'latex');
    lgd.FontSize = cfg.fontSize - 1; lgd.Box = 'off';

    base = fullfile(outDir, sprintf('rs_actually_useful_matrix_mindeg%d', md));
    for f = cfg.formats
        ext = f{1};
        file = [base '.' ext];
        if strcmp(ext, 'png')
            exportgraphics(fig, file, 'Resolution', 600);
        else
            exportgraphics(fig, file, 'ContentType', 'vector');
        end
    end
    close(fig);
end
end

function drawCells(ax, xc, rc, col, fill)
% one patch call for many ticks (vector graphics, no raster image).
% No edge lines: the gap between ticks comes from the geometry (width = fill
% of the slot), so it scales with the tick size instead of swamping it.
if isempty(xc), return; end
hw = fill / 2;
X = [xc - hw, xc + hw, xc + hw, xc - hw]';
Y = [rc - 0.5, rc - 0.5, rc + 0.5, rc + 0.5]';
patch(ax, 'XData', X, 'YData', Y, 'FaceColor', col, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function ok = readOkColumn(file)
% Read the "ok" column (0/1 or true/false) of an rs_usefulness_structs CSV.
% Returns a column vector (NaN where a value cannot be parsed).
ok = NaN(0, 1);
try
    T = readtable(file);
catch ME
    warning('Could not read %s (%s).', file, ME.message);
    return
end
if ~ismember('ok', T.Properties.VariableNames)
    warning('%s has no "ok" column (found: %s).', file, ...
        strjoin(T.Properties.VariableNames, ', '));
    return
end
c = T.ok;
if iscell(c) || isstring(c) || iscategorical(c)
    s = lower(string(c));
    c = str2double(s);
    c(s == "true")  = 1;
    c(s == "false") = 0;
end
ok = double(c(:));
end