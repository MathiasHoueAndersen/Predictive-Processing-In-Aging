%% ================================================================================ %%
%% ================== 1. Robustness check: Sensor-level Analysis ==================
%% ================================================================================ %%

%% A: Compute average ERF by group and condition
clear; clc;

analysis_dir = '/mainpath/SensorAnalyis';
matfile   = fullfile(analysis_dir, 'Block_2_MMNsubtracted.mat');
groupfile = '/mainpath/groups.mat';

group_names = {'Young','Older'};
cond_names  = {'Global','Local'};

left_set1  = {'MEG0131','MEG0121'};
left_set2  = {'MEG0231','MEG0241'};
right_set1 = {'MEG1331','MEG1341'};
right_set2 = {'MEG1611','MEG1621'};

flip_set2 = true;

alpha = 0.05;
cluster_alpha = 0.05;
n_perm = 5000;

combo_colors = [ ...
    0.1 0.2 0.5;   % Young-Local
    0.5 0.1 0.1;   % Older-Local
    0.4 0.6 0.85;  % Young-Global
    0.8 0.2 0.2];  % Older-Global

clr_group = [0 0.8 0.2];
clr_cond  = [0.6 0 0.8];
clr_int   = [1 0.6 0];

ribbon_alpha = 0.20;
shade_gray = [0.7 0.7 0.7];
shade_alpha = 0.15;
ymin = -100;
ymax = 70;

S = load(matfile);
G = load(groupfile);

data_mat   = S.data_mat;
chanlabels = S.chanlabels;

young = G.young(:);
older = G.older(:);

group_id = nan(size(data_mat,3),1);
group_id(young) = 1;
group_id(older) = 2;

if isrow(chanlabels)
    chanlabels = chanlabels';
end

[nChan, nTime, nSubj, nCond] = size(data_mat);

sample_interval_s = 0.004;
paradigm_onset_s  = 0.400;

time_full = ((0:nTime-1) * sample_interval_s) - paradigm_onset_s;

xlim_target = [-0.1 0.8];
t_idx = find(time_full >= xlim_target(1) & time_full <= xlim_target(2));

if isempty(t_idx)
    error('No timepoints found within [%g, %g] s.', xlim_target(1), xlim_target(2));
end

if t_idx(1) > 1
    t_idx = [t_idx(1)-1, t_idx];
end

data_mat = data_mat(:, t_idx, :, :);
time = time_full(t_idx);
nTime = numel(time);

fprintf('Loaded data_mat: %d channels x %d displayed time points x %d subjects x %d conditions\n', ...
    nChan, nTime, nSubj, nCond);

if numel(group_id) ~= nSubj
    error('group_id must have length %d, but has length %d.', nSubj, numel(group_id));
end

if any(isnan(group_id))
    warning('Some participants are not assigned to young or older.');
end

if nCond ~= 2
    error('Expected 2 conditions, but data_mat has %d conditions.', nCond);
end

left_idx_set1  = find_channels(chanlabels, left_set1);
left_idx_set2  = find_channels(chanlabels, left_set2);
right_idx_set1 = find_channels(chanlabels, right_set1);
right_idx_set2 = find_channels(chanlabels, right_set2);

fprintf('\nLEFT set 1 channels:\n');
disp(chanlabels(left_idx_set1));

fprintf('LEFT set 2 channels:\n');
disp(chanlabels(left_idx_set2));

fprintf('RIGHT set 1 channels:\n');
disp(chanlabels(right_idx_set1));

fprintf('RIGHT set 2 channels:\n');
disp(chanlabels(right_idx_set2));

set1_idx = [left_idx_set1; right_idx_set1];
set2_idx = [left_idx_set2; right_idx_set2];
all_selected_idx = [set1_idx; set2_idx];

left_data_set1  = data_mat(left_idx_set1, :, :, :);
left_data_set2  = data_mat(left_idx_set2, :, :, :);
right_data_set1 = data_mat(right_idx_set1, :, :, :);
right_data_set2 = data_mat(right_idx_set2, :, :, :);

if flip_set2
    left_data_set2  = -left_data_set2;
    right_data_set2 = -right_data_set2;
end

left_data  = cat(1, left_data_set1, left_data_set2);
right_data = cat(1, right_data_set1, right_data_set2);

left_erf  = squeeze(mean(left_data, 1, 'omitnan'));
right_erf = squeeze(mean(right_data, 1, 'omitnan'));

subject_erf = squeeze(mean(cat(4, left_erf, right_erf), 4, 'omitnan'));

if ndims(subject_erf) == 2
    subject_erf = reshape(subject_erf, nTime, nSubj, nCond);
end

nGroups = 2;

group_erf = nan(nTime, nGroups, nCond);
group_sem = nan(nTime, nGroups, nCond);

for g = 1:nGroups
    subj_idx = group_id == g;

    if ~any(subj_idx)
        warning('No participants found for group %d.', g);
        continue;
    end

    for c = 1:nCond
        tmp = subject_erf(:, subj_idx, c);

        group_erf(:, g, c) = mean(tmp, 2, 'omitnan');
        group_sem(:, g, c) = std(tmp, 0, 2, 'omitnan') ./ sqrt(sum(subj_idx));
    end
end

YG = squeeze(subject_erf(:, young, 1))';
YL = squeeze(subject_erf(:, young, 2))';
OG = squeeze(subject_erf(:, older, 1))';
OL = squeeze(subject_erf(:, older, 2))';

Y_mean = (YL + YG) ./ 2;
O_mean = (OL + OG) ./ 2;

D_y   = YG - YL;
D_o   = OG - OL;
D_all = [D_y; D_o];

obs_G = mean(Y_mean, 1, 'omitnan') - mean(O_mean, 1, 'omitnan');
obs_C = mean(D_all, 1, 'omitnan');
obs_I = mean(D_y, 1, 'omitnan') - mean(D_o, 1, 'omitnan');

[sig_G, thr_G, clu_G] = perm_between_time(Y_mean, O_mean, n_perm, cluster_alpha, alpha);
[sig_I, thr_I, clu_I] = perm_between_time(D_y, D_o, n_perm, cluster_alpha, alpha);
[sig_C, thr_C, clu_C] = perm_within_time(D_all, n_perm, cluster_alpha, alpha);

RES = struct();
RES.time = time;
RES.alpha = alpha;
RES.cluster_alpha = cluster_alpha;
RES.n_perm = n_perm;
RES.obs = struct('Group', obs_G, 'Condition', obs_C, 'Interaction', obs_I);
RES.thresholds = struct('Group', thr_G, 'Condition', thr_C, 'Interaction', thr_I);
RES.sig = struct('Group', sig_G, 'Condition', sig_C, 'Interaction', sig_I);
RES.windows = struct( ...
    'Group', mask_to_windows(sig_G, time), ...
    'Condition', mask_to_windows(sig_C, time), ...
    'Interaction', mask_to_windows(sig_I, time));
RES.N = struct('Young', numel(young), 'Older', numel(older));
RES.clusters = struct('Group', clu_G, 'Condition', clu_C, 'Interaction', clu_I);
RES.time_full = time_full;
RES.selected_time_idx = t_idx;
RES.sample_interval_s = sample_interval_s;
RES.paradigm_onset_s = paradigm_onset_s;

save(fullfile(analysis_dir, 'ERF_group_averages_selected_channels.mat'), ...
    'group_erf', 'group_sem', 'subject_erf', ...
    'left_erf', 'right_erf', ...
    'group_id', 'group_names', 'cond_names', ...
    'all_selected_idx', 'set1_idx', 'set2_idx', ...
    'left_idx_set1', 'left_idx_set2', ...
    'right_idx_set1', 'right_idx_set2', ...
    'chanlabels', 'flip_set2', 'RES');

fprintf('\nSaved: %s\n', fullfile(analysis_dir, 'ERF_group_averages_selected_channels.mat'));

out_xlsx = fullfile(analysis_dir, 'ERF_statistical_test_outputs.xlsx');
export_stats_to_xlsx(out_xlsx, RES);
fprintf('Saved: %s\n', out_xlsx);

h = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 300 500]);

set(h, 'DefaultAxesFontName', 'Helvetica', ...
       'DefaultTextFontName', 'Helvetica', ...
       'DefaultLegendFontName', 'Helvetica', ...
       'DefaultAxesFontSize', 8, ...
       'DefaultTextFontSize', 8, ...
       'DefaultLegendFontSize', 8, ...
       'DefaultAxesTitleFontSizeMultiplier', 1, ...
       'DefaultAxesLabelFontSizeMultiplier', 1);

tlo = tiledlayout(2,1, 'Padding','compact','TileSpacing','compact');

ax = nexttile(tlo, 1);
hold(ax,'on');

wins_G = mask_to_windows(sig_G, time);
wins_C = mask_to_windows(sig_C, time);
wins_I = mask_to_windows(sig_I, time);

draw_gray_windows_7(ax, wins_G, ymin, ymax, shade_gray, shade_alpha);
draw_gray_windows_7(ax, wins_C, ymin, ymax, shade_gray, shade_alpha);
draw_gray_windows_7(ax, wins_I, ymin, ymax, shade_gray, shade_alpha);

for c = 1:nCond
    if c == 1
        young_color = combo_colors(3,:);
        older_color = combo_colors(4,:);
    else
        young_color = combo_colors(1,:);
        older_color = combo_colors(2,:);
    end

    shaded_sem(ax, time, group_erf(:,1,c), group_sem(:,1,c), young_color, ribbon_alpha);
    shaded_sem(ax, time, group_erf(:,2,c), group_sem(:,2,c), older_color, ribbon_alpha);

    plot(ax, time, group_erf(:,1,c), 'LineWidth', 1.2, 'Color', young_color);
    plot(ax, time, group_erf(:,2,c), 'LineWidth', 1.2, 'Color', older_color);
end

draw_bottom_lines_7(ax, ...
    {wins_G, wins_C, wins_I}, ...
    {clr_group, clr_cond, clr_int}, ...
    xlim_target, ymin, ymax);

ylim(ax, [ymin ymax]);
xlim(ax, xlim_target);
xticks(-0.1:0.1:0.8)
grid(ax,'on'); box(ax,'on');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Activation (a.u.)');
title(ax, 'Global and Local');

xline(ax, 0, ':', 'Color', [0 0 0], 'LineWidth', 0.75, 'HandleVisibility','off');

axL = nexttile(tlo, 2);
cla(axL); hold(axL,'on'); axis(axL,'off');

h_leg = gobjects(1,4);
for j = 1:4
    h_leg(j) = plot(axL, NaN, NaN, '-', 'LineWidth', 3, 'Color', combo_colors(j,:));
end

hE = gobjects(1,3);
hE(1) = plot(axL, NaN, NaN, '-', 'LineWidth', 3.0, 'Color', clr_group);
hE(2) = plot(axL, NaN, NaN, '-', 'LineWidth', 3.0, 'Color', clr_cond);
hE(3) = plot(axL, NaN, NaN, '-', 'LineWidth', 3.0, 'Color', clr_int);

lgd = legend(axL, ...
    [h_leg, hE], ...
    {'Younger - Local','Older - Local','Younger - Global','Older - Global', ...
     'Group main effect','Condition main effect','Interaction effect'}, ...
    'NumColumns', 4, ...
    'Orientation','horizontal', ...
    'Location','northwest');

set(lgd, 'Interpreter','none', 'Units','normalized', 'Location','none');
set(lgd, 'Position', [0.20 0.05 0.60 0.35]);
set(lgd, 'ItemTokenSize', [15, 12]);
lgd.TextColor = [0 0 0];
lgd.Color     = [1 1 1];
lgd.EdgeColor = 'none';

sgtitle(tlo, 'Average ERF by group and condition');

out_png = fullfile(analysis_dir, 'ERF_group_averages_selected_channels.png');
print(h, out_png, '-dpng', '-r200');
close(h);

fprintf('Saved: %s\n', out_png);

%% ===================== HELPER FUNCTIONS =====================

function idx = find_channels(chanlabels, channel_suffixes)
    idx = [];

    for i = 1:numel(channel_suffixes)
        suffix = channel_suffixes{i};
        pattern = [suffix '$'];
        hit = find(~cellfun(@isempty, regexp(chanlabels, pattern, 'once')));

        if isempty(hit)
            error('Could not find channel ending with "%s". Check chanlabels.', suffix);
        elseif numel(hit) > 1
            fprintf('\nAmbiguous match for channel suffix %s:\n', suffix);
            disp(chanlabels(hit));
            error('More than one channel matched suffix "%s". Use full labels instead.');
        end

        idx = [idx; hit];
    end
end

function shaded_sem(ax, t, mu, se, color_rgb, alpha_val)
    mu = mu(:);
    se = se(:);
    t  = t(:);

    lo = mu - se;
    hi = mu + se;

    patch(ax, [t; flipud(t)], [lo; flipud(hi)], color_rgb, ...
        'FaceAlpha', alpha_val, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off', ...
        'Clipping', 'on');
end

function [sig_mask, thr, clusters] = perm_between_time(A, B, nperm, c_alpha, fwer_alpha)
    Na = size(A,1);
    Nb = size(B,1);
    T = size(A,2);

    if size(B,2) ~= T
        error('Time dim mismatch A/B.');
    end

    obs = mean(A,1,'omitnan') - mean(B,1,'omitnan');
    X = [A; B];

    pooled_abs = nan(nperm*T,1);
    for p = 1:nperm
        idx = randperm(Na+Nb);
        Ap  = X(idx(1:Na), :);
        Bp  = X(idx(Na+1:end), :);
        dpp = mean(Ap,1,'omitnan') - mean(Bp,1,'omitnan');
        pooled_abs((p-1)*T+1 : p*T) = abs(dpp(:));
    end

    cft = quantile_fast(pooled_abs, 1 - c_alpha);

    max_masses = nan(nperm,1);
    for p = 1:nperm
        idx = randperm(Na+Nb);
        Ap  = X(idx(1:Na), :);
        Bp  = X(idx(Na+1:end), :);
        dpp = mean(Ap,1,'omitnan') - mean(Bp,1,'omitnan');
        max_masses(p) = max_cluster_mass(abs(dpp), cft);
    end

    crit_mass = quantile_fast(max_masses, 1 - fwer_alpha);

    [clu_idx, clu_mass] = find_clusters(abs(obs), cft);
    sig_flags = clu_mass > crit_mass;

    sig_mask = false(1,T);
    for k = 1:numel(clu_idx)
        if sig_flags(k)
            sig_mask(clu_idx{k}) = true;
        end
    end

    thr = struct('cluster_forming', cft, 'cluster_mass', crit_mass);
    clusters = struct('idx_list', {clu_idx}, 'masses', clu_mass, 'sig_flags', sig_flags);
end

function [sig_mask, thr, clusters] = perm_within_time(D, nperm, c_alpha, fwer_alpha)
    [Ns, T] = size(D);
    obs = mean(D,1,'omitnan');

    flips = (randi(2, Ns, nperm)*2 - 3);

    pooled_abs = nan(nperm*T,1);
    for p = 1:nperm
        Dp  = D .* flips(:,p);
        dpp = mean(Dp,1,'omitnan');
        pooled_abs((p-1)*T+1 : p*T) = abs(dpp(:));
    end

    cft = quantile_fast(pooled_abs, 1 - c_alpha);

    max_masses = nan(nperm,1);
    for p = 1:nperm
        Dp  = D .* flips(:,p);
        dpp = mean(Dp,1,'omitnan');
        max_masses(p) = max_cluster_mass(abs(dpp), cft);
    end

    crit_mass = quantile_fast(max_masses, 1 - fwer_alpha);

    [clu_idx, clu_mass] = find_clusters(abs(obs), cft);
    sig_flags = clu_mass > crit_mass;

    sig_mask = false(1,T);
    for k = 1:numel(clu_idx)
        if sig_flags(k)
            sig_mask(clu_idx{k}) = true;
        end
    end

    thr = struct('cluster_forming', cft, 'cluster_mass', crit_mass);
    clusters = struct('idx_list', {clu_idx}, 'masses', clu_mass, 'sig_flags', sig_flags);
end

function m = max_cluster_mass(abs_series, cft)
    [~, masses] = find_clusters(abs_series, cft);
    if isempty(masses)
        m = 0;
    else
        m = max(masses);
    end
end

function [idx_list, masses] = find_clusters(abs_series, cft)
    mask = abs_series(:)' > cft;
    idx  = find(mask);

    idx_list = {};
    masses = [];

    if isempty(idx)
        return;
    end

    br = [1, find(diff(idx) > 1)+1, numel(idx)+1];

    for b = 1:numel(br)-1
        seg = idx(br(b):br(b+1)-1);
        idx_list{end+1} = seg;
        masses(end+1) = sum(abs_series(seg));
    end
end

function W = mask_to_windows(mask, tvec)
    mask = mask(:)';
    on = find(mask);
    W = [];

    if isempty(on)
        return;
    end

    br = [1, find(diff(on) > 1)+1, numel(on)+1];

    for b = 1:numel(br)-1
        seg = on(br(b):br(b+1)-1);
        W(end+1, :) = [tvec(seg(1)), tvec(seg(end))];
    end
end

function q = quantile_fast(x, p)
    x = x(isfinite(x));

    if isempty(x)
        q = NaN;
        return;
    end

    q = prctile(x, p*100);
end

function draw_gray_windows_7(ax, wins, ymin, ymax, shade_gray, shade_alpha)
    if isempty(wins)
        return;
    end

    for i = 1:size(wins,1)
        t1 = wins(i,1);
        t2 = wins(i,2);

        X = [t1 t2 t2 t1];
        Y = [ymin ymin ymax ymax];

        patch('XData', X, 'YData', Y, ...
            'FaceColor', shade_gray, ...
            'FaceAlpha', shade_alpha, ...
            'EdgeColor', 'none', ...
            'Parent', ax, ...
            'HitTest', 'off', ...
            'Clipping', 'on', ...
            'HandleVisibility', 'off');
    end
end

function draw_bottom_lines_7(ax, wins_cell, effect_colors, xlim_ok, ymin, ymax)
    oldUnits = get(ax,'Units');
    set(ax,'Units','pixels');
    axPos = get(ax,'Position');
    set(ax,'Units',oldUnits);

    yRange = ymax - ymin;
    spacing_factor = 5;
    dy_data = (axPos(4) > 0) * (spacing_factor / max(axPos(4),1)) * yRange;

    if dy_data == 0
        dy_data = 0.003*yRange;
    end

    y_base = ymin + 0.03 * yRange;
    n_levels = 3;

    for lvl = 1:3
        wins = wins_cell{lvl};

        if isempty(wins)
            continue;
        end

        y_line = y_base + (n_levels - lvl)*dy_data;

        for i = 1:size(wins,1)
            t1 = max(wins(i,1), xlim_ok(1));
            t2 = min(wins(i,2), xlim_ok(2));

            if t2 <= t1
                continue;
            end

            plot(ax, [t1 t2], [y_line y_line], '-', ...
                'LineWidth', 3.0, ...
                'Color', effect_colors{lvl}, ...
                'HitTest','off');
        end
    end
end

function export_stats_to_xlsx(out_xlsx, RES)
    if exist(out_xlsx, 'file')
        delete(out_xlsx);
    end

    effects = {'Group','Condition','Interaction'};

    settings = { ...
        'alpha', RES.alpha; ...
        'cluster_alpha', RES.cluster_alpha; ...
        'n_perm', RES.n_perm; ...
        'N_Young', RES.N.Young; ...
        'N_Older', RES.N.Older; ...
        'sample_interval_s', RES.sample_interval_s; ...
        'paradigm_onset_s', RES.paradigm_onset_s};

    writecell([{'Parameter','Value'}; settings], out_xlsx, 'Sheet', 'Settings');

    effect_col = {};
    cft_col = [];
    cm_col = [];

    for e = 1:numel(effects)
        eff = effects{e};
        effect_col{end+1,1} = eff;
        cft_col(end+1,1) = RES.thresholds.(eff).cluster_forming;
        cm_col(end+1,1) = RES.thresholds.(eff).cluster_mass;
    end

    T_thr = table(effect_col, cft_col, cm_col, ...
        'VariableNames', {'Effect','ClusterFormingThreshold','CriticalClusterMass'});
    writetable(T_thr, out_xlsx, 'Sheet', 'Thresholds');

    T_time = table( ...
        RES.time(:), ...
        RES.obs.Group(:), logical(RES.sig.Group(:)), ...
        RES.obs.Condition(:), logical(RES.sig.Condition(:)), ...
        RES.obs.Interaction(:), logical(RES.sig.Interaction(:)), ...
        'VariableNames', { ...
        'Time_s', ...
        'Observed_Group', 'Significant_Group', ...
        'Observed_Condition', 'Significant_Condition', ...
        'Observed_Interaction', 'Significant_Interaction'});

    writetable(T_time, out_xlsx, 'Sheet', 'TimepointResults');

    win_effect = {};
    win_number = [];
    win_start = [];
    win_end = [];

    for e = 1:numel(effects)
        eff = effects{e};
        W = RES.windows.(eff);

        if isempty(W)
            win_effect{end+1,1} = eff;
            win_number(end+1,1) = NaN;
            win_start(end+1,1) = NaN;
            win_end(end+1,1) = NaN;
        else
            for k = 1:size(W,1)
                win_effect{end+1,1} = eff;
                win_number(end+1,1) = k;
                win_start(end+1,1) = W(k,1);
                win_end(end+1,1) = W(k,2);
            end
        end
    end

    T_win = table(win_effect, win_number, win_start, win_end, ...
        'VariableNames', {'Effect','WindowIndex','Start_s','End_s'});
    writetable(T_win, out_xlsx, 'Sheet', 'SignificantWindows');

    clu_effect = {};
    clu_number = [];
    clu_start = [];
    clu_end = [];
    clu_n_timepoints = [];
    clu_mass = [];
    clu_significant = [];

    for e = 1:numel(effects)
        eff = effects{e};
        C = RES.clusters.(eff);

        if isempty(C.idx_list)
            clu_effect{end+1,1} = eff;
            clu_number(end+1,1) = NaN;
            clu_start(end+1,1) = NaN;
            clu_end(end+1,1) = NaN;
            clu_n_timepoints(end+1,1) = NaN;
            clu_mass(end+1,1) = NaN;
            clu_significant(end+1,1) = false;
        else
            for k = 1:numel(C.idx_list)
                idx = C.idx_list{k};

                clu_effect{end+1,1} = eff;
                clu_number(end+1,1) = k;
                clu_start(end+1,1) = RES.time(idx(1));
                clu_end(end+1,1) = RES.time(idx(end));
                clu_n_timepoints(end+1,1) = numel(idx);
                clu_mass(end+1,1) = C.masses(k);
                clu_significant(end+1,1) = logical(C.sig_flags(k));
            end
        end
    end

    T_clu = table(clu_effect, clu_number, clu_start, clu_end, ...
        clu_n_timepoints, clu_mass, logical(clu_significant), ...
        'VariableNames', {'Effect','ClusterIndex','Start_s','End_s', ...
        'N_Timepoints','ClusterMass','Significant'});
    writetable(T_clu, out_xlsx, 'Sheet', 'Clusters');
end

%% ================================================================================ %%
%% ================== 2. Robustness check: Behavioural Analysis ===================
%% ================================================================================ %%

% Purpose:
% Parse PsychoPy .log files and compute target detection performance.
%
% Target tone:     Sound 100_1 started
% Response:        Keypress: 1
% Hit window:      0–2 seconds after target onset

% Participant files in different format/naming convention:
% 0002bis.log, 0008bis.log (without practice), 0016 does not exist, 0019bis.log (without practice),
% 0021bis.log (without practice), 0026bis.log (without practice), 0028bis.log (without practice), 0029bis.log (without practice), 
% 0034bis.log (without practice)
% 0044bis.log (without practice), 0044bis.log (without practice), 0049bis.log (without practice)
% 0069bis.log (with practice), 0071bis.log (without pratice)
% 0077tris.log (without practice)

clear; clc; close all;

set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');
set(groot, 'defaultLegendFontName', 'Helvetica');

% ---------------- USER SETTINGS ----------------

rootDir = '/mainpath/BehaviouralDataAnalysis/TSA2021_BehavioralMEG';

resultsRootDir = '/mainpath/BehaviouralDataAnalysis/Results';

outDir = fullfile(resultsRootDir, 'TargetDetection_Tables');
groupPlotDir = fullfile(resultsRootDir, 'Group_Level_Plots');
ageGroupPlotDir = fullfile(resultsRootDir, 'AgeGroup_Plots');
participantPlotDir = fullfile(resultsRootDir, 'Participant_Plots');

dirsToMake = {outDir, groupPlotDir, ageGroupPlotDir, participantPlotDir};

for d = 1:numel(dirsToMake)
    if ~exist(dirsToMake{d}, 'dir')
        mkdir(dirsToMake{d});
    end
end

responseWindow = [0 3];   % seconds after target onset
targetPattern   = 'Sound 100_1 started';
responsePattern = 'Keypress: 1';
expStartPattern = sprintf('EXP \tSound  started');

% ---------------- LOAD AGE GROUPS ----------------

load('/mainpath/groups.mat');

older = double(older(:));
young = double(young(:));

% ---------------- FIND LOG FILES ----------------

logFilesAll = dir(fullfile(rootDir, '**', '*.log'));

specialLogNames = { ...
    '0002bis.log', ...
    '0008bis.log', ...
    '0011bis.log', ...
    '0019bis.log', ...
    '0021bis.log', ...
    '0026bis.log', ...
    '0028bis.log', ...
    '0029bis.log', ...
    '0034bis.log', ...
    '0042bis.log', ...
    '0044bis.log', ...
    '0049bis.log', ...
    '0069bis.log', ...
    '0071bis.log', ...
    '0077tris.log'};

noPracticeLogNames = { ...
    '0002bis.log', ...
    '0008bis.log', ...
    '0011bis.log', ...
    '0019bis.log', ...
    '0021bis.log', ...
    '0026bis.log', ...
    '0028bis.log', ...
    '0029bis.log', ...
    '0034bis.log', ...
    '0042bis.log', ...
    '0044bis.log', ...
    '0049bis.log', ...
    '0065.log', ...
    '0071bis.log', ...
    '0077tris.log'};

validLogName = false(numel(logFilesAll), 1);

excludedLogs = { ...
    '0011.log'};

for lf = 1:numel(logFilesAll)

    validLogName(lf) = ...
        (~isempty(regexp(logFilesAll(lf).name, ...
        '^(000[1-9]|00[1-9][0-9])\.log$', 'once')) || ...
        ismember(logFilesAll(lf).name, specialLogNames)) && ...
        ~ismember(logFilesAll(lf).name, excludedLogs);

end


logFiles = logFilesAll(validLogName);

fprintf('Found %d log files.\n', numel(logFiles));

allTargetRows = table();
summaryRows   = table();

% ---------------- PROCESS FILES ----------------

for f = 1:numel(logFiles)

    logPath = fullfile(logFiles(f).folder, logFiles(f).name);
    fprintf('\nProcessing: %s\n', logPath);

    txt = fileread(logPath);
    lines = regexp(txt, '\r\n|\n|\r', 'split')';

    targetTimes = [];
    responseTimes = [];
    expStartCount = 0;

    if strcmp(logFiles(f).name, '0032.log') || strcmp(logFiles(f).name, '0033.log') || ...
            ismember(logFiles(f).name, noPracticeLogNames)
        requiredExpStartCount = 1;
    else
        requiredExpStartCount = 2;
    end

    for i = 1:numel(lines)

        line = strtrim(lines{i});
        if isempty(line)
            continue
        end

        tok = regexp(line, '^([0-9]+\.?[0-9]*)', 'tokens', 'once');
        if isempty(tok)
            continue
        end

        t = str2double(tok{1});

        if contains(line, expStartPattern)
            expStartCount = expStartCount + 1;
            continue
        end

        if expStartCount < requiredExpStartCount
            continue
        end

        if contains(line, targetPattern)
            targetTimes(end+1,1) = t; %#ok<SAGROW>
        end

        if contains(line, responsePattern)
            responseTimes(end+1,1) = t; %#ok<SAGROW>
        end
    end

    subjTok = regexp(logFiles(f).name, '(\d{1,4})', 'tokens', 'once');

    if isempty(subjTok)
        subjTok = regexp(logFiles(f).folder, '(\d{1,4})', 'tokens', 'once');
    end

    if isempty(subjTok)
        subjID = sprintf('unknown_%03d', f);
        subjNum = NaN;
    else
        subjID = subjTok{1};
        subjNum = str2double(subjID);
    end

    if subjNum == 1
        fprintf('Skipping participant 1: %s\n', logFiles(f).name);
        continue
    end

    if ismember(subjNum, older)
        ageGroup = "older";
    elseif ismember(subjNum, young)
        ageGroup = "young";
    else
        ageGroup = "unknown";
        warning('Subject %s / %d not found in older or young group arrays.', subjID, subjNum);
    end

    if strcmp(logFiles(f).name, '0032.log') || strcmp(logFiles(f).name, '0033.log')
        blockID = NaN;
    else
        blockTok = regexp(logFiles(f).name, 'Block[_ ]?(\d+)', 'tokens', 'once');
        if isempty(blockTok)
            blockID = NaN;
        else
            blockID = str2double(blockTok{1});
        end
    end

    usedResponses = false(size(responseTimes));

    nTargets = numel(targetTimes);
    hit       = false(nTargets,1);
    rt        = nan(nTargets,1);
    matchedResponseTime = nan(nTargets,1);

    for tIdx = 1:nTargets

        t0 = targetTimes(tIdx);

        candidateIdx = find( ...
            responseTimes >= t0 + responseWindow(1) & ...
            responseTimes <= t0 + responseWindow(2) & ...
            ~usedResponses);

        if ~isempty(candidateIdx)
            chosenIdx = candidateIdx(1);

            hit(tIdx) = true;
            matchedResponseTime(tIdx) = responseTimes(chosenIdx);
            rt(tIdx) = responseTimes(chosenIdx) - t0;

            usedResponses(chosenIdx) = true;
        end
    end

    falseAlarmTimes = responseTimes(~usedResponses);
    nFalseAlarms = numel(falseAlarmTimes);

    fileCol  = repmat(string(logFiles(f).name), nTargets, 1);
    subjCol  = repmat(string(subjID), nTargets, 1);
    subjNumCol = repmat(subjNum, nTargets, 1);
    groupCol = repmat(ageGroup, nTargets, 1);
    blockCol = repmat(blockID, nTargets, 1);

    targetTrial = (1:nTargets)';

    thisTargets = table( ...
        subjCol, subjNumCol, groupCol, fileCol, blockCol, targetTrial, ...
        targetTimes, matchedResponseTime, rt, hit, ...
        'VariableNames', { ...
        'subj', 'subj_num', 'age_group', 'file', 'block', 'target_trial', ...
        'target_time', 'response_time', 'RT', 'hit'});

    allTargetRows = [allTargetRows; thisTargets]; %#ok<AGROW>

    nHits   = sum(hit);
    nMisses = sum(~hit);

    meanRT   = mean(rt(hit), 'omitnan');
    medianRT = median(rt(hit), 'omitnan');
    sdRT     = std(rt(hit), 'omitnan');

    thisSummary = table( ...
        string(subjID), subjNum, ageGroup, string(logFiles(f).name), blockID, ...
        nTargets, nHits, nMisses, nFalseAlarms, ...
        meanRT, medianRT, sdRT, ...
        'VariableNames', { ...
        'subj', 'subj_num', 'age_group', 'file', 'block', ...
        'n_targets', 'n_hits', 'n_misses', 'n_false_alarms', ...
        'mean_RT', 'median_RT', 'sd_RT'});

    summaryRows = [summaryRows; thisSummary]; %#ok<AGROW>

    if nTargets < 5 || nTargets > 20
        warning('Unexpected target count in %s: %d targets', logFiles(f).name, nTargets);
    end
end

% ---------------- SAVE TABLES ----------------

targetTablePath  = fullfile(outDir, 'target_level_results.csv');
summaryTablePath = fullfile(outDir, 'summary_results.csv');

writetable(allTargetRows, targetTablePath);
writetable(summaryRows, summaryTablePath);

fprintf('\nSaved:\n%s\n%s\n', targetTablePath, summaryTablePath);

% ---------------- AGGREGATE BY PARTICIPANT ----------------

subjs = unique(allTargetRows.subj);

participantRows = table();

for s = 1:numel(subjs)

    subj = subjs(s);
    idx = allTargetRows.subj == subj;

    nTargets = sum(idx);
    nHits = sum(allTargetRows.hit(idx));
    nMisses = nTargets - nHits;

    subjRTs = allTargetRows.RT(idx & allTargetRows.hit);

    subjSummaryIdx = summaryRows.subj == subj;

    nFalseAlarms = sum(summaryRows.n_false_alarms(subjSummaryIdx));

    tmp = table( ...
        subj, ...
        allTargetRows.subj_num(find(idx,1)), ...
        allTargetRows.age_group(find(idx,1)), ...
        nTargets, nHits, nMisses, nFalseAlarms, ...
        nHits / nTargets, ...
        mean(subjRTs, 'omitnan'), ...
        median(subjRTs, 'omitnan'), ...
        std(subjRTs, 'omitnan'), ...
        'VariableNames', { ...
        'subj', 'subj_num', 'age_group', ...
        'n_targets', 'n_hits', 'n_misses', 'n_false_alarms', ...
        'hit_rate', 'mean_RT', 'median_RT', 'sd_RT'});

    participantRows = [participantRows; tmp]; %#ok<AGROW>
end

participantSummaryPath = fullfile(outDir, 'participant_summary.csv');
writetable(participantRows, participantSummaryPath);

fprintf('Participants recognized in total: %d\n', height(participantRows));

% ---------------- GROUP-LEVEL PLOTS ----------------

validRT = allTargetRows.RT(allTargetRows.hit);

fig = figure('Color','w', 'Visible','off');
histogram(validRT, 20);
xlabel('Reaction time (s)');
ylabel('Count');
title('Reaction-time distribution: all participants');
moveTitle1x();
grid on;
saveas(fig, fullfile(groupPlotDir, 'GROUP_RT_distribution.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
boxchart(categorical(allTargetRows.subj(allTargetRows.hit)), allTargetRows.RT(allTargetRows.hit), 'MarkerStyle', 'none');
xlabel('Participant');
ylabel('Reaction time (s)');
title('Reaction times by participant');
moveTitle1x();
grid on;
saveas(fig, fullfile(groupPlotDir, 'GROUP_RT_by_participant.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
bar(categorical(participantRows.subj), participantRows.hit_rate);
ylim([0 1]);
xlabel('Participant');
ylabel('Hit rate');
title('Target detection accuracy by participant');
moveTitle1x();
grid on;
saveas(fig, fullfile(groupPlotDir, 'GROUP_hit_rate_by_participant.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
counts = [sum(allTargetRows.hit), sum(~allTargetRows.hit)];
bar(categorical({'Hits','Misses'}), counts);
ylabel('Count');
title('Correctness distribution: all targets');
moveTitle1x();
grid on;
saveas(fig, fullfile(groupPlotDir, 'GROUP_correctness_counts.png'));
close(fig);

% ---------------- OLDER VS YOUNGER PLOTS ----------------

knownGroupIdx = participantRows.age_group == "older" | participantRows.age_group == "young";
participantRowsKnown = participantRows(knownGroupIdx,:);

rng(1);
jitterAmount = 0.18;
ageGroups = ["young"; "older"];
ageGroupLabels = ["Younger"; "Older"];
youngerColor = [0.2 0.2 0.8];
olderColor   = [0.8 0.2 0.2];
ageGroupColors = [youngerColor; olderColor];
ageGroupCat = categorical(participantRowsKnown.age_group, ageGroups, ageGroupLabels);


% ---------------- HIT / FALSE-ALARM OUTLIER CHECK ----------------
% Outliers are defined across all known-age participants combined
% (younger + older), not separately within each age group.

hitMean = mean(participantRowsKnown.n_hits, 'omitnan');
hitSD   = std(participantRowsKnown.n_hits, 'omitnan');
hitLower = hitMean - 3 * hitSD;
hitUpper = hitMean + 3 * hitSD;

falseAlarmMean = mean(participantRowsKnown.n_false_alarms, 'omitnan');
falseAlarmSD   = std(participantRowsKnown.n_false_alarms, 'omitnan');
falseAlarmLower = falseAlarmMean - 3 * falseAlarmSD;
falseAlarmUpper = falseAlarmMean + 3 * falseAlarmSD;

hitOutlierIdx = participantRowsKnown.n_hits < hitLower | ...
                participantRowsKnown.n_hits > hitUpper;

falseAlarmOutlierIdx = participantRowsKnown.n_false_alarms < falseAlarmLower | ...
                       participantRowsKnown.n_false_alarms > falseAlarmUpper;

outlierSubjsHits = participantRowsKnown.subj(hitOutlierIdx);
outlierSubjsFalseAlarms = participantRowsKnown.subj(falseAlarmOutlierIdx);

excludeOutlierIdx = hitOutlierIdx | falseAlarmOutlierIdx;
outlierSubjsAny = participantRowsKnown.subj(excludeOutlierIdx);

fprintf('OUTLIER CHECK: HITS');
fprintf('Mean hits = %.4f, SD = %.4f, lower cutoff = %.4f, upper cutoff = %.4f', ...
    hitMean, hitSD, hitLower, hitUpper);

if isempty(outlierSubjsHits)
    fprintf('No hit outliers found using +/- 3 SD.');
else
    fprintf('Hit outlier subject IDs (+/- 3 SD):');
    disp(outlierSubjsHits);
end

fprintf('OUTLIER CHECK: FALSE ALARMS');
fprintf('Mean false alarms = %.4f, SD = %.4f, lower cutoff = %.4f, upper cutoff = %.4f', ...
    falseAlarmMean, falseAlarmSD, falseAlarmLower, falseAlarmUpper);

if isempty(outlierSubjsFalseAlarms)
    fprintf('No false-alarm outliers found using +/- 3 SD.');
else
    fprintf('False-alarm outlier subject IDs (+/- 3 SD):');
    disp(outlierSubjsFalseAlarms);
end

fprintf('Subjects excluded from additional outlier-excluded age-group plots:');
if isempty(outlierSubjsAny)
    fprintf('No subjects excluded.');
else
    disp(outlierSubjsAny);
end

participantRowsKnown_noOutliers = participantRowsKnown(~excludeOutlierIdx,:);
ageGroupCat_noOutliers = categorical(participantRowsKnown_noOutliers.age_group, ageGroups, ageGroupLabels);

% ---------------- AGE-GROUP MEAN ± SD BAR PLOTS ----------------

barPlotDir = fullfile(resultsRootDir, 'AgeGroup_BarPlots');

if ~exist(barPlotDir, 'dir')
    mkdir(barPlotDir);
end

barMetrics = { ...
    'mean_RT', ...
    'n_false_alarms', ...
    'n_misses', ...
    'n_hits'};

barMetricLabels = { ...
    'Reaction time (s)', ...
    'False alarms', ...
    'Misses', ...
    'Hits'};

barMetricTitles = { ...
    'Mean RT by age group', ...
    'Mean false alarms by age group', ...
    'Mean misses by age group', ...
    'Mean hits by age group'};

barMetricFiles = { ...
    'AGEGROUP_BAR_RT_mean_SD.png', ...
    'AGEGROUP_BAR_falsealarms_mean_SD.png', ...
    'AGEGROUP_BAR_misses_mean_SD.png', ...
    'AGEGROUP_BAR_hits_mean_SD.png'};

for m = 1:numel(barMetrics)

    metricName  = barMetrics{m};
    yLabelText  = barMetricLabels{m};
    titleText   = barMetricTitles{m};
    outFileName = barMetricFiles{m};

    groupMeans = nan(numel(ageGroups),1);

    for g = 1:numel(ageGroups)

        grp = ageGroups(g);
        idx = participantRowsKnown.age_group == grp & ...
              ~isnan(participantRowsKnown.(metricName));

        vals = participantRowsKnown.(metricName)(idx);

        groupMeans(g) = mean(vals, 'omitnan');
    end

    fig = figure('Color','w', 'Visible','off');

    b = bar(1:numel(ageGroups), groupMeans, 'FaceColor', 'flat');
    b.CData = ageGroupColors;

    hold on;


    xlim([0.5 numel(ageGroups)+0.5]);
    xticks(1:numel(ageGroups));
    xticklabels(cellstr(ageGroupLabels));

    xlabel('Age group');
    ylabel(yLabelText);
    title(titleText);
moveTitle1x();

    grid on;
    box off;

    set(gca, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12, ...
        'LineWidth', 1);

    saveas(fig, fullfile(barPlotDir, outFileName));
    close(fig);
end

fig = figure('Color','w', 'Visible','off');
boxchart(ageGroupCat, participantRowsKnown.mean_RT, 'MarkerStyle', 'none');
hold on;

rtPlotX = [];
rtPlotY = [];
rtPlotSubj = strings(0,1);

for g = 1:numel(ageGroups)
    grp = ageGroups(g);
    idx = participantRowsKnown.age_group == grp & ~isnan(participantRowsKnown.mean_RT);
    x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
    y = participantRowsKnown.mean_RT(idx);
    scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');

    rtPlotX = [rtPlotX; x]; %#ok<AGROW>
    rtPlotY = [rtPlotY; y]; %#ok<AGROW>
    rtPlotSubj = [rtPlotSubj; participantRowsKnown.subj(idx)]; %#ok<AGROW>
end

xlabel('Age group');
ylabel('Mean reaction time (s)');
title('Mean RT across participants: Older vs Younger');
moveTitle1x();
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'AGEGROUP_RT_boxplot.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
boxchart(ageGroupCat, participantRowsKnown.mean_RT, 'MarkerStyle', 'none');
hold on;

rng(1);
for g = 1:numel(ageGroups)
    grp = ageGroups(g);
    idx = participantRowsKnown.age_group == grp & ~isnan(participantRowsKnown.mean_RT);
    x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
    y = participantRowsKnown.mean_RT(idx);
    scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
    text(x + 0.03, y, cellstr(participantRowsKnown.subj(idx)), ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end

xlabel('Age group');
ylabel('Mean reaction time (s)');
title('Mean RT across participants: Older vs Younger');
moveTitle1x();
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'AGEGROUP_RT_boxplot_with_IDs.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
histogram(participantRowsKnown.mean_RT(participantRowsKnown.age_group == "young"), 15, 'FaceColor', youngerColor);
hold on;
histogram(participantRowsKnown.mean_RT(participantRowsKnown.age_group == "older"), 15, 'FaceColor', olderColor);
xlabel('Mean reaction time (s)');
ylabel('Participant count');
title('Mean RT across participants: Older vs Younger');
moveTitle1x();
legend({'Younger','Older'}, 'Location','best');
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'AGEGROUP_RT_histogram.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
histogram(participantRowsKnown.mean_RT(participantRowsKnown.age_group == "young"), 15, 'FaceColor', youngerColor);
hold on;
histogram(participantRowsKnown.mean_RT(participantRowsKnown.age_group == "older"), 15, 'FaceColor', olderColor);

youngIdx = participantRowsKnown.age_group == "young" & ~isnan(participantRowsKnown.mean_RT);
olderIdx = participantRowsKnown.age_group == "older" & ~isnan(participantRowsKnown.mean_RT);

scatter(participantRowsKnown.mean_RT(youngIdx), zeros(sum(youngIdx),1), ...
    45, youngerColor, 'filled', 'o', 'MarkerEdgeColor', 'k');
text(participantRowsKnown.mean_RT(youngIdx), zeros(sum(youngIdx),1), ...
    cellstr(participantRowsKnown.subj(youngIdx)), ...
    'FontSize', 8, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

scatter(participantRowsKnown.mean_RT(olderIdx), zeros(sum(olderIdx),1), ...
    45, olderColor, 'filled', 'o', 'MarkerEdgeColor', 'k');
text(participantRowsKnown.mean_RT(olderIdx), zeros(sum(olderIdx),1), ...
    cellstr(participantRowsKnown.subj(olderIdx)), ...
    'FontSize', 8, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

xlabel('Mean reaction time (s)');
ylabel('Participant count');
title('Mean RT across participants: Older vs Younger');
moveTitle1x();
legend({'Younger','Older'}, 'Location','best');
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'AGEGROUP_RT_histogram_with_IDs.png'));
close(fig);

metrics = { ...
    'n_hits', ...
    'n_misses', ...
    'n_false_alarms'};

metricLabels = { ...
    'Hits', ...
    'Misses', ...
    'False alarms'};

metricFiles = { ...
    'AGEGROUP_hits_boxplot.png', ...
    'AGEGROUP_misses_boxplot.png', ...
    'AGEGROUP_falsealarms_boxplot.png'};

metricFilesWithIDs = { ...
    'AGEGROUP_hits_boxplot_with_IDs.png', ...
    'AGEGROUP_misses_boxplot_with_IDs.png', ...
    'AGEGROUP_falsealarms_boxplot_with_IDs.png'};

for m = 1:numel(metrics)

    metricName = metrics{m};
    metricLabel = metricLabels{m};

    fig = figure('Color','w', 'Visible','off');

    yAll = participantRowsKnown.(metricName);
    boxchart(ageGroupCat, yAll, 'MarkerStyle', 'none');
    hold on;

    for g = 1:numel(ageGroups)
        grp = ageGroups(g);
        idx = participantRowsKnown.age_group == grp & ~isnan(participantRowsKnown.(metricName));

        x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
        y = participantRowsKnown.(metricName)(idx);

        scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
    end

    xlabel('Age group');
    ylabel('Count per participant');
    title(sprintf('%s by age group', metricLabel));
moveTitle1x();
    grid on;

    saveas(fig, fullfile(ageGroupPlotDir, metricFiles{m}));
    close(fig);

    fig = figure('Color','w', 'Visible','off');

    yAll = participantRowsKnown.(metricName);
    boxchart(ageGroupCat, yAll, 'MarkerStyle', 'none');
    hold on;

    rng(1);
    for g = 1:numel(ageGroups)
        grp = ageGroups(g);
        idx = participantRowsKnown.age_group == grp & ~isnan(participantRowsKnown.(metricName));

        x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
        y = participantRowsKnown.(metricName)(idx);

        scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
        text(x + 0.03, y, cellstr(participantRowsKnown.subj(idx)), ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    xlabel('Age group');
    ylabel('Count per participant');
    title(sprintf('%s by age group', metricLabel));
moveTitle1x();
    grid on;

    saveas(fig, fullfile(ageGroupPlotDir, metricFilesWithIDs{m}));
    close(fig);
end

% ---------------- OUTLIER-EXCLUDED AGE-GROUP PLOTS ----------------
% These plots exclude any participant identified as an outlier for either
% hits or false alarms using the +/- 3 SD criterion above.

fig = figure('Color','w', 'Visible','off');
boxchart(ageGroupCat_noOutliers, participantRowsKnown_noOutliers.mean_RT, 'MarkerStyle', 'none');
hold on;

rng(1);
for g = 1:numel(ageGroups)
    grp = ageGroups(g);
    idx = participantRowsKnown_noOutliers.age_group == grp & ...
          ~isnan(participantRowsKnown_noOutliers.mean_RT);
    x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
    y = participantRowsKnown_noOutliers.mean_RT(idx);
    scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
end

xlabel('Age group');
ylabel('Mean reaction time (s)');
title('Mean RT across participants: Older vs Younger, outliers excluded');
moveTitle1x();
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'AGEGROUP_RT_boxplot_outliers_excluded.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
boxchart(ageGroupCat_noOutliers, participantRowsKnown_noOutliers.mean_RT, 'MarkerStyle', 'none');
hold on;

rng(1);
for g = 1:numel(ageGroups)
    grp = ageGroups(g);
    idx = participantRowsKnown_noOutliers.age_group == grp & ...
          ~isnan(participantRowsKnown_noOutliers.mean_RT);
    x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
    y = participantRowsKnown_noOutliers.mean_RT(idx);
    scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
    text(x + 0.03, y, cellstr(participantRowsKnown_noOutliers.subj(idx)), ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end

xlabel('Age group');
ylabel('Mean reaction time (s)');
title('Mean RT across participants: Older vs Younger, outliers excluded');
moveTitle1x();
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'AGEGROUP_RT_boxplot_with_IDs_outliers_excluded.png'));
close(fig);

metricFilesOutliersExcluded = { ...
    'AGEGROUP_hits_boxplot_outliers_excluded.png', ...
    'AGEGROUP_misses_boxplot_outliers_excluded.png', ...
    'AGEGROUP_falsealarms_boxplot_outliers_excluded.png'};

metricFilesWithIDsOutliersExcluded = { ...
    'AGEGROUP_hits_boxplot_with_IDs_outliers_excluded.png', ...
    'AGEGROUP_misses_boxplot_with_IDs_outliers_excluded.png', ...
    'AGEGROUP_falsealarms_boxplot_with_IDs_outliers_excluded.png'};

for m = 1:numel(metrics)

    metricName = metrics{m};
    metricLabel = metricLabels{m};

    fig = figure('Color','w', 'Visible','off');

    yAll = participantRowsKnown_noOutliers.(metricName);
    boxchart(ageGroupCat_noOutliers, yAll, 'MarkerStyle', 'none');
    hold on;

    rng(1);
    for g = 1:numel(ageGroups)
        grp = ageGroups(g);
        idx = participantRowsKnown_noOutliers.age_group == grp & ...
              ~isnan(participantRowsKnown_noOutliers.(metricName));

        x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
        y = participantRowsKnown_noOutliers.(metricName)(idx);

        scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
    end

    xlabel('Age group');
    ylabel('Count per participant');
    title(sprintf('%s by age group, outliers excluded', metricLabel));
moveTitle1x();
    grid on;

    saveas(fig, fullfile(ageGroupPlotDir, metricFilesOutliersExcluded{m}));
    close(fig);

    fig = figure('Color','w', 'Visible','off');

    yAll = participantRowsKnown_noOutliers.(metricName);
    boxchart(ageGroupCat_noOutliers, yAll, 'MarkerStyle', 'none');
    hold on;

    rng(1);
    for g = 1:numel(ageGroups)
        grp = ageGroups(g);
        idx = participantRowsKnown_noOutliers.age_group == grp & ...
              ~isnan(participantRowsKnown_noOutliers.(metricName));

        x = g + (rand(sum(idx),1) - 0.5) * 2 * jitterAmount;
        y = participantRowsKnown_noOutliers.(metricName)(idx);

        scatter(x, y, 45, ageGroupColors(g,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
        text(x + 0.03, y, cellstr(participantRowsKnown_noOutliers.subj(idx)), ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    xlabel('Age group');
    ylabel('Count per participant');
    title(sprintf('%s by age group, outliers excluded', metricLabel));
moveTitle1x();
    grid on;

    saveas(fig, fullfile(ageGroupPlotDir, metricFilesWithIDsOutliersExcluded{m}));
    close(fig);
end

fig = figure('Color','w', 'Visible','off');
bar(categorical(participantRowsKnown.subj), ...
    [participantRowsKnown.n_hits, participantRowsKnown.n_misses, participantRowsKnown.n_false_alarms]);
xlabel('Participant');
ylabel('Count');
title('Hits, misses, and false alarms per participant');
moveTitle1x();
legend({'Hits','Misses','False alarms'}, 'Location','best');
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'PARTICIPANT_hits_misses_falsealarms.png'));
close(fig);

fig = figure('Color','w', 'Visible','off');
bar(categorical(participantRowsKnown.subj), ...
    [participantRowsKnown.n_hits, participantRowsKnown.n_misses, participantRowsKnown.n_false_alarms]);
xlabel('Participant');
ylabel('Count');
title('Hits, misses, and false alarms per participant');
moveTitle1x();
legend({'Hits','Misses','False alarms'}, 'Location','best');
grid on;
saveas(fig, fullfile(ageGroupPlotDir, 'PARTICIPANT_hits_misses_falsealarms_with_IDs.png'));
close(fig);

% ---------------- PARTICIPANT-LEVEL PLOTS ----------------

for s = 1:numel(subjs)

    subj = subjs(s);
    idx = allTargetRows.subj == subj;

    T = allTargetRows(idx,:);
    T = sortrows(T, {'file','target_trial'});

    x = (1:height(T))';

    % RT across target trials

    fig = figure('Color','w', 'Visible','off');

    hitIdx = T.hit;
    missIdx = ~T.hit;

    plot(x(hitIdx), T.RT(hitIdx), 'o-', 'LineWidth', 1.5);
    hold on;

    if any(missIdx)
        yMiss = repmat(responseWindow(2), sum(missIdx), 1);
        plot(x(missIdx), yMiss, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end

    xlabel('Target trial number');
    ylabel('Reaction time (s)');
    title(sprintf('Participant %s: RT across target trials', subj));
moveTitle1x();
    ylim([0 responseWindow(2) + 0.2]);
    grid on;

    legend({'Hit RT','Miss'}, 'Location','best');

    saveas(fig, fullfile(participantPlotDir, sprintf('Subj_%s_RT_trials.png', subj)));
    close(fig);

    % Correctness across target trials

    fig = figure('Color','w', 'Visible','off');

    stem(x, double(T.hit), 'filled');
    ylim([-0.1 1.1]);
    yticks([0 1]);
    yticklabels({'Miss','Hit'});
    xlabel('Target trial number');
    ylabel('Correctness');
    title(sprintf('Participant %s: correctness across target trials', subj));
moveTitle1x();
    grid on;

    saveas(fig, fullfile(participantPlotDir, sprintf('Subj_%s_correctness_trials.png', subj)));
    close(fig);

    % RT histogram per participant

    fig = figure('Color','w', 'Visible','off');

    histogram(T.RT(T.hit), 10);
    xlabel('Reaction time (s)');
    ylabel('Count');
    title(sprintf('Participant %s: RT distribution', subj));
moveTitle1x();
    grid on;

    saveas(fig, fullfile(participantPlotDir, sprintf('Subj_%s_RT_distribution.png', subj)));
    close(fig);
end

fprintf('\nAll results saved in:\n%s\n', resultsRootDir);

% ---------------- SANITY CHECK ----------------

fprintf('\nSANITY CHECK\n');
fprintf('Total targets found: %d\n', height(allTargetRows));
fprintf('Expected roughly: 3080\n');


disp(participantRows);


fprintf('\nFirst 20 response lines:\n');
respLines = lines(contains(lines, responsePattern));
disp(respLines(1:min(20,numel(respLines))));

fprintf('Targets: %d\n', numel(targetTimes));
fprintf('Responses: %d\n', numel(responseTimes));
fprintf('Hits: %d\n', sum(hit));
fprintf('False alarms: %d\n', nFalseAlarms);


figure;
histogram(diff(responseTimes), 100);
xlabel('Time between Keypress: 1 events (s)');
ylabel('Count');
title(sprintf('Participant %s: response intervals', subjID));
moveTitle1x();
grid on;


function moveTitle1x()

    t = get(gca,'Title');
    t.Units = 'normalized';
    pos = t.Position;
    pos(2) = 1.03;
    t.Position = pos;

end

%% ---------------- SIMPLE AGE-GROUP TESTS ----------------

testRows = participantRowsKnown( ...
    participantRowsKnown.age_group == "young" | ...
    participantRowsKnown.age_group == "older", :);

metricsToTest = {'mean_RT', 'n_hits', 'n_false_alarms'};
metricNames   = {'RT', 'Hits', 'False alarms'};

testResults = table();

fprintf('\nSIMPLE YOUNGER VS OLER ADULT TESTS\n');
fprintf('Independent-samples t-tests: younger vs older\n\n');

for m = 1:numel(metricsToTest)

    metric = metricsToTest{m};

    youngVals = testRows.(metric)(testRows.age_group == "young");
    olderVals = testRows.(metric)(testRows.age_group == "older");

    youngVals = youngVals(~isnan(youngVals));
    olderVals = olderVals(~isnan(olderVals));

    [~, p, ci, stats] = ttest2(youngVals, olderVals);

    youngMean = mean(youngVals, 'omitnan');
    olderMean = mean(olderVals, 'omitnan');

    youngSD = std(youngVals, 'omitnan');
    olderSD = std(olderVals, 'omitnan');

    fprintf('%s\n', metricNames{m});
    fprintf('Younger: n = %d, mean = %.4f, SD = %.4f\n', numel(youngVals), youngMean, youngSD);
    fprintf('Older:   n = %d, mean = %.4f, SD = %.4f\n', numel(olderVals), olderMean, olderSD);
    fprintf('t(%d) = %.4f, p = %.6f, 95%% CI = [%.4f, %.4f]\n\n', ...
        stats.df, stats.tstat, p, ci(1), ci(2));

    tmp = table( ...
        string(metricNames{m}), ...
        numel(youngVals), youngMean, youngSD, ...
        numel(olderVals), olderMean, olderSD, ...
        stats.tstat, stats.df, p, ci(1), ci(2), ...
        'VariableNames', { ...
        'measure', ...
        'n_young', 'mean_young', 'sd_young', ...
        'n_older', 'mean_older', 'sd_older', ...
        't_value', 'df', 'p_value', 'ci_lower', 'ci_upper'});

    testResults = [testResults; tmp]; %#ok<AGROW>
end

% FDR correction across the tested metrics
pvals = testResults.p_value;
[pSorted, sortIdx] = sort(pvals);
nTests = numel(pvals);

qSorted = pSorted .* nTests ./ (1:nTests)';
qSorted = min(qSorted, 1);

for i = nTests-1:-1:1
    qSorted(i) = min(qSorted(i), qSorted(i+1));
end

p_fdr = nan(size(pvals));
p_fdr(sortIdx) = qSorted;

testResults.p_fdr = p_fdr;
testResults.significant_fdr_05 = testResults.p_fdr < 0.05;

testResultsPath = fullfile(outDir, 'simple_age_group_ttests.csv');
writetable(testResults, testResultsPath);

testResultsXlsxPath = fullfile(outDir, 'simple_age_group_ttests.xlsx');
writetable(testResults, testResultsXlsxPath);

disp(testResults);

fprintf('Saved simple test results to:\n%s\n', testResultsPath);
fprintf('Saved simple test results to:\n%s\n', testResultsXlsxPath);
