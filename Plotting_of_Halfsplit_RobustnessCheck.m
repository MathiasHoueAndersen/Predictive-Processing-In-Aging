%% ================================================================================ %%
%% =============================== Halfsplit Analysis =============================
%% ================================================================================ %%

%% A: Plot BROADNESS time series: Older vs Younger, PC1-PC3
clear; clc;

data_dir = '/mainpath/HalfSplit';

broadness_files = { ...
    fullfile(data_dir, 'BROADNESS_Output_Half1.mat'), ...
    fullfile(data_dir, 'BROADNESS_Output_Half2.mat')};

halfsplit_names = {'HalfSplit1', 'HalfSplit2'};
split_nums = [1 2];

pc_list = 1:3;

groups_file = fullfile(data_dir, 'groups.mat');

out_dir = fullfile(data_dir, 'BROADNESS_Output', 'TimeSeries_GroupPlots_PC1_PC2_PC3');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_xlsx = fullfile(out_dir, 'BROADNESS_Timewise_ANOVA_ClusterPerm_AllSplits_PC1_PC2_PC3.xlsx');
if exist(out_xlsx, 'file')
    delete(out_xlsx);
end

G = load(groups_file);

if ~isfield(G, 'older') || ~isfield(G, 'young')
    error('groups.mat must contain variables named "older" and "young".');
end

older_labels = G.older(:)';
young_labels = G.young(:)';

tags   = {'Young_Local','Old_Local','Young_Global','Old_Global'};
labels = {'Younger - Local','Older - Local','Younger - Global','Older - Global'};

combo_colors = [ ...
    0.1 0.2 0.5;
    0.5 0.1 0.1;
    0.4 0.6 0.85;
    0.8 0.2 0.2];

ribbon_alpha = 0.20;

shade_gray  = [0.7 0.7 0.7];
shade_alpha = 0.15;

clr_group = [0 0.8 0.2];
clr_cond  = [0.6 0 0.8];
clr_int   = [1 0.6 0];

alpha          = 0.05;
cluster_alpha  = 0.05;
n_perm         = 5000;

sample_interval_s = 0.004;
zero_time_s       = 0.500;

xlim_target = [-0.1 0.8];
stats_xlim  = [0 0.8];
xtick_vals  = -0.1:0.1:0.8;

for hs = 1:numel(broadness_files)

    B = load(broadness_files{hs});

    if isfield(B, 'BROADNESS')
        BROADNESS = B.BROADNESS;
    else
        fn = fieldnames(B);
        BROADNESS = B.(fn{1});
    end

    TS = BROADNESS.TimeSeries_BrainNetworks;

    if ndims(TS) ~= 4
        error('Expected TimeSeries_BrainNetworks to be [time x PC x condition x participant], but got: %s', mat2str(size(TS)));
    end

    n_pc_available = size(TS, 2);

    for pcNumber = pc_list

        if pcNumber > n_pc_available
            error('Requested PC%d, but only %d PCs exist in %s.', ...
                pcNumber, n_pc_available, broadness_files{hs});
        end

        fprintf('\nProcessing %s, PC%d\n', halfsplit_names{hs}, pcNumber);

        % Select principal component
        X = squeeze(TS(:, pcNumber, :, :));

        if pcNumber == 2 || (pcNumber == 3 && hs == 1)
            X = -X;
        end

        if ndims(X) ~= 3
            error('Expected selected PC data to be [time x condition x participant], but got: %s', mat2str(size(X)));
        end

        n_time = size(X,1);
        n_cond = size(X,2);

        if n_cond ~= 2
            error('Expected 2 conditions: global and local. Found %d.', n_cond);
        end

        % Insert empty participant after participant 15
        empty_participant = nan(size(X,1), size(X,2), 1);
        X = cat(3, X(:,:,1:15), empty_participant, X(:,:,16:end));

        n_subj = size(X,3);

        older_idx = older_labels;
        young_idx = young_labels;

        if any(older_idx < 1 | older_idx > n_subj) || any(young_idx < 1 | young_idx > n_subj)
            error('One or more group labels exceed the participant dimension after inserting empty participant 16.');
        end

        time = ((0:n_time-1)' * sample_interval_s) - zero_time_s;

        t_plot_idx = find(time >= xlim_target(1) & time <= xlim_target(2));
        t_plot = time(t_plot_idx);

        if isempty(t_plot_idx)
            error('No samples found in requested time window %.3f to %.3f s.', xlim_target(1), xlim_target(2));
        end

        X_plot = X(t_plot_idx,:,:);
        dt_plot = median(diff(t_plot));

        t_stat_idx = find(time >= stats_xlim(1) & time <= stats_xlim(2));
        t_stat = time(t_stat_idx);

        if isempty(t_stat_idx)
            error('No samples found in stats window %.3f to %.3f s.', stats_xlim(1), stats_xlim(2));
        end

        X_stat = X(t_stat_idx,:,:);

        % condition 1 = global
        % condition 2 = local
        YG = squeeze(X_stat(:,1,young_idx))';
        YL = squeeze(X_stat(:,2,young_idx))';
        OG = squeeze(X_stat(:,1,older_idx))';
        OL = squeeze(X_stat(:,2,older_idx))';

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

        win_G = mask_to_windows(sig_G, t_stat);
        win_C = mask_to_windows(sig_C, t_stat);
        win_I = mask_to_windows(sig_I, t_stat);

        RES = struct();
        RES.HalfSplit     = halfsplit_names{hs};
        RES.Split         = split_nums(hs);
        RES.PC            = pcNumber;
        RES.time          = t_stat;
        RES.alpha         = alpha;
        RES.cluster_alpha = cluster_alpha;
        RES.n_perm        = n_perm;
        RES.obs           = struct('Group', obs_G, 'Condition', obs_C, 'Interaction', obs_I);
        RES.thresholds    = struct('Group', thr_G, 'Condition', thr_C, 'Interaction', thr_I);
        RES.sig           = struct('Group', sig_G, 'Condition', sig_C, 'Interaction', sig_I);
        RES.windows       = struct('Group', win_G, 'Condition', win_C, 'Interaction', win_I);
        RES.N             = struct('Young', size(YL,1), 'Old', size(OL,1));
        RES.clusters      = struct('Group', clu_G, 'Condition', clu_C, 'Interaction', clu_I);

        out_stats = fullfile(out_dir, sprintf( ...
            'BROADNESS_Timewise_ANOVA_ClusterPerm_Split%d_PC%d.mat', ...
            split_nums(hs), pcNumber));

        save(out_stats, 'RES', '-v7.3');

        sheet_name = sprintf('Split%d_PC%d', split_nums(hs), pcNumber);
        export_stats_to_xlsx(RES, out_xlsx, sheet_name);

        maskG_plot = sig_to_plot_mask(sig_G, t_stat, t_plot);
        maskC_plot = sig_to_plot_mask(sig_C, t_stat, t_plot);
        maskI_plot = sig_to_plot_mask(sig_I, t_stat, t_plot);

        wins_G_plot = mask_to_windows_on_grid(maskG_plot, t_plot, dt_plot, xlim_target);
        wins_C_plot = mask_to_windows_on_grid(maskC_plot, t_plot, dt_plot, xlim_target);
        wins_I_plot = mask_to_windows_on_grid(maskI_plot, t_plot, dt_plot, xlim_target);

        fig = figure('Visible','off', ...
                     'Color','w', ...
                     'Units','pixels', ...
                     'Position',[100 100 400 750]);

        set(fig, ...
            'DefaultAxesFontName', 'Helvetica', ...
            'DefaultTextFontName', 'Helvetica', ...
            'DefaultLegendFontName', 'Helvetica', ...
            'DefaultAxesFontSize', 8, ...
            'DefaultTextFontSize', 8, ...
            'DefaultLegendFontSize', 8, ...
            'DefaultAxesTitleFontSizeMultiplier', 1, ...
            'DefaultAxesLabelFontSizeMultiplier', 1);

        tlo = tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

        ax = nexttile(tlo, 1, [1 2]);
        hold(ax, 'on');

        xlim(ax, xlim_target);
        ylim(ax, [-1200 1000]);

        draw_gray_windows(ax, wins_G_plot, -1200, 1000, shade_gray, shade_alpha);
        draw_gray_windows(ax, wins_C_plot, -1200, 1000, shade_gray, shade_alpha);
        draw_gray_windows(ax, wins_I_plot, -1200, 1000, shade_gray, shade_alpha);

        plot_order = { ...
            young_idx, 2, combo_colors(1,:), labels{1}; ...
            older_idx, 2, combo_colors(2,:), labels{2}; ...
            young_idx, 1, combo_colors(3,:), labels{3}; ...
            older_idx, 1, combo_colors(4,:), labels{4}};

        h_lines = gobjects(1, size(plot_order,1));

        for j = 1:size(plot_order,1)

            subj_idx = plot_order{j,1};
            c        = plot_order{j,2};
            col      = plot_order{j,3};

            X_group = squeeze(X_plot(:, c, subj_idx));

            mu = mean(X_group, 2, 'omitnan');
            se = std(X_group, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(X_group), 2));

            patch(ax, ...
                [t_plot(:); flipud(t_plot(:))], ...
                [mu - se; flipud(mu + se)], ...
                col, ...
                'FaceAlpha', ribbon_alpha, ...
                'EdgeColor', 'none', ...
                'HandleVisibility', 'off', ...
                'Clipping', 'on');

            h_lines(j) = plot(ax, t_plot, mu, ...
                'LineWidth', 1.8, ...
                'Color', col);
        end

        draw_bottom_lines(ax, {wins_G_plot, wins_C_plot, wins_I_plot}, ...
            {clr_group, clr_cond, clr_int}, xlim_target, -1200, 1000);

        xline(ax, 0, ':', ...
            'Color', [0 0 0], ...
            'LineWidth', 0.75, ...
            'HandleVisibility','off');

        xlabel(ax, 'Time (s)');
        ylabel(ax, sprintf('BROADNESS activation, PC%d', pcNumber));

        xlim(ax, xlim_target);
        ylim(ax, [-1200 1000]);
        xticks(ax, xtick_vals);

        grid(ax, 'on');
        box(ax, 'on');

        axL = nexttile(tlo, 3, [1 2]);
        cla(axL);
        hold(axL, 'on');
        axis(axL, 'off');

        h_leg = gobjects(1, numel(tags));
        for j = 1:numel(tags)
            h_leg(j) = plot(axL, NaN, NaN, '-', ...
                'LineWidth', 3, ...
                'Color', combo_colors(j,:));
        end

        h_eff = gobjects(1,3);
        h_eff(1) = plot(axL, NaN, NaN, '-', 'LineWidth', 3, 'Color', clr_group);
        h_eff(2) = plot(axL, NaN, NaN, '-', 'LineWidth', 3, 'Color', clr_cond);
        h_eff(3) = plot(axL, NaN, NaN, '-', 'LineWidth', 3, 'Color', clr_int);

        lgd = legend(axL, ...
            [h_leg h_eff], ...
            {'Younger - Local', ...
             'Older - Local', ...
             'Younger - Global', ...
             'Older - Global', ...
             'Group main effect', ...
             'Condition main effect', ...
             'Interaction effect'}, ...
            'NumColumns', 4, ...
            'Orientation', 'horizontal', ...
            'Location', 'northwest');

        set(lgd, 'Interpreter','none', 'Units','normalized', 'Location','none');
        set(lgd, 'Position', [0.08 0.08 0.84 0.14]);
        set(lgd, 'ItemTokenSize', [15, 12]);
        lgd.TextColor = [0 0 0];
        lgd.Color     = [1 1 1];
        lgd.EdgeColor = 'none';

        sgtitle(tlo, sprintf('Older vs Younger — Split %d — PC%d', split_nums(hs), pcNumber));

        out_png = fullfile(out_dir, sprintf( ...
            'BROADNESS_TimeSeries_Older_vs_Younger_Local_Global_Split%d_PC%d.png', ...
            split_nums(hs), pcNumber));

        out_fig = fullfile(out_dir, sprintf( ...
            'BROADNESS_TimeSeries_Older_vs_Younger_Local_Global_Split%d_PC%d.fig', ...
            split_nums(hs), pcNumber));

        print(fig, out_png, '-dpng', '-r300');
        savefig(fig, out_fig);

        fprintf('Saved stats: %s\n', out_stats);
        fprintf('Saved XLSX:  %s [%s]\n', out_xlsx, sheet_name);
        fprintf('Saved PNG:   %s\n', out_png);
        fprintf('Saved FIG:   %s\n', out_fig);

        close(fig);
    end
end
%% B: Plot variance explained for first 10 PCs, Split 1 and Split 2

% This section reloads each BROADNESS file, extracts BROADNESS.Variance_BrainNetworks,
% computes effective dimensionality, and plots the first 10 PCs separately for Split 1 and Split 2.

variance_out_dir = fullfile(data_dir, 'BROADNESS_Output', 'VarianceExplained_PC1_PC10');
if ~exist(variance_out_dir, 'dir')
    mkdir(variance_out_dir);
end

% Styling controls
legendFontSize = 10;
gridAlpha      = 0.15;
gridColor      = [0 0 0];

darkPurple = [88 24 124] / 255;

for hs = 1:numel(broadness_files)

    B = load(broadness_files{hs});

    if isfield(B, 'BROADNESS')
        BROADNESS = B.BROADNESS;
    else
        fn = fieldnames(B);
        BROADNESS = B.(fn{1});
    end

    if ~isfield(BROADNESS, 'Variance_BrainNetworks')
        error('Could not find BROADNESS.Variance_BrainNetworks in %s.', broadness_files{hs});
    end

    % Extract eigenvalue / variance vector
    lam = BROADNESS.Variance_BrainNetworks(:);

    % Enforce non-negativity
    lam(lam < 0) = 0;

    % Compute effective dimensionality
    ED = (sum(lam)^2) / max(sum(lam.^2), realmin);

    % Sanity checks
    fprintf('\n%s / Split %d\n', halfsplit_names{hs}, split_nums(hs));
    fprintf('Effective Dimensionality (ED): %.3f\n', ED);

    nComp = numel(lam);
    if ED < 1 || ED > nComp
        warning('ED = %.3f outside theoretical bounds [1, %d]. Check eigenvalues.', ED, nComp);
    end

    % First 10 PCs
    nPC = min(10, numel(lam));
    x   = 1:nPC;
    y   = lam(1:nPC);

    % Figure + axes
    fig = figure('Visible','off', ...
                 'Color','w', ...
                 'Units','pixels', ...
                 'Position',[100 100 500 400]);

    set(fig, ...
        'DefaultTextFontName','Helvetica', ...
        'DefaultAxesFontName','Helvetica', ...
        'DefaultLegendFontName','Helvetica');

    ax = axes('Parent', fig);
    hold(ax, 'on');

    grid(ax, 'on');
    ax.GridLineStyle = '-';
    ax.GridAlpha     = gridAlpha;
    ax.GridColor     = gridColor;

    hVar = plot(ax, x, y, '-o', ...
        'Color', darkPurple, ...
        'MarkerFaceColor', darkPurple, ...
        'MarkerEdgeColor', darkPurple, ...
        'LineWidth', 1.8, ...
        'MarkerSize', 6, ...
        'DisplayName', 'Variance Explained');

    hEff = xline(ax, ED, ':k', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Effective Dimensionality');

    xlabel(ax, 'Principal Component', 'FontName','Helvetica');
    ylabel(ax, 'Variance Explained',  'FontName','Helvetica');
    ylim(ax, [0, 60]);


    xlim(ax, [0.5, nPC + 0.5]);
    xticks(ax, 1:nPC);

    set(ax, ...
        'FontName','Helvetica', ...
        'Box','on', ...
        'LineWidth', 1.0, ...
        'TickDir','out');

    lgd = legend(ax, hEff, ...
        sprintf('Effective Dimensionality = %.3f', ED), ...
        'Location','northeast');

    set(lgd, ...
        'Box','on', ...
        'FontName','Helvetica', ...
        'FontSize', legendFontSize, ...
        'ItemTokenSize',[12 8]);

    hold(ax, 'off');

    % Save outputs
    out_png = fullfile(variance_out_dir, sprintf( ...
        'BROADNESS_VarianceExplained_First10PCs_Split%d.png', split_nums(hs)));

    out_fig = fullfile(variance_out_dir, sprintf( ...
        'BROADNESS_VarianceExplained_First10PCs_Split%d.fig', split_nums(hs)));

    out_mat = fullfile(variance_out_dir, sprintf( ...
        'BROADNESS_VarianceExplained_First10PCs_Split%d_ED.mat', split_nums(hs)));

    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    VAR_RES = struct();
    VAR_RES.HalfSplit = halfsplit_names{hs};
    VAR_RES.Split     = split_nums(hs);
    VAR_RES.Variance_BrainNetworks = lam;
    VAR_RES.First10PCs = y;
    VAR_RES.ED = ED;

    save(out_mat, 'VAR_RES', '-v7.3');

    fprintf('Saved PNG: %s\n', out_png);
    fprintf('Saved FIG: %s\n', out_fig);
    fprintf('Saved MAT: %s\n', out_mat);

    close(fig);
end

%% C: Generate .nii files for PC1-PC3, Split 1 and Split 2

clear; clc;

data_dir = '/mainpath/';

matFiles = { ...
    fullfile(data_dir, 'BROADNESS_Output_Half1.mat'), ...
    fullfile(data_dir, 'BROADNESS_Output_Half2.mat')};

split_nums = [1 2];
pc_list = 1:3;

coordFile = '/mainpath/BROADNESS_MEG_AuditoryRecognition-main/BROADNESS_Toolbox/BROADNESS_External/MNI152_8mm_coord_dyi.mat';

refNii = '/mainpath/BROADNESS_MEG_AuditoryRecognition-main/BROADNESS_Toolbox/BROADNESS_External/MNI152_8mm_brain_diy.nii.gz';

outNiiDir = fullfile(data_dir, 'BROADNESS_Output', 'NIfTI_PC1_PC2_PC3');
if ~exist(outNiiDir, 'dir')
    mkdir(outNiiDir);
end

dataField = 'ActivationPatterns_BrainNetworks';

C = load(coordFile);

if ~isfield(C, 'MNI8')
    error('Could not find variable MNI8 in %s', coordFile);
end

xyz = C.MNI8;

info = niftiinfo(refNii);

for hs = 1:numel(matFiles)

    S = load(matFiles{hs});

    if ~isfield(S, 'BROADNESS')
        error('Could not find variable BROADNESS in %s', matFiles{hs});
    end

    B = S.BROADNESS;

    if ~isfield(B, dataField)
        error('Could not find BROADNESS.%s in %s', dataField, matFiles{hs});
    end

    X = B.(dataField);

    fprintf('\nSplit %d data size: %d sources x %d PCs\n', ...
        split_nums(hs), size(X,1), size(X,2));

    if size(xyz,1) ~= size(X,1) || size(xyz,2) ~= 3
        error('MNI8 must be %d x 3, but it is %d x %d.', ...
            size(X,1), size(xyz,1), size(xyz,2));
    end

    T = info.Transform.T;

    vox = [xyz, ones(size(xyz,1),1)] / T;
    vox = round(vox(:,1:3));

    for pcNumber = pc_list

        if pcNumber > size(X,2)
            error('Requested PC%d, but only %d PCs exist in %s.', ...
                pcNumber, size(X,2), matFiles{hs});
        end

        pcValues = X(:, pcNumber);

        % Flip ONLY PC3 in Split 2 before writing the .nii file.
        if split_nums(hs) == 2 && pcNumber == 3
            pcValues = -pcValues;
            fprintf('Flipped sign for Split 2, PC3 only.\n');
        end

        vol = zeros(info.ImageSize, 'single');

        valid = ...
            vox(:,1) >= 1 & vox(:,1) <= info.ImageSize(1) & ...
            vox(:,2) >= 1 & vox(:,2) <= info.ImageSize(2) & ...
            vox(:,3) >= 1 & vox(:,3) <= info.ImageSize(3) & ...
            ~isnan(pcValues) & ...
            ~isinf(pcValues);

        fprintf('Split %d, PC%d: valid sources written: %d / %d\n', ...
            split_nums(hs), pcNumber, sum(valid), numel(pcValues));

        if sum(valid) == 0
            error('No valid voxels found for Split %d, PC%d. Coordinate system may not match reference NIfTI.', ...
                split_nums(hs), pcNumber);
        end

        voxValid = vox(valid,:);
        pcValid = pcValues(valid);

        linIdx = sub2ind(info.ImageSize, ...
            voxValid(:,1), voxValid(:,2), voxValid(:,3));

        [uniqueIdx, ~, groupIdx] = unique(linIdx);
        avgVals = accumarray(groupIdx, pcValid, [], @mean);

        vol(uniqueIdx) = single(avgVals);

        infoOut = info;
        infoOut.Datatype = 'single';
        infoOut.BitsPerPixel = 32;

        outNii = fullfile(outNiiDir, sprintf( ...
            'BROADNESS_PC%d_Split%d.nii', pcNumber, split_nums(hs)));

        niftiwrite(vol, outNii, infoOut, 'Compressed', false);

        fprintf('Saved NIfTI: %s\n', outNii);
    end
end

%% ===================== LOCAL HELPERS =====================

function export_stats_to_xlsx(RES, out_xlsx, sheet_name)

    nT = numel(RES.time);

    group_windows       = windows_to_string(RES.windows.Group);
    condition_windows   = windows_to_string(RES.windows.Condition);
    interaction_windows = windows_to_string(RES.windows.Interaction);

    group_cluster_idx       = cluster_idx_to_string(RES.clusters.Group.idx_list, RES.time);
    condition_cluster_idx   = cluster_idx_to_string(RES.clusters.Condition.idx_list, RES.time);
    interaction_cluster_idx = cluster_idx_to_string(RES.clusters.Interaction.idx_list, RES.time);

    T = table( ...
        repmat(string(RES.HalfSplit), nT, 1), ...
        repmat(RES.Split, nT, 1), ...
        repmat(RES.PC, nT, 1), ...
        RES.time(:), ...
        RES.obs.Group(:), ...
        logical(RES.sig.Group(:)), ...
        repmat(RES.thresholds.Group.cluster_forming, nT, 1), ...
        repmat(RES.thresholds.Group.cluster_mass, nT, 1), ...
        repmat(string(group_windows), nT, 1), ...
        repmat(string(group_cluster_idx), nT, 1), ...
        repmat(string(numvec_to_string(RES.clusters.Group.masses)), nT, 1), ...
        repmat(string(logvec_to_string(RES.clusters.Group.sig_flags)), nT, 1), ...
        RES.obs.Condition(:), ...
        logical(RES.sig.Condition(:)), ...
        repmat(RES.thresholds.Condition.cluster_forming, nT, 1), ...
        repmat(RES.thresholds.Condition.cluster_mass, nT, 1), ...
        repmat(string(condition_windows), nT, 1), ...
        repmat(string(condition_cluster_idx), nT, 1), ...
        repmat(string(numvec_to_string(RES.clusters.Condition.masses)), nT, 1), ...
        repmat(string(logvec_to_string(RES.clusters.Condition.sig_flags)), nT, 1), ...
        RES.obs.Interaction(:), ...
        logical(RES.sig.Interaction(:)), ...
        repmat(RES.thresholds.Interaction.cluster_forming, nT, 1), ...
        repmat(RES.thresholds.Interaction.cluster_mass, nT, 1), ...
        repmat(string(interaction_windows), nT, 1), ...
        repmat(string(interaction_cluster_idx), nT, 1), ...
        repmat(string(numvec_to_string(RES.clusters.Interaction.masses)), nT, 1), ...
        repmat(string(logvec_to_string(RES.clusters.Interaction.sig_flags)), nT, 1), ...
        repmat(RES.N.Young, nT, 1), ...
        repmat(RES.N.Old, nT, 1), ...
        repmat(RES.alpha, nT, 1), ...
        repmat(RES.cluster_alpha, nT, 1), ...
        repmat(RES.n_perm, nT, 1), ...
        'VariableNames', { ...
            'HalfSplit', ...
            'Split', ...
            'PC', ...
            'Time_s', ...
            'Group_ObservedDifference', ...
            'Group_Significant', ...
            'Group_ClusterFormingThreshold', ...
            'Group_CriticalClusterMass', ...
            'Group_SignificantWindows_s', ...
            'Group_AllClusterWindows_s', ...
            'Group_AllClusterMasses', ...
            'Group_AllClusterSignificantFlags', ...
            'Condition_ObservedDifference', ...
            'Condition_Significant', ...
            'Condition_ClusterFormingThreshold', ...
            'Condition_CriticalClusterMass', ...
            'Condition_SignificantWindows_s', ...
            'Condition_AllClusterWindows_s', ...
            'Condition_AllClusterMasses', ...
            'Condition_AllClusterSignificantFlags', ...
            'Interaction_ObservedDifference', ...
            'Interaction_Significant', ...
            'Interaction_ClusterFormingThreshold', ...
            'Interaction_CriticalClusterMass', ...
            'Interaction_SignificantWindows_s', ...
            'Interaction_AllClusterWindows_s', ...
            'Interaction_AllClusterMasses', ...
            'Interaction_AllClusterSignificantFlags', ...
            'N_Young', ...
            'N_Old', ...
            'Alpha', ...
            'ClusterAlpha', ...
            'NPerm'});

    writetable(T, out_xlsx, 'Sheet', sheet_name);
end

function s = windows_to_string(W)

    if isempty(W)
        s = "";
        return;
    end

    parts = strings(size(W,1),1);
    for i = 1:size(W,1)
        parts(i) = sprintf('%.6f to %.6f', W(i,1), W(i,2));
    end

    s = strjoin(parts, '; ');
end

function s = cluster_idx_to_string(idx_list, tvec)

    if isempty(idx_list)
        s = "";
        return;
    end

    parts = strings(numel(idx_list),1);
    for i = 1:numel(idx_list)
        idx = idx_list{i};
        parts(i) = sprintf('%.6f to %.6f', tvec(idx(1)), tvec(idx(end)));
    end

    s = strjoin(parts, '; ');
end

function s = numvec_to_string(x)

    if isempty(x)
        s = "";
        return;
    end

    x = x(:)';
    parts = strings(1,numel(x));
    for i = 1:numel(x)
        parts(i) = sprintf('%.10g', x(i));
    end

    s = strjoin(parts, '; ');
end

function s = logvec_to_string(x)

    if isempty(x)
        s = "";
        return;
    end

    x = logical(x(:)');
    parts = strings(1,numel(x));
    for i = 1:numel(x)
        parts(i) = string(x(i));
    end

    s = strjoin(parts, '; ');
end

function [sig_mask, thr, clusters] = perm_between_time(A, B, nperm, c_alpha, fwer_alpha)

    A = A(~all(isnan(A),2), :);
    B = B(~all(isnan(B),2), :);

    Na = size(A,1);
    Nb = size(B,1);
    T  = size(A,2);

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

    D = D(~all(isnan(D),2), :);

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
        idx_list{end+1} = seg; %#ok<AGROW>
        masses(end+1) = sum(abs_series(seg)); %#ok<AGROW>
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
        W(end+1, :) = [tvec(seg(1)), tvec(seg(end))]; %#ok<AGROW>
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

function mask_plot = sig_to_plot_mask(sig_mask, t_stat, t_plot)

    t_stat = t_stat(:);
    t_plot = t_plot(:);
    sig_mask = logical(sig_mask(:))';

    if isempty(t_stat) || isempty(t_plot) || numel(sig_mask) ~= numel(t_stat)
        mask_plot = false(size(t_plot));
        return;
    end

    dt_plot = median(diff(t_plot));
    mask_plot = false(size(t_plot));

    for k = 1:numel(t_stat)
        if ~sig_mask(k)
            continue;
        end

        [d, j] = min(abs(t_plot - t_stat(k)));

        if d <= dt_plot/2 + 10*eps(max(abs(t_plot)))
            mask_plot(j) = true;
        end
    end

    mask_plot = mask_plot(:)';
end

function wins = mask_to_windows_on_grid(mask_plot, t_plot, dt_plot, xlim_ok)

    wins = [];

    if ~any(mask_plot)
        return;
    end

    on = find(mask_plot(:)');
    br = [1, find(diff(on) > 1)+1, numel(on)+1];

    for b = 1:numel(br)-1
        seg = on(br(b):br(b+1)-1);
        t1 = t_plot(seg(1)) - dt_plot/2;
        t2 = t_plot(seg(end)) + dt_plot/2;

        t1 = max(t1, xlim_ok(1));
        t2 = min(t2, xlim_ok(2));

        if t2 > t1
            wins(end+1,:) = [t1, t2]; %#ok<AGROW>
        end
    end
end

function draw_gray_windows(ax, wins, ymin, ymax, shade_gray, shade_alpha)

    if isempty(wins)
        return;
    end

    for i = 1:size(wins,1)
        t1 = wins(i,1);
        t2 = wins(i,2);

        patch('XData', [t1 t2 t2 t1], ...
              'YData', [ymin ymin ymax ymax], ...
              'FaceColor', shade_gray, ...
              'FaceAlpha', shade_alpha, ...
              'EdgeColor', 'none', ...
              'Parent', ax, ...
              'HitTest', 'off', ...
              'Clipping', 'on');
    end
end

function draw_bottom_lines(ax, wins_cell, effect_colors, xlim_ok, ymin, ymax)

    oldUnits = get(ax,'Units');
    set(ax,'Units','pixels');
    axPos = get(ax,'Position');
    set(ax,'Units',oldUnits);

    yRange = ymax - ymin;
    spacing_factor = 5;
    dy_data = (axPos(4) > 0) * (spacing_factor / max(axPos(4),1)) * yRange;

    if dy_data == 0
        dy_data = 0.003 * yRange;
    end

    y_base = ymin + 0.03 * yRange;
    n_levels = 3;

    for lvl = 1:3
        wins = wins_cell{lvl};

        if isempty(wins)
            continue;
        end

        y_line = y_base + (n_levels - lvl) * dy_data;

        for i = 1:size(wins,1)
            t1 = max(wins(i,1), xlim_ok(1));
            t2 = min(wins(i,2), xlim_ok(2));

            if t2 <= t1
                continue;
            end

            plot(ax, [t1 t2], [y_line y_line], '-', ...
                'LineWidth', 3.0, ...
                'Color', effect_colors{lvl}, ...
                'HitTest', 'off');
        end
    end
end

%% D: Correlations between halfsplits and all data used in the main analysis
% Export SAP and TS from correlations.mat to Excel

clear; clc;

% Paths
infile  = '/mmainpath/correlations.mat';
outdir  = '/mainpath/';
outfile = fullfile(outdir, 'correlations_export.xlsx');

% Load data
S = load(infile);

SAP = S.SAP;   % expected: 11 x 3
TS  = S.TS;    % expected: 11 x 3 x 2

% Basic checks
assert(size(SAP,2) == 3, 'SAP should have 3 columns.');
assert(size(TS,2) == 3, 'TS should have 3 network columns.');
assert(size(TS,3) == 2, 'TS should have 2 conditions in the 3rd dimension.');

nRows = size(SAP,1);

% Row labels
Split = strings(nRows,1);
for i = 1:nRows
    if i == nRows
        Split(i) = "Row_" + i + "_possibly_main_or_fixed_split";
    else
        Split(i) = "Row_" + i;
    end
end

% SAP table
SAP_table = table( ...
    Split, ...
    SAP(:,1), SAP(:,2), SAP(:,3), ...
    'VariableNames', {'Split','SAP_PC1','SAP_PC2','SAP_PC3'} ...
);

% TS Global table: TS(:,:,1)
TS_Global_table = table( ...
    Split, ...
    TS(:,1,1), TS(:,2,1), TS(:,3,1), ...
    'VariableNames', {'Split','TS_Global_PC1','TS_Global_PC2','TS_Global_PC3'} ...
);

% TS Local table: TS(:,:,2)
TS_Local_table = table( ...
    Split, ...
    TS(:,1,2), TS(:,2,2), TS(:,3,2), ...
    'VariableNames', {'Split','TS_Local_PC1','TS_Local_PC2','TS_Local_PC3'} ...
);

% Combined wide table
Combined_table = table( ...
    Split, ...
    SAP(:,1), SAP(:,2), SAP(:,3), ...
    TS(:,1,1), TS(:,2,1), TS(:,3,1), ...
    TS(:,1,2), TS(:,2,2), TS(:,3,2), ...
    'VariableNames', { ...
        'Split', ...
        'SAP_PC1','SAP_PC2','SAP_PC3', ...
        'TS_Global_PC1','TS_Global_PC2','TS_Global_PC3', ...
        'TS_Local_PC1','TS_Local_PC2','TS_Local_PC3'} ...
);

% Long-format TS table
rows = {};
cnt = 0;

condition_names = {'Global','Local'};

for rr = 1:size(TS,1)
    for pc = 1:3
        for cc = 1:2
            cnt = cnt + 1;
            rows(cnt,:) = { ...
                char(Split(rr)), ...
                pc, ...
                condition_names{cc}, ...
                TS(rr,pc,cc) ...
            };
        end
    end
end

TS_Long_table = cell2table(rows, ...
    'VariableNames', {'Split','PC','Condition','Correlation'} ...
);

% Write Excel file
if exist(outfile, 'file')
    delete(outfile);
end

writetable(Combined_table,  outfile, 'Sheet', 'Combined');
writetable(SAP_table,       outfile, 'Sheet', 'SAP');
writetable(TS_Global_table, outfile, 'Sheet', 'TS_Global');
writetable(TS_Local_table,  outfile, 'Sheet', 'TS_Local');
writetable(TS_Long_table,   outfile, 'Sheet', 'TS_Long');

fprintf('Saved Excel file to:\n%s\n', outfile);