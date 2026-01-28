%% Analysis of data seperated by group and condition (for figure 6)

%% 0) STARTUP

% Simply download the BROADNESS Toolbox folder and place it in the working directory,
% making sure not to alter the structure of its functions, subfolders, or files.

clear 
close all
clc

% Setup directories
path_home = '/main_path/BROADNESS_MEG_AuditoryRecognition-main/BROADNESS_Toolbox';
addpath(path_home)
BROADNESS_Startup(path_home);
addpath(fullfile(matlabroot,'toolbox','stats','stats'),'-begin') % makes sure the pca function is the standard function in matlab

%%

%% 1a) PERFORM BROADNESS (ONLY ESSENTIAL INPUTS)

%%% ------------------- USER SETTINGS ------------------- %%%

load('/main_path/Data/Mathias_MMN/MMNSubtracted_Average_SignFixed.mat');
DATA = dum;

% Remove first participant and define time vector

% 100 ms baseline
time = -0.100:0.004:0.8;
DATA = DATA(:,101:326,:,2:78); 

size(DATA)

% Load groups and adjust subid after removing first participant
load('/main_path/Data/Mathias_MMN/groups.mat'); % older, young
older = older - 1;
young = young - 1;

older_subj = older(:);   % row-wise versions
young_subj = young(:);

%% 1b) Split data into 4 different groups

cond_global = 1;
cond_local  = 2;

% Divide the struct
DATA_young_local = DATA(:,:,cond_local,young_subj);
DATA_older_local = DATA(:,:,cond_local,older_subj);
DATA_young_global = DATA(:,:,cond_global,young_subj);
DATA_older_global = DATA(:,:,cond_global,older_subj);

% Check sizes
size(DATA_young_local)
size(DATA_older_local)
size(DATA_young_global)
size(DATA_older_global)

% RUN BROADNESS
BROADNESS_young_local = BROADNESS_NetworkEstimation(DATA_young_local,time);
BROADNESS_older_local = BROADNESS_NetworkEstimation(DATA_older_local,time);
BROADNESS_young_global = BROADNESS_NetworkEstimation(DATA_young_global,time);
BROADNESS_older_global = BROADNESS_NetworkEstimation(DATA_older_global,time);

%% Force them back to 4D
[n1, n2, n3] = size(BROADNESS_young_local.TimeSeries_BrainNetworks);   % n1=226, n2=225, n3=40
BROADNESS_young_local.TimeSeries_BrainNetworks = reshape(BROADNESS_young_local.TimeSeries_BrainNetworks, n1, n2, 1, n3);
BROADNESS_young_local.OriginalData = reshape(BROADNESS_young_local.OriginalData, 3559, 226, 1, 37);

[nn1, nn2, nn3] = size(BROADNESS_older_local.TimeSeries_BrainNetworks);   % n1=226, n2=225, n3=40
BROADNESS_older_local.TimeSeries_BrainNetworks = reshape(BROADNESS_older_local.TimeSeries_BrainNetworks, nn1, nn2, 1, nn3);
BROADNESS_older_local.OriginalData = reshape(BROADNESS_older_local.OriginalData, 3559, 226, 1, 40);

[nnn1, nnn2, nnn3] = size(BROADNESS_young_global.TimeSeries_BrainNetworks);   % n1=226, n2=225, n3=40
BROADNESS_young_global.TimeSeries_BrainNetworks = reshape(BROADNESS_young_global.TimeSeries_BrainNetworks, nnn1, nnn2, 1, nnn3);
BROADNESS_young_global.OriginalData = reshape(BROADNESS_young_global.OriginalData, 3559, 226, 1, 37);

[nnnn1, nnnn2, nnnn3] = size(BROADNESS_older_global.TimeSeries_BrainNetworks);   % n1=226, n2=225, n3=40
BROADNESS_older_global.TimeSeries_BrainNetworks = reshape(BROADNESS_older_global.TimeSeries_BrainNetworks, nnnn1, nnnn2, 1, nnnn3);
BROADNESS_older_global.OriginalData = reshape(BROADNESS_older_global.OriginalData, 3559, 226, 1, 40);

%% STD BROADNESS

% Run BROADNESS network estimation (default parameters)
BROADNESS = BROADNESS_NetworkEstimation(DATA,time);

%% 1c) %% Calculate ED on BROADNESS, PCA applied to all data

% Extract eigenvalue / variance vector

% Run 1 of these 4 lines in seperate

lam = BROADNESS_young_local.Variance_BrainNetworks(:);
%lam = BROADNESS_young_global.Variance_BrainNetworks(:);
%lam = BROADNESS_older_local.Variance_BrainNetworks(:);
%lam = BROADNESS_older_global.Variance_BrainNetworks(:);

% Enforce non-negativity (ED is undefined otherwise)
lam(lam < 0) = 0;

% Compute effective dimensionality
ED = (sum(lam)^2) / max(sum(lam.^2), realmin);

% ---- Sanity checks ----
fprintf('Effective Dimensionality (ED): %.3f\n', ED);

% Theoretical bounds
nComp = numel(lam);
if ED < 1 || ED > nComp
    warning('ED = %.3f outside theoretical bounds [1, %d]. Check eigenvalues.', ED, nComp);
end

% ED young local = 2.682
% ED young global = 2.835
% ED older local = 1.986
% ED older global = 3.123

%% 1d) Plot variance explained 

% Plot variance explained for first 10 PCs from:
%   BROADNESS.Variance_BrainNetworks

% Styling:
% - dark purple line + markers
% - black dotted vertical line at x = 2.835
% - Helvetica everywhere
% - boxed legend inside plot

% -----------------------------
% Data prep (robust to row/col)
% -----------------------------

% Run 1 of these 4 lines in seperate
ve = BROADNESS_young_local.Variance_BrainNetworks(:);   %xEff = 2.682; % young local
%ve = BROADNESS_older_local.Variance_BrainNetworks(:);   %xEff = 1.986; % older local
%ve = BROADNESS_young_global.Variance_BrainNetworks(:);   %xEff = 2.835; % young global
%ve = BROADNESS_older_global.Variance_BrainNetworks(:);   %xEff = 3.123; % older global

nPC = min(10, numel(ve));
x   = 1:nPC;
y   = ve(1:nPC);

% If values are proportions and one wants percentages:
% y = 100*y;

% -----------------------------
% Define dark purple color
% -----------------------------
darkPurple = [88 24 124] / 255;   % perceptually dark, print-safe purple

% -----------------------------
% Figure + axes (Helvetica)
% -----------------------------
fig = figure('Color','w');
ax  = axes('Parent',fig);
hold(ax,'on');

set(fig, 'DefaultTextFontName','Helvetica');
set(fig, 'DefaultAxesFontName','Helvetica');

% -----------------------------
% Plot variance explained
% -----------------------------
hVar = plot(ax, x, y, '-o', ...
    'Color', darkPurple, ...
    'MarkerFaceColor', darkPurple, ...
    'MarkerEdgeColor', darkPurple, ...
    'LineWidth', 1.8, ...
    'MarkerSize', 6, ...
    'DisplayName', 'Variance Explained');

% -----------------------------
% Vertical dotted line
% -----------------------------
%hEff = xline(ax, xEff, ':k', ...
%    'LineWidth', 1.5, ...
%    'DisplayName', 'Effective Dimensionality');

% -----------------------------
% Labels / axes styling
% -----------------------------
xlabel(ax, 'Principal Component', 'FontName','Helvetica');
ylabel(ax, 'Variance Explained',  'FontName','Helvetica');

xlim(ax, [0.5, nPC + 0.5]);
ylim(ax, [0 70]);
xticks(ax, 1:nPC);

set(ax, ...
    'FontName','Helvetica', ...
    'Box','on', ...
    'LineWidth', 1.0, ...
    'TickDir','out');

grid(ax, 'on');
ax.GridLineStyle  = ':';
ax.GridAlpha      = 0.25;
ax.MinorGridAlpha = 0.15;

% Optional y-limits (disable if inappropriate for data)
% ylim(ax, [min(0, min(y)*0.95), max(y)*1.05]);

% -----------------------------
% Legend (boxed, inside plot)
% -----------------------------
lgd = legend(ax, hEff, 'Effective Dimensionality', ...
    'Location','northeast');
set(lgd, ...
    'Box','on', ...
    'FontName','Helvetica');

hold(ax,'off');


%% 2a) BROADNESS VISUALIZATION

% NOTE: YOU ONLY NEED TO RUN 2a) or 2b), as 2b) simply has more advanced
% settings.

%This section generates 5 plots, described as follows:
%   #1) Dynamic brain activity map of the original data
%   #2) Variance explained by the networks
%   #3) Time series of the networks
%   #4) Activation patterns of the networks (3D)
%   #5) Activation patterns of the networks (nifti images)

%%% ------------------- USER SETTINGS ------------------- %%%

% Common options
Options = [];
Options.name_nii = '/main_path/Output/';
load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); 
Options.MNI_coords = MNI8;

% Run each of the 5 lines below in seperate

% All data
%Options.Labels = {'All data'}; BROADNESS_Visualizer(BROADNESS, Options);

% Group/condition specific
%Options.Labels = {'Younger adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_young_local, Options_yl);
%Options.Labels = {'Older adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_older_local, Options); 
%Options.Labels = {'Younger adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_young_global, Options);
Options.Labels = {'Older adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_older_global, Options);

%% 2b - i) Flip sign of a single NIfTI - Older global - as it has been reversed via PCA (multiply voxel intensities by -1)

% ONLY DO THIS FOR OLDER GLOBAL PC #1 FILES

inFile  = '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/PCA_ActivationPattern_BrainNetwork_#1.nii';
outFile = '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/PCA_ActivationPattern_BrainNetwork_#1_FLIPPED.nii';

% Load
nii = load_nii(inFile);

% Ensure numeric
img = double(nii.img);

% Flip sign
img = -img;

% Put back (keep header/affine exactly as-is)
nii.img = img;

% Save
save_nii(nii, outFile);

fprintf('Saved flipped NIfTI:\n  %s\n', outFile);


%% 2b - ii) ALSO flip sign of Activation values in the corresponding Excel table (5th column)

excelFile = '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/ActivationTable_BrainNetwork_1.xlsx';

% Read table
T = readtable(excelFile);

% Flip sign of Activation column (5th column)
T{:,5} = -T{:,5};

% Write back (overwrite same file)
writetable(T, excelFile);

fprintf('Flipped Activation values in Excel table:\n  %s\n', excelFile);


%% 2c) - Convert MNI coords into brain labels for positively and negatively contributing voxels independently
%
% For each network ActivationTable_*.xlsx:
%   1) split voxels into positive vs negative contributions
%   2) map each voxel's MNI (mm) coordinate to an atlas label
%   3) count voxels per label and output sorted lists
%   4) ALSO: if a voxel lands on background (0) or unlabeled, assign it the
%            nearest label within a growing voxel-radius (0..maxSearchRadius)
%   5) Output per-area counts split by the radius used (exact, 1-away, 2-away, ...)
%

% ----------------------- USER SETTINGS -----------------------
networkFiles = {
    '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/ActivationTable_BrainNetwork_1.xlsx'
    '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/ActivationTable_BrainNetwork_2.xlsx'
    '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/ActivationTable_BrainNetwork_3.xlsx'
};

% Choose an atlas in MNI space (NIfTI) + a label lookup (CSV/TXT).
atlasNiftiPath  = '/main_path/AAL3/AAL3v1.nii';
atlasLabelsPath = '/main_path/AAL3/AAL3v1.nii.txt';

% If you KNOW your contribution column name, set it here (recommended).
% Otherwise leave as "" to let the script guess.
contributionColumnOverride = "";  % e.g., "Contribution", "Weight", "Loading", "Beta", "Z", etc.

% What to do with atlas voxels that map to 0 / background / unknown
unknownLabelName = "Unknown/Background";

% How far (in ATLAS VOXELS) to search for nearest non-zero label
maxSearchRadius = 3;   % 0 = exact only; 1 = allow shell radius 1; etc.

% Output folder (will be created if missing)
outDir = fullfile('/main_path/', 'voxel_counts_by_area');
if ~exist(outDir, 'dir'); mkdir(outDir); end
% ------------------------------------------------------------

% Load atlas volume
assert(exist(atlasNiftiPath, 'file')==2, 'Atlas NIfTI not found: %s', atlasNiftiPath);
V = spm_vol(atlasNiftiPath);
A = spm_read_vols(V);

% Load atlas labels LUT (robust for AAL-style txt)
lut = loadAtlasLUT(atlasLabelsPath, unknownLabelName);

% Collect results for all PCs (so we can write ONE combined file)
results = struct([]);

% Process each network file
for iF = 1:numel(networkFiles)
    xlsxPath = networkFiles{iF};
    assert(exist(xlsxPath, 'file')==2, 'Network table not found: %s', xlsxPath);

    T = readtable(xlsxPath);

    % Find MNI coordinate columns + contribution column
    [xCol, yCol, zCol] = guessMNIColumns(T);
    contribCol = guessContributionColumn(T, xCol, yCol, zCol, contributionColumnOverride);

    XYZ = [T.(xCol), T.(yCol), T.(zCol)];
    contrib = T.(contribCol);

    % Basic sanity checks
    if ~isnumeric(XYZ) || size(XYZ,2)~=3
        error('MNI columns are not numeric 3D coords. Detected columns: %s, %s, %s', xCol, yCol, zCol);
    end
    if ~isnumeric(contrib)
        error('Contribution column is not numeric. Detected column: %s', contribCol);
    end

    % Map each MNI coord -> atlas label (WITH radius-based snapping)
    [labels, usedRadius] = mniToAtlasLabelWithSnap(XYZ, V, A, lut, unknownLabelName, maxSearchRadius);

    % Split positive/negative (zero excluded)
    posMask = contrib > 0;
    negMask = contrib < 0;

    % Count labels by distance radius used
    posCountsByR = countLabelsByRadius(labels(posMask), usedRadius(posMask), maxSearchRadius);
    negCountsByR = countLabelsByRadius(labels(negMask), usedRadius(negMask), maxSearchRadius);

    % Convert to tables (Area + Dist0..DistR + Total)
    posTableByR = countsByRadiusToTable(posCountsByR, maxSearchRadius, "Area");
    negTableByR = countsByRadiusToTable(negCountsByR, maxSearchRadius, "Area");

    % Also keep the old simple totals (optional)
    posSimple = countLabels(labels(posMask));
    negSimple = countLabels(labels(negMask));
    posTableSimple = countsToTable(posSimple, "Area", "VoxelCount");
    negTableSimple = countsToTable(negSimple, "Area", "VoxelCount");

    % Summary: how many voxels labeled at each radius overall (not per area)
    posRadiusSummary = radiusSummaryTable(usedRadius(posMask), maxSearchRadius, "Positive");
    negRadiusSummary = radiusSummaryTable(usedRadius(negMask), maxSearchRadius, "Negative");

    % Store results (instead of writing 3 separate files)
    [~, baseName, ~] = fileparts(xlsxPath);

    results(iF).baseName = baseName;
    results(iF).xCol = xCol; results(iF).yCol = yCol; results(iF).zCol = zCol; results(iF).contribCol = contribCol;
    results(iF).posTableByR = posTableByR;
    results(iF).negTableByR = negTableByR;
    results(iF).posRadiusSummary = posRadiusSummary;
    results(iF).negRadiusSummary = negRadiusSummary;
    results(iF).posTableSimple = posTableSimple;
    results(iF).negTableSimple = negTableSimple;

    % Print key info
    fprintf('\n=== %s ===\n', baseName);
    fprintf('MNI columns: %s %s %s | Contribution: %s\n', xCol, yCol, zCol, contribCol);
    fprintf('Atlas snapping max radius: %d voxels\n', maxSearchRadius);

    fprintf('\n+ Positive radius summary\n');
    disp(posRadiusSummary);

    fprintf('\n- Negative radius summary\n');
    disp(negRadiusSummary);

    fprintf('\n+ Positive areas by radius (top 15 by Total)\n');
    disp(posTableByR(1:min(15,height(posTableByR)),:));

    fprintf('\n- Negative areas by radius (top 15 by Total)\n');
    disp(negTableByR(1:min(15,height(negTableByR)),:));
end

% ---- Write ONE combined workbook with PC1/PC2/PC3 side-by-side in each sheet ----
outXlsx = fullfile(outDir, 'PC1_PC2_PC3_voxelCountsByArea.xlsx');

% Sheets to write (same content as before, but combined)
sheetNames = {
    'Positive_ByRadius'
    'Negative_ByRadius'
    'Positive_RadiusSummary'
    'Negative_RadiusSummary'
    'Positive_SimpleTotals'
    'Negative_SimpleTotals'
};

% Build and write each sheet as a combined cell array
for s = 1:numel(sheetNames)
    sh = sheetNames{s};

    blocks = cell(1, numel(results));
    for iF = 1:numel(results)
        pcLabel = sprintf('PC #%.0f', iF);

        switch sh
            case 'Positive_ByRadius'
                blocks{iF} = tableToPCBlockCell(results(iF).posTableByR, pcLabel);
            case 'Negative_ByRadius'
                blocks{iF} = tableToPCBlockCell(results(iF).negTableByR, pcLabel);
            case 'Positive_RadiusSummary'
                blocks{iF} = tableToPCBlockCell(results(iF).posRadiusSummary, pcLabel);
            case 'Negative_RadiusSummary'
                blocks{iF} = tableToPCBlockCell(results(iF).negRadiusSummary, pcLabel);
            case 'Positive_SimpleTotals'
                blocks{iF} = tableToPCBlockCell(results(iF).posTableSimple, pcLabel);
            case 'Negative_SimpleTotals'
                blocks{iF} = tableToPCBlockCell(results(iF).negTableSimple, pcLabel);
        end
    end

    combined = combinePCBlocksHorizontally(blocks);

    writecell(combined, outXlsx, 'Sheet', sh, 'WriteMode', 'overwritesheet');
end

fprintf('\nDone. Combined output saved to: %s\n', outXlsx);

%% ----------------------- FUNCTIONS -----------------------

function lut = loadAtlasLUT(labelsPath, unknownLabelName)
    % Robust LUT loader that supports:
    % 1) CSV with headers (id + name/label)
    % 2) AAL-style text where each line is whitespace-delimited, e.g.:
    %       1 Precentral_L 1
    %       2 Precentral_R 2
    %    and possibly some lines missing the final numeric id:
    %       35 Cingulate_Ant_L
    assert(exist(labelsPath,'file')==2, 'Atlas label file not found: %s', labelsPath);

    [~,~,ext] = fileparts(labelsPath);

    if strcmpi(ext,'.csv')
        L = readtable(labelsPath, 'TextType','string');
        varNames = lower(string(L.Properties.VariableNames));

        idIdx   = find(contains(varNames,"id") | contains(varNames,"index") | contains(varNames,"value"), 1);
        nameIdx = find(contains(varNames,"name") | contains(varNames,"label") | contains(varNames,"region") | contains(varNames,"structure"), 1);

        if isempty(idIdx) || isempty(nameIdx)
            error('Could not identify id/name columns in CSV LUT. Need columns like id,name.');
        end

        ids = L{:, idIdx};
        names = L{:, nameIdx};

        if ~isnumeric(ids)
            ids = str2double(string(ids));
        end
        if any(isnan(ids))
            error('Atlas LUT ids contain NaNs after parsing. Check your CSV LUT format.');
        end

        lut = containers.Map('KeyType','double','ValueType','char');
        for k = 1:numel(ids)
            lut(double(ids(k))) = char(names(k));
        end

        if ~isKey(lut, 0)
            lut(0) = char(unknownLabelName);
        end
        return;
    end

    lines = readlines(labelsPath);
    lut = containers.Map('KeyType','double','ValueType','char');

    for i = 1:numel(lines)
        s = strtrim(lines(i));
        if strlength(s)==0
            continue;
        end

        parts = regexp(s, '\s+', 'split');
        if numel(parts) < 2
            continue;
        end

        id1 = str2double(parts{1});
        if isnan(id1)
            continue;
        end

        idLast = str2double(parts{end});
        if ~isnan(idLast) && numel(parts) >= 3
            atlasId = idLast;
            nameTokens = parts(2:end-1);
        else
            atlasId = id1;
            nameTokens = parts(2:end);
        end

        if isempty(nameTokens)
            continue;
        end

        regionName = strjoin(nameTokens, ' ');
        lut(double(atlasId)) = char(regionName);
    end

    if ~isKey(lut, 0)
        lut(0) = char(unknownLabelName);
    end

    if lut.Count <= 1
        warning('Atlas LUT parsing produced very few entries. Check LUT file formatting: %s', labelsPath);
    end
end

function [xCol, yCol, zCol] = guessMNIColumns(T)
    vn = string(T.Properties.VariableNames);
    vnl = lower(vn);

    xCandidates = vn( ((contains(vnl,"mni") & contains(vnl,"x"))) | strcmp(vnl,"x") | endsWith(vnl,"_x") );
    yCandidates = vn( ((contains(vnl,"mni") & contains(vnl,"y"))) | strcmp(vnl,"y") | endsWith(vnl,"_y") );
    zCandidates = vn( ((contains(vnl,"mni") & contains(vnl,"z"))) | strcmp(vnl,"z") | endsWith(vnl,"_z") );

    if isempty(xCandidates), xCandidates = vn(strcmp(vnl,"x") | contains(vnl,"coordx") | contains(vnl,"x_mm") | contains(vnl,"xmm")); end
    if isempty(yCandidates), yCandidates = vn(strcmp(vnl,"y") | contains(vnl,"coordy") | contains(vnl,"y_mm") | contains(vnl,"ymm")); end
    if isempty(zCandidates), zCandidates = vn(strcmp(vnl,"z") | contains(vnl,"coordz") | contains(vnl,"z_mm") | contains(vnl,"zmm")); end

    if isempty(xCandidates) || isempty(yCandidates) || isempty(zCandidates)
        numMask = varfun(@isnumeric, T, 'OutputFormat','uniform');
        numCols = vn(numMask);
        if numel(numCols) < 3
            error('Could not find MNI columns. No obvious X/Y/Z or >=3 numeric columns.');
        end
        xCol = numCols(1); yCol = numCols(2); zCol = numCols(3);
        return;
    end

    xCol = xCandidates(1);
    yCol = yCandidates(1);
    zCol = zCandidates(1);
end

function contribCol = guessContributionColumn(T, xCol, yCol, zCol, override)
    if strlength(override) > 0
        assert(any(strcmp(string(T.Properties.VariableNames), override)), ...
            'Override contribution column not found: %s', override);
        contribCol = override;
        return;
    end

    vn = string(T.Properties.VariableNames);
    vnl = lower(vn);

    exclude = ismember(vn, [xCol,yCol,zCol]);

    preferred = contains(vnl,"contrib") | contains(vnl,"weight") | contains(vnl,"loading") | ...
                contains(vnl,"beta") | strcmp(vnl,"z") | contains(vnl,"zscore") | ...
                contains(vnl,"t") | contains(vnl,"stat") | contains(vnl,"value") | contains(vnl,"coef");

    idx = find(preferred & ~exclude, 1);
    if ~isempty(idx)
        contribCol = vn(idx);
        return;
    end

    numMask = varfun(@isnumeric, T, 'OutputFormat','uniform') & ~exclude;
    numCols = vn(numMask);
    if isempty(numCols)
        error('No numeric column found for contribution (besides coordinates).');
    end

    vars = zeros(numel(numCols),1);
    for k = 1:numel(numCols)
        x = T.(numCols(k));
        vars(k) = var(x(~isnan(x)));
    end
    [~,kmax] = max(vars);
    contribCol = numCols(kmax);
end

function [labels, usedRadius] = mniToAtlasLabelWithSnap(XYZmm, V, A, lut, unknownLabelName, maxR)
    % Returns:
    %  labels: assigned region name for each coordinate
    %  usedRadius:
    %    0..maxR = snapped using that voxel shell radius (0 means exact voxel non-zero label)
    %    -1      = could not find any non-zero label within maxR or out-of-bounds
    n = size(XYZmm,1);
    labels = strings(n,1);
    usedRadius = -1 * ones(n,1);

    M = V.mat;    % voxel->mm
    iM = inv(M);  % mm->voxel

    volSize = size(A);

    for r = 1:n
        mm = [XYZmm(r,:), 1]';
        vox = iM * mm;
        ijk0 = round(vox(1:3))';

        % If out of bounds, keep Unknown
        if any(isnan(ijk0)) || any(ijk0 < 1) || ...
                ijk0(1) > volSize(1) || ijk0(2) > volSize(2) || ijk0(3) > volSize(3)
            labels(r) = unknownLabelName;
            usedRadius(r) = -1;
            continue;
        end

        val0 = A(ijk0(1), ijk0(2), ijk0(3));
        if isnan(val0); val0 = 0; end
        key0 = double(round(val0));

        % Exact hit must be non-zero AND in LUT
        if key0 ~= 0 && isKey(lut, key0)
            labels(r) = string(lut(key0));
            usedRadius(r) = 0;
            continue;
        end

        % Otherwise, search shells r=1..maxR for nearest non-zero labels
        found = false;
        for rad = 1:maxR
            keysInShell = collectNonzeroKeysInShell(A, ijk0, rad);

            if ~isempty(keysInShell)
                % Choose the most frequent label key in this shell (mode)
                chosenKey = modeTiebreakSmallest(keysInShell);

                if isKey(lut, chosenKey)
                    labels(r) = string(lut(chosenKey));
                else
                    % If LUT is missing this key (unexpected), mark unknown
                    labels(r) = unknownLabelName;
                end
                usedRadius(r) = rad;
                found = true;
                break;
            end
        end

        if ~found
            labels(r) = unknownLabelName;
            usedRadius(r) = -1;
        end
    end
end

function keys = collectNonzeroKeysInShell(A, ijk0, rad)
    % Collect non-zero atlas integer keys from the Chebyshev shell at radius rad.
    % Shell means max(|dx|,|dy|,|dz|) == rad (not <= rad).
    volSize = size(A);

    xs = (ijk0(1)-rad):(ijk0(1)+rad);
    ys = (ijk0(2)-rad):(ijk0(2)+rad);
    zs = (ijk0(3)-rad):(ijk0(3)+rad);

    % Clip to bounds
    xs = xs(xs>=1 & xs<=volSize(1));
    ys = ys(ys>=1 & ys<=volSize(2));
    zs = zs(zs>=1 & zs<=volSize(3));

    vals = [];

    for xi = 1:numel(xs)
        for yi = 1:numel(ys)
            for zi = 1:numel(zs)
                dx = abs(xs(xi) - ijk0(1));
                dy = abs(ys(yi) - ijk0(2));
                dz = abs(zs(zi) - ijk0(3));

                if max([dx,dy,dz]) ~= rad
                    continue; % not in shell boundary
                end

                v = A(xs(xi), ys(yi), zs(zi));
                if isnan(v); v = 0; end
                k = double(round(v));
                if k ~= 0
                    vals(end+1,1) = k; %#ok<AGROW>
                end
            end
        end
    end

    keys = vals;
end

function chosen = modeTiebreakSmallest(x)
    % Returns the most frequent value in x; ties broken by choosing smallest key.
    if isempty(x)
        chosen = [];
        return;
    end
    ux = unique(x);
    counts = zeros(size(ux));
    for i = 1:numel(ux)
        counts(i) = sum(x == ux(i));
    end
    maxc = max(counts);
    candidates = ux(counts == maxc);
    chosen = min(candidates);
end

function countsMap = countLabels(lbls)
    countsMap = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(lbls)
        k = char(lbls(i));
        if isKey(countsMap, k)
            countsMap(k) = countsMap(k) + 1;
        else
            countsMap(k) = 1;
        end
    end
end

function countsByR = countLabelsByRadius(lbls, usedRadius, maxR)
    % countsByR: containers.Map where each key is an area name (char),
    % value is a 1x(maxR+2) vector:
    %   [dist0, dist1, ..., distMaxR, distUnresolved]
    % unresolved is usedRadius == -1 (out-of-bounds or no label found within maxR)
    assert(numel(lbls) == numel(usedRadius), 'lbls and usedRadius must have same length.');

    countsByR = containers.Map('KeyType','char','ValueType','any');

    for i = 1:numel(lbls)
        area = char(lbls(i));
        r = usedRadius(i);

        if ~isKey(countsByR, area)
            countsByR(area) = zeros(1, (maxR+1) + 1); % dist0..distMaxR plus unresolved
        end
        vec = countsByR(area);

        if r >= 0 && r <= maxR
            vec(r+1) = vec(r+1) + 1;
        else
            % unresolved bucket
            vec(end) = vec(end) + 1;
        end

        countsByR(area) = vec;
    end
end

function T = countsByRadiusToTable(countsByR, maxR, areaColName)
    % Output table columns:
    %   Area | Dist0 | Dist1 | ... | DistMaxR | Unresolved | Total
    keysC = string(countsByR.keys);
    nK = numel(keysC);

    distMat = zeros(nK, (maxR+1) + 1); % dist0..distMaxR + Unresolved
    for i = 1:nK
        v = countsByR(char(keysC(i)));
        distMat(i,:) = v(:)';
    end

    total = sum(distMat, 2);

    varNames = strings(1, 1 + (maxR+1) + 1 + 1); % area + dists + unresolved + total
    varNames(1) = areaColName;
    for r = 0:maxR
        varNames(1 + (r+1)) = "Dist" + string(r);
    end
    varNames(1 + (maxR+1) + 1) = "Unresolved";
    varNames(end) = "Total";

    T = table(keysC(:), 'VariableNames', varNames(1));
    for r = 0:maxR
        T.("Dist"+string(r)) = distMat(:, r+1);
    end
    T.Unresolved = distMat(:, end);
    T.Total = total;

    T = sortrows(T, "Total", "descend");
end

function T = radiusSummaryTable(usedRadius, maxR, labelPrefix)
    % Overall counts by radius (how many voxels were exact vs snapped by r=1..maxR vs unresolved)
    % Columns: Type | Radius | Count
    radii = [0:maxR, -1];
    counts = zeros(size(radii));

    for i = 1:numel(radii)
        rr = radii(i);
        counts(i) = sum(usedRadius == rr);
    end

    radiusLabel = strings(numel(radii),1);
    for i = 1:numel(radii)
        if radii(i) == -1
            radiusLabel(i) = "Unresolved";
        else
            radiusLabel(i) = "Dist" + string(radii(i));
        end
    end

    Type = repmat(string(labelPrefix), numel(radii), 1);
    T = table(Type, radiusLabel, counts(:), 'VariableNames', ["Type","Radius","Count"]);
end

function T = countsToTable(countsMap, nameCol, countCol)
    keysC = string(countsMap.keys);
    valsC = cell2mat(countsMap.values);

    T = table(keysC(:), valsC(:), 'VariableNames',[nameCol, countCol]);
    T = sortrows(T, countCol, 'descend');
end

function block = tableToPCBlockCell(T, pcLabel)
    % Convert a table into a cell "block" with an extra first column:
    % Header: [pcLabel, <table variable names...>]
    % Rows:   ["", <table row values...>]
    hdr = string(T.Properties.VariableNames);
    nRows = height(T);
    nCols = width(T);

    block = cell(nRows + 1, nCols + 1);

    % Header row
    block{1,1} = pcLabel;
    for c = 1:nCols
        block{1,c+1} = char(hdr(c));
    end

    % Data rows
    data = table2cell(T);
    for r = 1:nRows
        block{r+1,1} = ""; % blank under "PC #"
        for c = 1:nCols
            block{r+1,c+1} = data{r,c};
        end
    end
end

function combined = combinePCBlocksHorizontally(blocks)
    % Combine blocks side-by-side with ONE blank spacer column between them.
    % Also pads shorter blocks with blanks so all align.
    nB = numel(blocks);
    heights = zeros(1,nB);
    widths  = zeros(1,nB);
    for i = 1:nB
        heights(i) = size(blocks{i},1);
        widths(i)  = size(blocks{i},2);
    end
    H = max(heights);

    spacer = 1; % one blank column between PC blocks

    totalW = sum(widths) + spacer*(nB-1);
    combined = cell(H, totalW);
    combined(:) = {""};

    colStart = 1;
    for i = 1:nB
        b = blocks{i};
        h = size(b,1);
        w = size(b,2);

        combined(1:h, colStart:(colStart+w-1)) = b;

        colStart = colStart + w;
        if i < nB
            colStart = colStart + spacer;
        end
    end
end


%% 2d) BROADNESS VISUALIZATION (ALTERNATIVE SCENARIO WITH OPTIONAL INPUTS)

% This section demonstrates the same function as 2a)  
% but with optional settings provided. Any missing arguments  
% will automatically use their default values.  

%%% ------------------- USER SETTINGS ------------------- %%%

% Minimal user settings: output folder
Options = [];
Options.name_nii = '/main_path/Output'; %output folder
load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;
Options.WhichPlots = [1 1 1 1 1]; %which plots to be generated
Options.ncomps = [1 2 3]; %indices of PCs to be plotted (all plots)
Options.ncomps_var = 60; %number of PCs to be plotted (only in Variance plot)
Options.Labels = {'Global', 'Local'}; %experimental condition labels
Options.color_PCs = [
    0.4,    0.761,  0.647;   % teal-green
    0.988,  0.553,  0.384;   % coral
    0.553,  0.627,  0.796;   % periwinkle
    0.906,  0.541,  0.765;   % pink-purple
    0.651,  0.847,  0.329;   % green
    1.000,  0.851,  0.184;   % yellow
    0.898,  0.769,  0.580;   % beige
    0.702,  0.702,  0.702    % gray
]; %RGB code colors for PCs
Options.color_conds = [
    0.106,  0.620,  0.467;   % green
    0.851,  0.373,  0.008;   % orange
    0.459,  0.439,  0.702;   % purple
    0.906,  0.161,  0.541;   % magenta
    0.400,  0.651,  0.118;   % lime green
    0.902,  0.671,  0.008;   % gold
    0.651,  0.463,  0.114;   % brown
    0.4,    0.4,    0.4      % gray
]; %RGB code colors for experimental conditions

% If you wish to remove the cerebellum voxels (not included in 3D brain template (#4)), please set 'remove_cerebellum_label' to 1
% NOTE: This removal works only for 8mm brain
remove_cerebellum_label = 0;
if remove_cerebellum_label == 1
    load([path_home '/BROADNESS_External/cerebellum_coords.mat']); %only cerebellar voxels
    % Remove cerebellar voxels since they are not included in the 3D brain template (#4)
    [~, idx_cerebellum] = ismember(MNI8, cerebellum_coords, 'rows');  % find cerebellum indexes in MNI coordinates matrix (all voxels)
    MNI8(idx_cerebellum~=0,:) = nan; %assigning nans to MNI coordinates matrix
    Options.MNI_coords = MNI8; %assigning the MNI coordinates of the data for visualization purposes (both 3D main template (#4) and nifti images (#5))
end

%%% ------------------ COMPUTATION --------------------- %%%

% Run each of the 5 lines below in seperate

% Visualize brain networks features
%BROADNESS_Visualizer(BROADNESS_young_local,Options)

%Options.Labels = {'Younger adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_young_local, Options);
%Options.Labels = {'Older adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_older_local, Options); 
%Options.Labels = {'Younger adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_young_global, Options);
%Options.Labels = {'Older adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_older_global, Options);

%%

%% *** ANALYSIS ON SPATIAL-TEMPORAL FEATURES OF BROADNESS BRAIN NETWORKS ***

%%

%% 3) PHASE SPACE EMBEDDING AND RECURRENCE QUANTIFICATION ANALYSIS 

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply use the structure outputted by the BROADNESS_NetworkEstimation function
% Additional optional inputs can be provided, as described in the function. 
% For 3D visualization: 'principalcomps',[1:3] below
% For 2D visualization: 'principalcomps',[1:2] below
% For switching between individual x, y, z limits and limits standard across all plots, edit 224-225 (2D) or line
% 257-259 (3D) in BROADNESS_PhaseSpace_RQA.m: Insert commented lines for
% individual limits, and defined limits for standard limits across all
% plots.

%%% ------------------ COMPUTATION --------------------- %%%


% Run each of these sections in isolation, either 2D or 3D

% 2D

%RQA_BROADNESS = BROADNESS_PhaseSpace_RQA(BROADNESS,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'All data');

%RQA_BROADNESS_young_local = BROADNESS_PhaseSpace_RQA(BROADNESS_young_local,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yl');
%RQA_BROADNESS_older_local = BROADNESS_PhaseSpace_RQA(BROADNESS_older_local,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Ol');
%RQA_BROADNESS_young_global = BROADNESS_PhaseSpace_RQA(BROADNESS_young_global,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yg');
%RQA_BROADNESS_older_global = BROADNESS_PhaseSpace_RQA(BROADNESS_older_global,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Og');

% 3D

%RQA_BROADNESS = BROADNESS_PhaseSpace_RQA(BROADNESS,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'All data');

%RQA_BROADNESS_young_local = BROADNESS_PhaseSpace_RQA(BROADNESS_young_local,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yl');
%RQA_BROADNESS_older_local = BROADNESS_PhaseSpace_RQA(BROADNESS_older_local,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Ol');
%RQA_BROADNESS_young_global = BROADNESS_PhaseSpace_RQA(BROADNESS_young_global,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yg');
%RQA_BROADNESS_older_global = BROADNESS_PhaseSpace_RQA(BROADNESS_older_global,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Og');


%% 4) SPATIAL GRADIENTS EMBEDDING AND CLUSTERING ANALYSIS

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply use the structure outputted by the BROADNESS_NetworkEstimation function
% Additional optional inputs can be provided, as described in the function. 
% Delete 'scatterplots','all' to only visualize optimal k in the cluster plot.

load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;

%%% ------------------ COMPUTATION --------------------- %%%
% All data
SPATIAL_GRADIENTS_BROADNESS = BROADNESS_SpatialGradients(BROADNESS,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);

% Group and condition specific
SPATIAL_GRADIENTS_BROADNESS_young_local = BROADNESS_SpatialGradients(BROADNESS_young_local,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);
SPATIAL_GRADIENTS_BROADNESS_older_local = BROADNESS_SpatialGradients(BROADNESS_older_local,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);
SPATIAL_GRADIENTS_BROADNESS_young_global = BROADNESS_SpatialGradients(BROADNESS_young_global,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);
SPATIAL_GRADIENTS_BROADNESS_older_global = BROADNESS_SpatialGradients(BROADNESS_older_global,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);


%% 5a) STATS - RQA ANOVA TEST

%% --- Minimal 2×2 mixed ANOVAs on RQA metrics ---------------------------
% Needs:
% - RQA_BROADNESS.RQA_metrics : {Nsubj x 1} cell, each a 2x8 table (1=Global, 2=Local)
% - young_subj, older_subj     : subject indices for each group (1..Nsubj)

RQAcell = RQA_BROADNESS.RQA_metrics;
N = numel(RQAcell);

% Group vector per subject
Group = strings(N,1);
Group(young_subj) = "Young";
Group(older_subj) = "Old";
Group = categorical(Group, ["Young","Old"]);

% Metric names (take from first non-empty entry)
idx1 = find(~cellfun(@isempty,RQAcell),1,'first');
metricNames = RQAcell{idx1}.Properties.VariableNames;

% Container for results
Res = table('Size',[numel(metricNames) 4], ...
            'VariableTypes',{'string','double','double','double'}, ...
            'VariableNames',{'Metric','p_Group','p_Condition','p_Interaction'});

fprintf('\n=== Mixed ANOVA (Group × Condition) per metric ===\n');
for m = 1:numel(metricNames)
    met = metricNames{m};

    % Extract Global/Local values for all subjects
    G = nan(N,1); L = nan(N,1);
    for s = 1:N
        T = RQAcell{s};
        if ~isempty(T) && all(ismember(met, T.Properties.VariableNames))
            G(s) = T{1,met};  % row 1 = Global
            L(s) = T{2,met};  % row 2 = Local
        end
    end

    % Build long table: Subject, Group, Condition, Value
    Subject   = repelem((1:N).',2);
    GroupLong = [Group; Group];
    Condition = categorical([repmat("Global",N,1); repmat("Local",N,1)], ["Global","Local"]);
    Value     = [G; L];

    % Drop missing rows (robust to NaNs)
    keep = isfinite(Value) & (GroupLong~="");
    Tlong = table(Subject(keep), GroupLong(keep), Condition(keep), Value(keep), ...
                  'VariableNames',{'Subject','Group','Condition','Value'});

    % Mixed effects model: fixed Group, Condition, Interaction effects. Random effects: intercept per subject
    lme = fitlme(Tlong, 'Value ~ Group*Condition + (1|Subject)');

    % F-tests for fixed effects
    % Get F-tests with Satterthwaite df:
    % Group main effect: Young vs Old (averaging over Global/Local): Do young and Old differ in the metric on average?
    % Condition main effect: Global vs Local (within-subject, averaged over groups): Is Global different from Local within subjects?
    % Interaction: the Global–Local difference differs between Young and Old: Interaction: Is the Global–Local contrast different between groups?
    
    A = anova(lme,'DFMethod','Satterthwaite');
    pG  = A.pValue(strcmp(A.Term,'Group'));
    pC  = A.pValue(strcmp(A.Term,'Condition'));
    pIx = A.pValue(strcmp(A.Term,'Group:Condition'));

   
    Res.Metric(m)       = string(met);
    Res.p_Group(m)      = pG;
    Res.p_Condition(m)  = pC;
    Res.p_Interaction(m)= pIx;

    fprintf('%-8s  Group p=%.4g | Condition p=%.4g | Interaction p=%.4g\n', ...
            met, pG, pC, pIx);
end

disp(Res);

%% 5b) Plotting RQA metrics — ONE FIGURE (2x4 grid), VIOLINS per metric with mean, ±SD
% Uses subject level data computed by RQA_BROADNESS.RQA_metrics -> table with rows: 1=Global, 2=Local)
% Requires variables young and older

plot_outdir = '/main_path/Output/RQA';

% ---- Check RQA structure ----
if ~isfield(RQA_BROADNESS, 'RQA_metrics')
    error('RQA_BROADNESS.RQA_metrics not found.');
end
RQAcell = RQA_BROADNESS.RQA_metrics; % {nPart} cell, each a table with rows: 1=Global, 2=Local

% ---- Metric names ----
metrics = {'RR','L','TT','V_max','ENTR','DET','LAM','DIV'};
nice    = {'Recurrence Rate (RR)','Mean Diagonal Length (L)','Trapping Time (TT)','Max Vertical (V_{max})', 'Entropy (ENTR)','Determinism (DET)','Laminarity (LAM)','Divergence (DIV)'};

% ---- 4 combos: Young-Local, Old-Local, Young-Global, Old-Global ----
COND_GLOBAL = 1; COND_LOCAL = 2;
combos = { ...
    struct('name','Younger adults - Local effect',  'condIdx',COND_LOCAL,  'ids',young_subj), ...
    struct('name','Older adults - Local effect',    'condIdx',COND_LOCAL,  'ids',older_subj), ...
    struct('name','Younger adults - Global effect', 'condIdx',COND_GLOBAL, 'ids',young_subj), ...
    struct('name','Older adults - Global effect',   'condIdx',COND_GLOBAL, 'ids',older_subj)  ...
};

% Colors per combo (consistent across plots)
comb_colors = [ ...
    0.1 0.2 0.5;  % Young-Local   (blue)
    0.5 0.1 0.1;  % Old-Local     (red)
    0.2 0.4 0.8;  % Young-Global  (light blue)
    0.8 0.2 0.2]; % Old-Global    (light red)

% Helper to extract a vector from RQAcell for given metric, ids, condition index
get_data = @(metric, ids, condIdx) local_collect_vals(RQAcell, metric, ids, condIdx);

% Seed for reproducible jitter
rng(42);

% ---- One big figure with 2x4 layout (last tile = legend) ----
H = figure('Color','w','Units','pixels','Position',[30 30 2600 1700]);
tl = tiledlayout(H, 2, 4, 'Padding','compact', 'TileSpacing','compact');

for m = 1:numel(metrics)
    met   = metrics{m};
    label = nice{m};

    % Gather data for the four combos
    D = cell(1,4);
    for k = 1:4
        D{k} = get_data(met, combos{k}.ids, combos{k}.condIdx);
        D{k} = D{k}(isfinite(D{k}));
    end

    ax = nexttile(tl, m); hold(ax,'on');

    % X positions for the 4 violins in this subplot
    xpositions = 1:4;

    % Draw each combo (mean + SD only; no SEM)
    for k = 1:4
        draw_violin_mean_SD(ax, D{k}, xpositions(k), comb_colors(k,:));
    end

    % --- Y-limits from rendered graphics (captures kernel tails & markers) ---
    kids = findall(ax, 'Type','patch','-or','Type','line','-or','Type','scatter','-or','Type','area');
    Y = [];
    for h = reshape(kids,1,[])
        if isprop(h,'YData')
            y = get(h,'YData');
            if ~isempty(y), Y = [Y; y(:)]; end %#ok<AGROW>
        end
    end
    finiteY = Y(isfinite(Y));
    if isempty(finiteY), finiteY = 0; end
    ymin = min(finiteY); ymax = max(finiteY);
    if ~isfinite(ymin), ymin = 0; end
    if ~isfinite(ymax), ymax = 1; end

    if ymax == ymin
        pad = max(0.1*max(1,abs(ymax)), 0.1);
    else
        % Larger pad than before to prevent clipping of violin caps
        pad = 0.12 * (ymax - ymin);
    end

    ylim(ax, [ymin - pad, ymax + pad]);


    xlim(ax, [0.5 4.5]);

    % Remove ALL x-axis labels/ticklabels as in style
    set(ax,'XTick', xpositions, 'XTickLabel', []);
    xlabel(ax, '');
    ylabel(ax, label);
    title(ax, label, 'FontWeight');
    grid(ax,'on'); box(ax,'on');
end

% ---- Legend tile in position 9 ----
axL = nexttile(tl, 9); cla(axL); hold(axL,'on'); axis(axL,'off');

legend_items = gobjects(0);
legend_labels = strings(0);

for k = 1:4
    legend_items(end+1) = plot(axL, NaN, NaN, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', comb_colors(k,:), 'MarkerEdgeColor','k'); %#ok<AGROW>
    legend_labels(end+1) = combos{k}.name; %#ok<AGROW>
end

lg = legend(axL, legend_items, legend_labels, 'Location','northwest');
set(lg,'Interpreter','none','FontSize',13,'Box','off');

% Overall title
sgtitle(tl, 'Recurrence Quantification Analysis', 'FontWeight','bold');

% Helvetica
set(findall(H, '-property', 'FontName'), 'FontName', 'Helvetica');

% Save outputs
out_png = fullfile(plot_outduntitledir, 'RQA_Violin_AllMetrics_2x4.png');
print(H, out_png, '-dpng', '-r300');
set(H, 'PaperPositionMode','auto');
close(H);
fprintf('Saved 3x3 figure:\n  %s\n  %s\n', out_png);


%% 6a) STATS — Timepoint-wise mixed ANOVA via CLUSTER-BASED permutation (over time)
% PURPOSE
% We test whether the time courses of a single principal component (PC) differ:
%   - between age groups (Young vs Old)  --> "Group main effect"
%   - between conditions (Global vs Local) within subjects --> "Condition main effect"
%   - in their difference across groups --> "Interaction"
%
% We do this separately for PC = 1, 2, 3 in the 0–0.8 s time window.
% Significance is controlled with a cluster-based permutation test:
%   1) Form clusters of adjacent time samples whose effect size exceeds a
%      "cluster-forming threshold" (CFT).
%   2) For each permutation, compute the maximum cluster "mass" (sum of |effect| in cluster).
%   3) An observed cluster is significant if its mass exceeds the (1 - alpha) quantile
%      of that max-mass null distribution (FWER control).
%
% Inputs expected in the workspace (each is a struct with a 4D array):
%   BROADNESS_young_local .TimeSeries_BrainNetworks : [time x PC x cond x subj]
%   BROADNESS_young_global.TimeSeries_BrainNetworks : [time x PC x cond x subj]
%   BROADNESS_older_local .TimeSeries_BrainNetworks : [time x PC x cond x subj]
%   BROADNESS_older_global.TimeSeries_BrainNetworks : [time x PC x cond x subj]
% NOTE: We only use cond index 1 (time series for this PC), so effectively [time x PC x 1 x subj].
%
% Also expected:
%   time : vector of time stamps in seconds (must cover 0–0.8 s)
%
% Outputs (per PC):
%   A .mat file with RES struct containing observed effects, thresholds, significant masks,
%   cluster lists, and convenience metadata; plus several diagnostic PNG plots.

stats_outdir = '/main_path/Output/ANOVA_Stats';
if ~exist(stats_outdir,'dir'); mkdir(stats_outdir); end  % ensure output folder exists

% ---- Parameters ----
alpha         = 0.05;   % family-wise error rate target for clusters (final significance level)
cluster_alpha = 0.05;   % per-sample cluster-forming threshold (liberal; used only to build clusters)
n_perm        = 5000;   % number of permutations for null estimation (increase = more stable, slower)
pcs_to_use    = 1:3;    % PCs to analyze in this script

% ---- Determine dims: we expect [time x PC x cond x subj] ----
sz   = size(BROADNESS_young_local.TimeSeries_BrainNetworks);
nT_bn = sz(1);          % number of time samples available in the data (from "broadness" arrays)
nPC   = sz(2);          % number of principal components available

% ---- align external 'time' vector to nT_bn ----
% We need a time vector aligned to the first dimension of the data arrays.
if ~exist('time','var') || isempty(time)
    error('Variable "time" not found. Define the time vector in seconds.');
end
if numel(time) < nT_bn
    error('Provided time vector (%d) is shorter than BROADNESS time dimension (%d).', numel(time), nT_bn);
end
time_use = time(1:nT_bn);  % trim any extra time points just in case

% ---- analysis window ----
% Restrict the analysis to [0, 0.8] seconds, as specified in the design.
t_mask = (time_use >= 0) & (time_use <= 0.8);
if ~any(t_mask)
    error('No timepoints within [0, 0.8] s after alignment.');
end
t_idx  = find(t_mask);       % indices of the selected time points
t_win  = time_use(t_idx);    % time stamps used in the stats
Twin   = numel(t_idx);       % number of samples in the window
dt_win = median(diff(t_win));% time step (used for info and plotting titles)

% ---- PC availability check ----
if max(pcs_to_use) > nPC
    error('Requested PC index exceeds available PCs (%d).', nPC);
end

% ===================== MAIN LOOP OVER PCs =====================
for pc = pcs_to_use
    % ---- Extract subject x time matrices for this PC ----
    % We want matrices with rows = subjects, columns = time.
    % squeeze(... )' transposes [time x subj] -> [subj x time].
    % Note: we use cond index 1 (the only condition dim used in these arrays).
    YL = squeeze(BROADNESS_young_local. TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_young x T]
    YG = squeeze(BROADNESS_young_global.TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_young x T]
    OL = squeeze(BROADNESS_older_local. TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_old   x T]
    OG = squeeze(BROADNESS_older_global.TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_old   x T]

    % ---- POLARITY FIX: flip Older-Global for PC1 (so stats + masks follow) ----
    if pc == 1
       OG = -OG;
    end

    % ---- Defensive shape checks (fail early if misaligned) ----
    if any([size(YL,2), size(YG,2), size(OL,2), size(OG,2)] ~= Twin)
        error('Time dimension mismatch after slicing for PC %d.', pc);
    end
    if size(YL,1) ~= size(YG,1)
        error('Young group Local/Global subject count mismatch for PC %d.', pc);
    end
    if size(OL,1) ~= size(OG,1)
        error('Old group Local/Global subject count mismatch for PC %d.', pc);
    end

    % ---- Compute OBSERVED effects (each is a 1 x T vector) ----
    % Group main (between-subject): average Local+Global within each subject, then Young minus Old.
    Y_mean = (YL + YG) ./ 2;  % [N_young x T]
    O_mean = (OL + OG) ./ 2;  % [N_old   x T]
    obs_G  = mean(Y_mean, 1, 'omitnan') - mean(O_mean, 1, 'omitnan');  % 1 x T

    % Condition main (within-subject): (Global - Local) pooled across both groups.
    D_y   = YG - YL;                  % [N_young x T] within-subject difference in Young
    D_o   = OG - OL;                  % [N_old   x T] within-subject difference in Old
    D_all = [D_y; D_o];               % [N_total x T] stack both groups
    obs_C = mean(D_all, 1, 'omitnan');% 1 x T pooled within-subject effect

    % Interaction: "difference of differences" = (Young (G-L)) - (Old (G-L)).
    obs_I = mean(D_y, 1, 'omitnan') - mean(D_o, 1, 'omitnan');  % 1 x T

    % ---- Cluster-based permutation tests (each returns a Boolean mask over time, thresholds, and clusters) ----
    % Between-subject permutation (reshuffle group labels) for Group main and Interaction:
    [sig_G, thr_G, clu_G] = perm_between_time(Y_mean, O_mean, n_perm, cluster_alpha, alpha); % Group main
    [sig_I, thr_I, clu_I] = perm_between_time(D_y,   D_o,   n_perm, cluster_alpha, alpha);   % Interaction

    % Within-subject permutation (sign flips) for Condition main:
    [sig_C, thr_C, clu_C] = perm_within_time(D_all, n_perm, cluster_alpha, alpha);           % Condition main

    % ---- Convert significant masks into readable time windows [start, end] in seconds ----
    % These are contiguous segments of significant samples.
    win_G = mask_to_windows(sig_G, t_win);
    win_C = mask_to_windows(sig_C, t_win);
    win_I = mask_to_windows(sig_I, t_win);

    % ---- Save structured results for this PC ----
    RES = struct();
    RES.PC            = pc;                   % which PC we analyzed
    RES.time          = t_win;                % time support used in the test
    RES.dt            = dt_win;               % median sample spacing
    RES.alpha         = alpha;                % final FWER alpha
    RES.cluster_alpha = cluster_alpha;        % cluster-forming per-sample alpha
    RES.n_perm        = n_perm;               % number of permutations
    RES.obs           = struct('Group', obs_G, 'Condition', obs_C, 'Interaction', obs_I);  % observed effect curves
    RES.thresholds    = struct('Group', thr_G, 'Condition', thr_C, 'Interaction', thr_I);  % CFT and critical mass
    RES.sig           = struct('Group', sig_G, 'Condition', sig_C, 'Interaction', sig_I);  % significant masks
    RES.windows       = struct('Group', win_G, 'Condition', win_C, 'Interaction', win_I);  % [start,end] windows
    RES.N             = struct('Young', size(YL,1), 'Old', size(OL,1));                    % sample sizes
    RES.clusters      = struct('Group', clu_G, 'Condition', clu_C, 'Interaction', clu_I);  % raw cluster lists

    save(fullfile(stats_outdir, sprintf('PC%d_Timewise_ANOVA_ClusterPerm.mat', pc)), 'RES', '-v7.3');

    % ---- Quick diagnostic plots (saved per PC) ----
    % Each plot draws the observed curve and shades significant clusters to help visual QC.

    try
        make_effect_plot(t_win, obs_G, sig_G, thr_G, sprintf('PC%d — Group main (Young-Old)', pc), ...
            fullfile(stats_outdir, sprintf('PC%d_GroupMain.png', pc)));
        make_effect_plot(t_win, obs_C, sig_C, thr_C, sprintf('PC%d — Condition main (Global-Local)', pc), ...
            fullfile(stats_outdir, sprintf('PC%d_ConditionMain.png', pc)));
        make_effect_plot(t_win, obs_I, sig_I, thr_I, sprintf('PC%d — Interaction (YoungΔ-OldΔ)', pc), ...
            fullfile(stats_outdir, sprintf('PC%d_Interaction.png', pc)));
    catch ME
        warning('Plotting failed for PC %d: %s', pc, ME.message);
    end

    % ===================== INSERTED PLOTS OF ALL RELEVANT COMPUTED ITEMS =====================
    % 1) Raw condition means (YL, YG, OL, OG), with ANY-effect significant windows shaded.

    try
        fig1 = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1100 450]);
        ax1 = axes(fig1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

        % Group-level mean time courses per condition
        muYL = mean(YL,1,'omitnan'); muYG = mean(YG,1,'omitnan');
        muOL = mean(OL,1,'omitnan'); muOG = mean(OG,1,'omitnan');

        % Colored lines for each condition/group combo
        plot(ax1, t_win, muYL, 'LineWidth',1.2, 'Color',[0.10 0.20 0.60]); % Young-Local
        plot(ax1, t_win, muYG, 'LineWidth',1.2, 'Color',[0.35 0.70 0.95]); % Young-Global
        plot(ax1, t_win, muOL, 'LineWidth',1.2, 'Color',[0.55 0.10 0.70]); % Old-Local
        plot(ax1, t_win, muOG, 'LineWidth',1.2, 'Color',[0.95 0.45 0.70]); % Old-Global

        yl1 = ylim(ax1);

        % Gray shading where ANY of the three effects is significant (union of masks)
        sig_any = sig_G | sig_C | sig_I;
        shade_clusters_union(ax1, t_win, sig_any, yl1, [0.7 0.7 0.7], 0.12);

        title(ax1, sprintf('PC%d — Raw condition means (YL,YG,OL,OG) with significant windows (gray)', pc), 'Interpreter','none');
        xlabel(ax1,'Time (s)'); ylabel(ax1,'Activation (a.u.)');
        legend(ax1, {'Younger-Local','Younger-Global','Older-Local','Older-Global'}, 'Location','northwest');

        print(fig1, fullfile(stats_outdir, sprintf('PC%d_RawConditionMeans.png', pc)), '-dpng','-r200');
        close(fig1);
    catch ME
        warning('Raw means plotting failed for PC %d: %s', pc, ME.message);
    end

    % 2) Group means (mean of Local+Global) with Group-effect clusters shaded.
    try
        fig2 = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1100 450]);
        ax2 = axes(fig2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');

        plot(ax2, t_win, mean(Y_mean,1,'omitnan'), 'LineWidth',1.5, 'Color',[0.2 0.2 0.8]); % Young mean
        plot(ax2, t_win, mean(O_mean,1,'omitnan'), 'LineWidth',1.5, 'Color',[0.8 0.2 0.2]); % Old mean

        yl2 = ylim(ax2);
        shade_clusters_union(ax2, t_win, sig_G, yl2, [0.7 0.7 0.7], 0.15);

        title(ax2, sprintf('PC%d — Group means (Young vs Old) with Group-effect clusters (gray)', pc), 'Interpreter','none');
        xlabel(ax2,'Time (s)'); ylabel(ax2,'Activation (a.u.)');
        legend(ax2, {'Young (mean of Local+Global)','Old (mean of Local+Global)'}, 'Location','northwest');

        print(fig2, fullfile(stats_outdir, sprintf('PC%d_GroupMeans_ClusterGray.png', pc)), '-dpng','-r200');
        close(fig2);
    catch ME
        warning('Group means plotting failed for PC %d: %s', pc, ME.message);
    end

    % 3) Condition differences (G-L) per group and pooled, with appropriate shading for Interaction/Condition.
    try
        fig3 = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1100 500]);
        tlo3 = tiledlayout(fig3, 2,1, 'Padding','compact','TileSpacing','compact');

        % (a) G-L per group with Interaction shading
        ax3a = nexttile(tlo3,1); hold(ax3a,'on'); grid(ax3a,'on'); box(ax3a,'on');
        dY = mean(D_y,1,'omitnan');  % Young (G-L)
        dO = mean(D_o,1,'omitnan');  % Old   (G-L)
        plot(ax3a, t_win, dY, 'LineWidth',1.4, 'Color',[0.2 0.6 0.9]);
        plot(ax3a, t_win, dO, 'LineWidth',1.4, 'Color',[0.9 0.4 0.2]);
        yl3a = ylim(ax3a);
        shade_clusters_union(ax3a, t_win, sig_I, yl3a, [0.7 0.7 0.7], 0.15);
        title(ax3a, sprintf('PC%d — (G-L) per group with Interaction clusters (gray)', pc), 'Interpreter','none');
        xlabel(ax3a,'Time (s)'); ylabel(ax3a,'Δ (G-L)');
        legend(ax3a, {'Young: G-L','Old: G-L'}, 'Location','northwest');

        % (b) Pooled Condition effect with Condition shading
        ax3b = nexttile(tlo3,2); hold(ax3b,'on'); grid(ax3b,'on'); box(ax3b,'on');
        plot(ax3b, t_win, obs_C, 'LineWidth',1.6, 'Color',[0.0 0.45 0.74]);
        yl3b = ylim(ax3b);
        shade_clusters_union(ax3b, t_win, sig_C, yl3b, [0.7 0.7 0.7], 0.15);
        title(ax3b, sprintf('PC%d — Condition main effect (pooled G-L) with clusters (gray)', pc), 'Interpreter','none');
        xlabel(ax3b,'Time (s)'); ylabel(ax3b,'Δ (G-L) pooled');

        print(fig3, fullfile(stats_outdir, sprintf('PC%d_ConditionAndInteractionDetails.png', pc)), '-dpng','-r200');
        close(fig3);
    catch ME
        warning('Condition/Interaction plotting failed for PC %d: %s', pc, ME.message);
    end

    % 4) Binary significance masks (three rows: Group, Condition, Interaction) to see where each effect is on.
    try
        fig4 = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1100 400]);
        ax4 = axes(fig4); hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');

        y0 = 0;
        plot_sig_mask(ax4, t_win, sig_G, y0+2, 'Group');
        plot_sig_mask(ax4, t_win, sig_C, y0+1, 'Condition');
        plot_sig_mask(ax4, t_win, sig_I, y0+0, 'Interaction');

        yticks(ax4, y0+[0 1 2]); yticklabels(ax4, {'Interaction','Condition','Group'});
        xlabel(ax4,'Time (s)'); ylim(ax4,[-0.5 2.5]); xlim(ax4,[t_win(1) t_win(end)]);
        title(ax4, sprintf('PC%d — Binary significance masks (cluster-corrected)', pc), 'Interpreter','none');

        print(fig4, fullfile(stats_outdir, sprintf('PC%d_SignificanceMasks.png', pc)), '-dpng','-r200');
        close(fig4);
    catch ME
        warning('Significance mask plotting failed for PC %d: %s', pc, ME.message);
    end
    % ===================== END INSERTED PLOTS =====================

end

fprintf('6 updated: cluster-based stats saved in %s\n', stats_outdir);


% ===================== LOCAL HELPERS =====================
% The helper functions below implement the permutation logic, cluster detection,
% mask->window conversion, robust quantiles, and plotting utilities.

function [sig_mask, thr, clusters] = perm_between_time(A, B, nperm, c_alpha, fwer_alpha)
% Between-subject cluster permutation across time on mean differences.
% INPUTS:
%   A : [N_A x T]   group A subjects x time
%   B : [N_B x T]   group B subjects x time
%   nperm : number of permutations
%   c_alpha : pointwise alpha to define clusters (CFT)
%   fwer_alpha : final cluster-level alpha (FWER)
% OUTPUTS:
%   sig_mask : 1 x T logical, true where any significant cluster spans
%   thr      : struct with fields 'cluster_forming' (CFT) and 'cluster_mass' (critical mass)
%   clusters : struct with cluster indices, masses, and which ones are significant (sig_flags)
    Na = size(A,1); Nb = size(B,1); T = size(A,2);
    if size(B,2) ~= T; error('Time dim mismatch A/B.'); end

    % Observed difference (Young-Old or D_y - D_o depending on call site)
    obs = mean(A,1,'omitnan') - mean(B,1,'omitnan'); % 1 x T

    % Concatenate to reshuffle labels in permutations
    X = [A; B];

    % --- Build the cluster-forming threshold (CFT) ---
    % Strategy: pool the absolute per-sample effects from permuted splits.
    % Then take the (1 - c_alpha) quantile as CFT.
    pooled_abs = nan(nperm*T,1);
    for p = 1:nperm
        idx = randperm(Na+Nb);
        Ap  = X(idx(1:Na),       :);
        Bp  = X(idx(Na+1:end),   :);
        dpp = mean(Ap,1,'omitnan') - mean(Bp,1,'omitnan');
        pooled_abs((p-1)*T+1 : p*T) = abs(dpp(:));
    end
    cft = quantile_fast(pooled_abs, 1 - c_alpha);

    % --- Build the null of MAX cluster mass (FWER control) ---
    % For each permutation, compute the largest cluster mass.
    max_masses = nan(nperm,1);
    for p = 1:nperm
        idx = randperm(Na+Nb);
        Ap  = X(idx(1:Na),       :);
        Bp  = X(idx(Na+1:end),   :);
        dpp = mean(Ap,1,'omitnan') - mean(Bp,1,'omitnan');
        max_masses(p) = max_cluster_mass(abs(dpp), cft);
    end
    crit_mass = quantile_fast(max_masses, 1 - fwer_alpha);

    % --- Observed clusters and significance decision ---
    [clu_idx, clu_mass] = find_clusters(abs(obs), cft);
    sig_flags = clu_mass > crit_mass;

    % Convert significant clusters into a 1 x T boolean mask
    sig_mask = false(1,T);
    for k = 1:numel(clu_idx)
        if sig_flags(k), sig_mask(clu_idx{k}) = true; end
    end

    thr      = struct('cluster_forming', cft, 'cluster_mass', crit_mass);
    clusters = struct('idx_list', {clu_idx}, 'masses', clu_mass, 'sig_flags', sig_flags);
end

function [sig_mask, thr, clusters] = perm_within_time(D, nperm, c_alpha, fwer_alpha)
% Within-subject cluster permutation by random sign flips (paired design).
% INPUT:
%   D : [N_subj x T] per-subject within-contrast (e.g., Global - Local)
% OUTPUT:
%   same as perm_between_time
    [Ns, T] = size(D);
    obs = mean(D,1,'omitnan');  % observed pooled within-subject effect

    % Pre-generate ±1 flip matrix (Ns x nperm) for speed.
    flips = (randi(2, Ns, nperm)*2 - 3);  % entries are +1 or -1

    % Build CFT via pooled absolute permuted effects
    pooled_abs = nan(nperm*T,1);
    for p = 1;nperm
        Dp  = D .* flips(:,p);            % random sign flip per subject
        dpp = mean(Dp,1,'omitnan');       % permuted effect
        pooled_abs((p-1)*T+1 : p*T) = abs(dpp(:));
    end
    cft = quantile_fast(pooled_abs, 1 - c_alpha);

    % Build max cluster mass null (FWER)
    max_masses = nan(nperm,1);
    for p = 1:nperm
        Dp  = D .* flips(:,p);
        dpp = mean(Dp,1,'omitnan');
        max_masses(p) = max_cluster_mass(abs(dpp), cft);
    end
    crit_mass = quantile_fast(max_masses, 1 - fwer_alpha);

    % Observed clusters
    [clu_idx, clu_mass] = find_clusters(abs(obs), cft);
    sig_flags = clu_mass > crit_mass;

    % Convert to mask
    sig_mask = false(1,T);
    for k = 1:numel(clu_idx)
        if sig_flags(k), sig_mask(clu_idx{k}) = true; end
    end

    thr      = struct('cluster_forming', cft, 'cluster_mass', crit_mass);
    clusters = struct('idx_list', {clu_idx}, 'masses', clu_mass, 'sig_flags', sig_flags);
end

function m = max_cluster_mass(abs_series, cft)
% Return the largest cluster mass in a 1 x T absolute effect series.
% If no clusters form, return 0.
    [~, masses] = find_clusters(abs_series, cft);
    if isempty(masses), m = 0; else, m = max(masses); end
end

function [idx_list, masses] = find_clusters(abs_series, cft)
% Find contiguous indices where abs_series > cft and compute their masses.
% INPUTS:
%   abs_series : 1 x T nonnegative vector (absolute effect)
%   cft        : scalar cluster-forming threshold
% OUTPUTS:
%   idx_list   : cell array; each cell is a vector of indices for one cluster
%   masses     : vector; sum(abs_series) within each cluster
    mask = abs_series(:)' > cft; % suprathreshold boolean vector
    idx  = find(mask);
    idx_list = {}; masses = [];
    if isempty(idx), return; end

    % Break points where clusters separate (difference > 1)
    br = [1, find(diff(idx) > 1)+1, numel(idx)+1];
    for b = 1:numel(br)-1
        seg = idx(br(b):br(b+1)-1);
        idx_list{end+1} = seg;             %#ok<AGROW>
        masses(end+1)   = sum(abs_series(seg)); %#ok<AGROW>
    end
end

function W = mask_to_windows(mask, tvec)
% Convert a 1 x T boolean mask into [nWin x 2] time windows using tvec.
% Each row: [start_time, end_time] in seconds for a contiguous "true" segment.
    mask = mask(:)'; on = find(mask);
    W = [];
    if isempty(on), return; end
    br = [1, find(diff(on) > 1)+1, numel(on)+1];
    for b = 1:numel(br)-1
        seg = on(br(b):br(b+1)-1);
        W(end+1, :) = [tvec(seg(1)), tvec(seg(end))]; %#ok<AGROW>
    end
end

function q = quantile_fast(x, p)
% Robust quantile that ignores NaNs using prctile.
% INPUTS:
%   x : vector of samples
%   p : desired quantile in [0,1]
% OUTPUT:
%   q : p-quantile of x
    x = x(isfinite(x));
    if isempty(x), q = NaN; return; end
    q = prctile(x, p*100);
end

function make_effect_plot(t, obs, sig, thr, ttl, outpng)
% Simple 1D effect plot:
%   - line = observed effect over time
%   - shaded areas = significant clusters (FWER controlled)
%   - title shows CFT and critical mass for transparency
    h = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1100 900]);
    plot(t, obs, 'LineWidth', 1.5); hold on; grid on; box on;
    yl = ylim;
    shade_clusters(t, sig, yl);
    title(sprintf('%s | CFT=%.3g | crit mass=%.3g', ttl, thr.cluster_forming, thr.cluster_mass), 'Interpreter','none');
    xlabel('Time (s)'); ylabel('\Delta');
    print(h, outpng, '-dpng', '-r200');
    close(h);
end

function shade_clusters(t, sig, yl)
% Fill rectangles over contiguous significant segments using y-limits "yl".
    on = find(sig(:)');
    if isempty(on), return; end
    br = [1, find(diff(on) > 1)+1, numel(on)+1];
    for b = 1:numel(br)-1
        seg = on(br(b):br(b+1)-1);
        x1 = t(seg(1)); x2 = t(seg(end));
        area([x1 x2], [yl(2) yl(2)], yl(1), 'FaceAlpha', 0.12, 'EdgeColor','none');
    end
end

% ===================== INSERTED LOCAL HELPERS FOR PLOTTING =====================
function shade_clusters_union(ax, t, sig_mask, yl, face_rgb, face_alpha)
% Shade regions where sig_mask==true on provided axes "ax".
% Used to show "any effect" or effect-specific masks.
    on = find(sig_mask(:)');
    if isempty(on), return; end
    br = [1, find(diff(on) > 1)+1, numel(on)+1];
    for b = 1:numel(br)-1
        seg = on(br(b):br(b+1)-1);
        x1 = t(seg(1)); x2 = t(seg(end));
        patch('XData',[x1 x2 x2 x1], 'YData',[yl(1) yl(1) yl(2) yl(2)], ...
              'FaceColor',face_rgb, 'FaceAlpha',face_alpha, 'EdgeColor','none', ...
              'Parent',ax, 'HitTest','off');
    end
end

function plot_sig_mask(ax, t, mask, ybase, labeltxt)
% Draw a thick line for each contiguous "true" region at vertical level ybase.
% Useful for a compact "barcode" visualization of significance.
    on = find(mask(:)');
    if isempty(on)
        % draw a faint reference line if there is no significance
        plot(ax, t([1 end]), [ybase ybase], ':', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
        return;
    end
    br = [1, find(diff(on) > 1)+1, numel(on)+1];
    for b = 1:numel(br)-1
        seg = on(br(b):br(b+1)-1);
        x1 = t(seg(1)); x2 = t(seg(end));
        plot(ax, [x1 x2], [ybase ybase], '-', 'LineWidth', 3);
    end
    % label the row in the margin
    text(mean(t), ybase+0.15, labeltxt, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize',9, 'Color',[0 0 0], 'Interpreter','none');
end

%% 6b) TIMESERIES 2×3 PC-WISE SUMMARY PLOT using ANOVA results from 6a)
% Layout:
%   [ PC1 | PC2 | PC3
%   | LEGEND ]

% Uses:
%   - BROADNESS_young_local/global, BROADNESS_older_local/global
%   - stats_outdir and RES files created in 6)
%   - 'time' vector in workspace

% Output:
%   - PNG saved in stats_outdir

% --------- sanity on 6) outputs ---------
need_files = {
    fullfile(stats_outdir,'PC1_Timewise_ANOVA_ClusterPerm.mat')
    fullfile(stats_outdir,'PC2_Timewise_ANOVA_ClusterPerm.mat')
    fullfile(stats_outdir,'PC3_Timewise_ANOVA_ClusterPerm.mat')
};
have_files = cellfun(@(p) exist(p,'file')==2, need_files);
if ~all(have_files)
    warning('7: Missing some 6 result files. Found=%s', mat2str(have_files));
end

% --------- require 'time' ---------
if ~exist('time','var') || isempty(time)
    error('7: time vector not found in workspace.');
end

% --------- figure aesthetics ---------
tags   = {'Young_Local','Old_Local','Young_Global','Old_Global'};
labels = {'Younger - Local','Older - Local','Younger - Global','Older - Global'};

comb_colors = [ ...
    0.1 0.2 0.7;  % Young-Local   (blue)
    0.5 0.1 0.1;  % Old-Local     (red)
    0.4 0.6 0.85;  % Young-Global  (light blue)
    0.8 0.2 0.2]; % Old-Global    (light red)

% Gray shaded areas for any significant effect
shade_gray  = [0.7 0.7 0.7];
shade_alpha = 0.15;
ribbon_alpha = 0.20;   % transparency for SE ribbons

% Bottom indicator colors per effect
clr_group = [0 0.8 0.2];  % green
clr_cond  = [0.6 0 0.8]; % purple
clr_int   = [1 0.6 0];   % orange

% Y-limits (adjust as required)
ymin = -1200; ymax = 1000;

% --------- requested full x-window + dotted vertical reference lines ---------
xlim_target = [-0.1 0.8];
vlines = [0.0];

% Indices and time vector for plotting (full requested window)
t_plot_idx = find(time >= xlim_target(1) & time <= xlim_target(2));
t_plot     = time(t_plot_idx);
if isempty(t_plot)
    error('7: No samples within [%g, %g] in the time vector.', xlim_target(1), xlim_target(2));
end
dt_plot = median(diff(t_plot));

% Pack the four BROADNESS structs in the same order as labels
BMAP = containers.Map( ...
    tags, ...
    {BROADNESS_young_local, BROADNESS_older_local, BROADNESS_young_global, BROADNESS_older_global} ...
);

% --------- precompute mu/se for each combo and each PC lazily (on FULL window) ---------
get_mu_se = @(tag, pc) extract_mu_se_proper_7(BMAP(tag), t_plot_idx, pc);

% --------- plotting ---------
pc_list = [1 2 3];
h  = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1400 1000]);

set(h, 'DefaultAxesFontName', 'Helvetica', ...
       'DefaultTextFontName', 'Helvetica', ...
       'DefaultLegendFontName', 'Helvetica', ...
       'DefaultAxesFontSize', 8, ...
       'DefaultTextFontSize', 8, ...
       'DefaultLegendFontSize', 8, ...
       'DefaultAxesTitleFontSizeMultiplier', 1, ...
       'DefaultAxesLabelFontSizeMultiplier', 1);
tlo = tiledlayout(2,3, 'Padding','compact','TileSpacing','compact');

for ip = 1:numel(pc_list)
    pc = pc_list(ip);

    % Load stats (for shading/indicators only)
    matfile = fullfile(stats_outdir, sprintf('PC%d_Timewise_ANOVA_ClusterPerm.mat', pc));
    if exist(matfile,'file')~=2
        warning('7: Missing %s. Skipping PC%d.', matfile, pc);
        nexttile; axis off; title(sprintf('PC%d (missing results)', pc));
        continue;
    end
    S  = load(matfile, 'RES');
    Rt = S.RES;                       % has RES.time, RES.sig.*

    ax = nexttile(tlo, ip); hold(ax,'on');

    % === A) Build grid-aligned masks on t_plot (guarantees perfect alignment) ===
    maskG_plot = sig_to_plot_mask(Rt.sig.Group,      Rt.time, t_plot);
    maskC_plot = sig_to_plot_mask(Rt.sig.Condition,  Rt.time, t_plot);
    maskI_plot = sig_to_plot_mask(Rt.sig.Interaction,Rt.time, t_plot);

    % Convert masks to windows on t_plot, expanding half-bin to cover last sample
    wins_G = mask_to_windows_on_grid(maskG_plot, t_plot, dt_plot, xlim_target);
    wins_C = mask_to_windows_on_grid(maskC_plot, t_plot, dt_plot, xlim_target);
    wins_I = mask_to_windows_on_grid(maskI_plot, t_plot, dt_plot, xlim_target);

    % === B) Draw shading (behind curves) using the grid-aligned windows ===
    draw_gray_windows_7(ax, wins_G, ymin, ymax, shade_gray, shade_alpha);
    draw_gray_windows_7(ax, wins_C, ymin, ymax, shade_gray, shade_alpha);
    draw_gray_windows_7(ax, wins_I, ymin, ymax, shade_gray, shade_alpha);

    % === C) Mean + SE ribbon for each combo (draw ribbon first, then line) ===
    for j = 1:numel(tags)
        [mu, se] = get_mu_se(tags{j}, pc);
        if isempty(mu) || all(isnan(mu)), continue; end
        mu = mu(:);  se = se(:);           % ensure column vectors
        
        % --- minimal fix: flip Older-Global only for PC1 ---
        if pc == 1 && strcmp(tags{j}, 'Old_Global')
            mu = -mu;   % sign flip
        end

        lo = mu - se;  hi = mu + se;       % SE band
    
        % Ribbon (shaded SE)
        patch(ax, [t_plot(:); flipud(t_plot(:))], ...
                  [lo;        flipud(hi)], ...
                  comb_colors(j,:), ...
                  'FaceAlpha', ribbon_alpha, ...
                  'EdgeColor', 'none', ...
                  'HandleVisibility', 'off', ...
                  'Clipping', 'on');
    
        % Mean line on top
        plot(ax, t_plot, mu, 'LineWidth', 1.2, 'Color', comb_colors(j,:));
    end

    % === D) Bottom indicator lines using the same grid-aligned windows ===
    draw_bottom_lines_7(ax, {wins_G, wins_C, wins_I}, {clr_group, clr_cond, clr_int}, ...
                        xlim_target, ymin, ymax);

    % Cosmetics
    ylim(ax, [ymin ymax]);
    xlim(ax, xlim_target);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Activation (a.u.)');
    title(ax, sprintf('PC%d', pc));

    % Dotted vertical reference lines
    for vv = 1:numel(vlines)
        if vlines(vv) >= xlim_target(1) && vlines(vv) <= xlim_target(2)
            xline(ax, vlines(vv), ':', 'Color', [0 0 0], 'LineWidth', 0.75, 'HandleVisibility','off');
        end
    end
end

% Legend tile
axL = nexttile(tlo, 4); cla(axL); hold(axL,'on'); axis(axL,'off');

% curve entries (4)
h_leg = gobjects(1, numel(tags));
for j = 1:numel(tags)
    h_leg(j) = plot(axL, NaN, NaN, '-', 'LineWidth', 3, 'Color', comb_colors(j,:));
end

% effect line entries (3)
hE = gobjects(1,3);
hE(1) = plot(axL, NaN, NaN, '-', 'LineWidth', 3.0, 'Color', clr_group);
hE(2) = plot(axL, NaN, NaN, '-', 'LineWidth', 3.0, 'Color', clr_cond);
hE(3) = plot(axL, NaN, NaN, '-', 'LineWidth', 3.0, 'Color', clr_int);

% ----------- KEY CHANGE: ordering + NumColumns = 4 -----------
lgd = legend(axL, ...
    [h_leg, hE], ...
    [labels, {'Group main effect','Condition main effect','Interaction effect'}], ...
    'NumColumns', 4, ...
    'Orientation','horizontal',...
    'Location','northwest');

set(lgd, 'Interpreter','none', 'Units','normalized', 'Location','none');
set(lgd, 'Position', [0.53 0.10 0.40 0.30]);
set(lgd, 'ItemTokenSize', [15, 12]);
lgd.TextColor = [0 0 0];
lgd.Color     = [1 1 1];
lgd.EdgeColor = 'none';


sgtitle(tlo, 'Brain Network Time Series — cluster-perm windows (gray), effect ticks (colored)');

out_png = fullfile(stats_outdir, 'Summary_1x3_PCwise_Using_ANOVA6_GRAY.png');
print(h, out_png, '-dpng','-r200');
close(h);
fprintf('7 saved: %s\n', out_png);

% ===================== LOCAL HELPERS (7) =====================

function [mu, se] = extract_mu_se_proper_7(B, t_idx, pc)
% Returns column vectors [numel(t_idx) x 1] for mu and se
% Expect data layout: [time x PC x cond x subj]
    X = B.TimeSeries_BrainNetworks(t_idx, pc, 1, :);  % [T x 1 x 1 x Nsubj]
    T = numel(t_idx);
    N = size(X, 4);
    X = reshape(X, T, N);                             % [T x Nsubj]
    mu = mean(X, 2, 'omitnan');                       % [T x 1]
    if N <= 1
        se = zeros(T,1);
    else
        se = std(X, 0, 2, 'omitnan') ./ sqrt(N);      % [T x 1]
    end
end

function mask_plot = sig_to_plot_mask(sig_mask, t_stat, t_plot)
% Map a logical mask defined on t_stat onto the plotting grid t_plot.
% We set a t_plot sample true if the nearest t_stat sample is significant,
% provided the distance <= half a t_plot bin.
    t_stat  = t_stat(:);
    sig_mask = logical(sig_mask(:))';
    if isempty(t_stat) || isempty(t_plot) || numel(sig_mask) ~= numel(t_stat)
        mask_plot = false(size(t_plot));
        return;
    end
    dt_plot = median(diff(t_plot));
    mask_plot = false(size(t_plot));
    % nearest neighbor assign
    for k = 1:numel(t_stat)
        if ~sig_mask(k), continue; end
        [d, j] = min(abs(t_plot - t_stat(k)));
        if d <= dt_plot/2 + 10*eps(max(abs(t_plot)))  % robust tolerance
            mask_plot(j) = true;
        end
    end
    % optional: bridge tiny gaps due to roundoff (two isolated falses between trues)
    % (kept conservative; comment out if undesired)
    % mask_plot = imclose(mask_plot, ones(1,3)); % requires Image Processing Toolbox
end

function wins = mask_to_windows_on_grid(mask_plot, t_plot, dt_plot, xlim_ok)
% Convert a logical mask on t_plot to [n x 2] windows in seconds,
% expanding by half-bin so the last marked sample is fully covered.
    wins = [];
    if ~any(mask_plot), return; end
    on = find(mask_plot(:)');
    br = [1, find(diff(on) > 1)+1, numel(on)+1];
    for b = 1:numel(br)-1
        seg = on(br(b):br(b+1)-1);
        t1  = t_plot(seg(1))  - dt_plot/2;
        t2  = t_plot(seg(end)) + dt_plot/2;
        % clip to requested x-limits
        t1 = max(t1, xlim_ok(1));
        t2 = min(t2, xlim_ok(2));
        if t2 > t1
            wins(end+1,:) = [t1, t2]; %#ok<AGROW>
        end
    end
end

function draw_gray_windows_7(ax, wins, ymin, ymax, shade_gray, shade_alpha)
% wins: [n x 2] in seconds on the plotting grid.
    if isempty(wins), return; end
    for i=1:size(wins,1)
        t1 = wins(i,1); t2 = wins(i,2);
        X = [t1 t2 t2 t1];
        Y = [ymin ymin ymax ymax];
        patch('XData',X,'YData',Y,'FaceColor',shade_gray,'FaceAlpha',shade_alpha, ...
              'EdgeColor','none','Parent',ax,'HitTest','off','Clipping','on');
    end
end

function draw_bottom_lines_7(ax, wins_cell, effect_colors, xlim_ok, ymin, ymax)
% wins_cell = {wins_GROUP, wins_COND, wins_INT}; each [n x 2] on the plotting grid.
    oldUnits = get(ax,'Units'); set(ax,'Units','pixels');
    axPos = get(ax,'Position'); set(ax,'Units',oldUnits);
    yRange = ymax - ymin;
    spacing_factor = 5;
    dy_data = (axPos(4) > 0) * (spacing_factor / max(axPos(4),1)) * yRange;
    if dy_data==0, dy_data = 0.003*yRange; end
    y_base = ymin + 0.03 * yRange;  % bottom baseline
    n_levels = 3;                   % top→bottom: Group, Condition, Interaction

    for lvl = 1:3
        wins = wins_cell{lvl};
        if isempty(wins), continue; end
        y_line = y_base + (n_levels - lvl)*dy_data;  % stack lines
        for i=1:size(wins,1)
            t1 = max(wins(i,1), xlim_ok(1));
            t2 = min(wins(i,2), xlim_ok(2));
            if t2 <= t1, continue; end
            plot(ax, [t1 t2], [y_line y_line], '-', 'LineWidth', 3.0, ...
                 'Color', effect_colors{lvl}, 'HitTest','off');
        end
    end
end


%% ===================== HELPERS (local functions) =====================

function vec = local_collect_vals(RQAcell, metric, ids, condIdx)
% Collect a vector of values for a given metric from RQA_BROADNESS.RQA_metrics
% RQAcell{p} is a table with rows: 1=Global, 2=Local
    vec = nan(numel(ids),1);
    for i = 1:numel(ids)
        p = ids(i);
        if p < 1 || p > numel(RQAcell) || isempty(RQAcell{p})
            vec(i) = NaN; continue;
        end
        T = RQAcell{p};
        if height(T) < condIdx || ~ismember(metric, T.Properties.VariableNames)
            vec(i) = NaN; continue;
        end
        val = T{condIdx, metric};
        if istable(val) || iscell(val), val = val{:}; end
        vec(i) = val;
    end
end

function draw_violin_mean_SD(ax, data, xpos, colorRGB)
% Draw a single violin at xpos with given color, overlay mean ± SD and jittered points
    data = data(:);
    data = data(isfinite(data));
    if isempty(data)
        return;
    end

    % Kernel density for violin
    try
        [f, yi] = ksdensity(data, 'Function','pdf');
    catch
        % Fallback if ksdensity not available
        yi = linspace(min(data), max(data), 100);
        f  = histcounts(data, [yi-1e-9, yi(end)+1e-9], 'Normalization','pdf');
        f  = [f(1), f, f(end)];
        yi = linspace(min(data), max(data), numel(f));
    end

    % Normalize width
    f = f ./ max(f);
    width = 0.35; % half-width of violin
    x_left  = xpos - f * width;
    x_right = xpos + f * width;

    % Violin patch
    xv = [x_left, fliplr(x_right)];
    yv = [yi, fliplr(yi)];
    ph = patch(ax, xv, yv, colorRGB, 'FaceAlpha', 0.25, 'EdgeColor', colorRGB*0.6, 'LineWidth', 1);

    % Jittered points
    n = numel(data);
    jitter = (rand(n,1)-0.5) * 0.20; % symmetric jitter
    scatter(ax, xpos + jitter, data, 14, 'MarkerFaceColor', colorRGB, 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.6);

    % Mean ± SD
    m  = mean(data);
    sd = std(data, 0);
    plot(ax, [xpos-0.22 xpos+0.22], [m m], '-', 'Color', colorRGB*0.8, 'LineWidth', 3);      % mean bar
    plot(ax, [xpos xpos], [m-sd m+sd], '-', 'Color', colorRGB*0.8, 'LineWidth', 2);          % SD whisker
    plot(ax, [xpos-0.10 xpos+0.10], [m-sd m-sd], '-', 'Color', colorRGB*0.8, 'LineWidth', 2);% SD caps
    plot(ax, [xpos-0.10 xpos+0.10], [m+sd m+sd], '-', 'Color', colorRGB*0.8, 'LineWidth', 2);
end



%% FINAL REMARKS

% Please check the BROADNESS GitHub repository for new releases.  
% https://github.com/leonardob92/BROADNESS_MEG_AuditoryRecognition/tree/main/BROADNESS_Toolbox
% Feel free to reach out to us if you need guidance or consultation.  
% Leonardo Bonetti: leonardo.bonetti@clin.au.dk
%                   leonardo.bonetti@psych.ox.ac.uk
% Mattia Rosso:     mattia.rosso@clin.au.dk
% Chiara Malvaso:   chiara.malvaso@studio.unibo.it
%
%  Please cite the first BROADNESS paper if using the BROADNESS toolbox:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science.
%  BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
%  https://doi.org/10.1002/advs.202507878