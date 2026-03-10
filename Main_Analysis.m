% 
% ========================================================================
%  The Ageing Brain Elicits Increased Hierarchical Prediction Error Signalling in Auditory and Cingulate Brain Networks
%
%  Please cite the preprint:
%  Authors. Archive server.
%  Title.
%  https://doi.org/X
%
% ========================================================================
%
%  This script demonstrates one application of BROADNESS 
%  (Broadband Brain Network Estimation via Source Separation) on
%  source-reconstructed MEG data from the auditory local-global paradigm.
%  
%  When using the BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  please cite the first BROAD-NESS paper:
%
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science.
%  BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
%  https://doi.org/10.1002/advs.202507878
%  
% ========================================================================
%  The dataset used in the present script is available at the following link:
%  https://doi.org/10.5281/zenodo.18231641
%
%  The code for the experimental paradigm, the auditory local-global paradigm, used in
%  this study is available at the following link:
%  https://doi.org/10.5281/zenodo.18231506
%
%  PLEASE NOTE
%  The current BROADNESS pipeline takes as input a single 4D matrix that
%  includes all individual participants’ data, computes the group average
%  for network estimation (using PCA), and then outputs the corresponding
%  time series for each participant for statistical analysis.
%  The spatial activation patterns of the networks is instead provided for
%  the group level.
%  Note that if you wish to only have a quick computation of the data averaged across
%  participants, you can use the data averaged across participants available
%  in "MMNSubtracted_Average_SignFixed.mat".
%
% ========================================================================

%
%  The script loads the MEG dataset into the MATLAB workspace and 
%  calls four core functions plus a start-up:
%
% ------------------------------------------------------------------------
%  FUNCTIONS OVERVIEW:
% ------------------------------------------------------------------------
%  - 0) BROADNESS_Startup() 
%       Initializes the environment.
%
%  - 1) BROADNESS_NetworkEstimation() 
%       Performs PCA on the event-related field/potentials to derive the
%       underlying brain networks.
%
%  - 2) BROADNESS_Visualizer() 
%       Takes as input selected outputs from `BROADNESS_NetworkEstimation` 
%       to visualize a number of features such as brain network time series
%       and topographies.
%
%  - 3) BROADNESS_PhaseSpace_RQA()
%       Computes phase space embedding and RQA on brain networks time
%       series [it works for both the time series extracted by
%       BROADNESS_NetworkEstimation (PCA) and BROADNESS_AlternativeNetworkEstimation_ICA (ICA)]
%
%  - 4) BROADNESS_SpatialGradients()
%       Computes spatial gradients embedding and clustering on the spatial
%       activation patterns of the brain networks

% -------------------------------------------------------------------------
%  AUTHOR OF THIS SCRIPT:
%  Mathias Houe Andersen
%  mathiasha@drcmr.dk
%  Danish Research Centre for Magnetic Resonance, Hvidovre University
%  Hospital
%  Faculty of Health and Medical Sciences, University of Copenhagen,
%  Copenhagen, Denmark.
%  Updated version 28/01/2026
% ========================================================================

% ------------------------------------------------------------------------
%  AUTHORS OF THE BROADNESS TOOLBOX:
%  Leonardo Bonetti, Chiara Malvaso, Mathias Houe Andersen & Mattia Rosso
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  chiara.malvaso@studio.unibo.it
%  mathiasha@drcmr.dk
%  mattia.rosso@clin.au.dk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Danish Research Centre for Magnetic Resonance, Hvidovre University
%  Hospital
%  Faculty of Health and Medical Sciences, University of Copenhagen
%  Department of Physics, University of Bologna
%  Aarhus (DK), Oxford (UK), Copenhagen (DK), Bologna (Italy)
% ========================================================================

%% 0) STARTUP

% Simply download the BROADNESS Toolbox folder and place it in the working directory,
% Make sure not to alter the structure of its functions, subfolders, or files.

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

% Load data
load('/main_path/Data/Mathias_MMN/MMNSubtracted_Average_SignFixed.mat');
DATA = dum;

% Remove first participant and define time vector with a 100 ms baseline
time = -0.100:0.004:0.8;
DATA = DATA(:,101:326,:,2:78); 

size(DATA)

% Load groups and adjust subid after removing first participant
load('/main_path/Data/Mathias_MMN/groups.mat');
older = older - 1;
young = young - 1;

older_subj = older(:);   % row-wise versions
young_subj = young(:);

%%% ------------------ COMPUTATION --------------------- %%%

% Run BROADNESS network estimation (default parameters)
BROADNESS = BROADNESS_NetworkEstimation(DATA, time);

%% 1b) Split data into 4 different groups

cond_global = 1;
cond_local  = 2;

% Copy the struct
BROADNESS_young_global = BROADNESS;
BROADNESS_young_local  = BROADNESS;
BROADNESS_older_global = BROADNESS;
BROADNESS_older_local  = BROADNESS;

% Slice TimeSeries
BROADNESS_young_global.TimeSeries_BrainNetworks = BROADNESS.TimeSeries_BrainNetworks(:, :, cond_global, young_subj);
BROADNESS_young_local.TimeSeries_BrainNetworks  = BROADNESS.TimeSeries_BrainNetworks(:, :, cond_local,  young_subj);
BROADNESS_older_global.TimeSeries_BrainNetworks = BROADNESS.TimeSeries_BrainNetworks(:, :, cond_global, older_subj);
BROADNESS_older_local.TimeSeries_BrainNetworks  = BROADNESS.TimeSeries_BrainNetworks(:, :, cond_local,  older_subj);

% Slice OriginalData to keep dimensions consistent
BROADNESS_young_global.OriginalData = BROADNESS.OriginalData(:, :, cond_global, young_subj);
BROADNESS_young_local.OriginalData  = BROADNESS.OriginalData(:, :, cond_local,  young_subj);
BROADNESS_older_global.OriginalData = BROADNESS.OriginalData(:, :, cond_global, older_subj);
BROADNESS_older_local.OriginalData  = BROADNESS.OriginalData(:, :, cond_local,  older_subj);

% Sanity checks by displaying sizes
assert(size(BROADNESS_young_local.TimeSeries_BrainNetworks,3) == size(BROADNESS_young_local.OriginalData,3), 'cond mismatch');
assert(size(BROADNESS_young_local.TimeSeries_BrainNetworks,4) == size(BROADNESS_young_local.OriginalData,4), 'subj mismatch');

% Display sizes
disp(size(BROADNESS.TimeSeries_BrainNetworks))
disp(size(BROADNESS_young_global.TimeSeries_BrainNetworks))
disp(size(BROADNESS_young_local.TimeSeries_BrainNetworks))
disp(size(BROADNESS_older_global.TimeSeries_BrainNetworks))
disp(size(BROADNESS_older_local.TimeSeries_BrainNetworks))

%% 1c) %% Calculate ED on BROADNESS, PCA applied to all data

% Assumes:
% BROADNESS.Variance_BrainNetworks : 225 x 1 double

% Extract eigenvalue / variance vector
lam = BROADNESS.Variance_BrainNetworks(:);

% Enforce non-negativity
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

% Effective Dimensionality (ED) = 2.835

%% 1d) Plot variance explained 

% Plot variance explained for first 10 PCs from BROADNESS.Variance_BrainNetworks

% Styling:
% - dark purple line + markers
% - black dotted vertical line at x = 2.835
% - Helvetica everywhere
% - boxed legend inside plot

% -----------------------------
% Data prep (robust to row/col)
% -----------------------------

% -----------------------------
% Styling controls
% -----------------------------
legendFontSize = 10;      % resize legend text
gridAlpha      = 0.15;    % transparency of grid lines (0–1)
gridColor      = [0 0 0]; % grid color (black, light via alpha)


ve = BROADNESS.Variance_BrainNetworks(:);   % force column vector

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

grid(ax,'on');
ax.GridLineStyle = '-';
ax.GridAlpha     = gridAlpha;
ax.GridColor     = gridColor;

% Optional: subtle minor grid (often too busy—use cautiously)
% ax.MinorGridLineStyle = ':';
% ax.MinorGridAlpha     = gridAlpha * 0.7;
% ax.XMinorGrid = 'on';
% ax.YMinorGrid = 'on';


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

xEff = 2.835;
hEff = xline(ax, xEff, ':k', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Effective Dimensionality');

% -----------------------------
% Labels / axes styling
% -----------------------------

xlabel(ax, 'Principal Component', 'FontName','Helvetica');
ylabel(ax, 'Variance Explained',  'FontName','Helvetica');

xlim(ax, [0.5, nPC + 0.5]);
xticks(ax, 1:nPC);

set(ax, ...
    'FontName','Helvetica', ...
    'Box','on', ...
    'LineWidth', 1.0, ...
    'TickDir','out');

% Optional y-limits
% ylim(ax, [min(0, min(y)*0.95), max(y)*1.05]);

% -----------------------------
% Legend (boxed, inside plot)
% -----------------------------
lgd = legend(ax, hEff, 'Effective Dimensionality', ...
    'Location','northeast');

set(lgd, ...
    'Box','on', ...
    'FontName','Helvetica', ...
    'FontSize', legendFontSize, ...
    'ItemTokenSize',[12 8]);   % controls marker/line spacing

hold(ax,'off');


%% 2a) BROADNESS VISUALIZATION

% NOTE: Either 2a) or 2b) should run, as 2b) is simply an extended version
% of 2a) with additional settings.

%   This section generates 5 plots, described as follows:
%   #1) Dynamic brain activity map of the original data
%   #2) Variance explained by the networks
%   #3) Time series of the networks
%   #4) Activation patterns of the networks (3D)
%   #5) Activation patterns of the networks (nifti images)

%%% ------------------- USER SETTINGS ------------------- %%%

% Common options
Options = [];
Options.name_nii = '/main_path/Output';
load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); 
Options.MNI_coords = MNI8;

% Run each of the 5 lines below in seperate

% All data with label
% Options.Labels = {'All data'}; BROADNESS_Visualizer(BROADNESS, Options);

% 4 Groups with different labels
Options.Labels = {'Younger adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_young_local, Options);
%Options.Labels = {'Older adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_older_local, Options); 
%Options.Labels = {'Younger adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_young_global, Options);
%Options.Labels = {'Older adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_older_global, Options);

%% 2b) BROADNESS VISUALIZATION (ALTERNATIVE SCENARIO WITH OPTIONAL INPUTS)

% This section demonstrates the same function as above,  
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

% If one wishes to remove the cerebellum voxels (not included in 3D brain template (#4)), please set 'remove_cerebellum_label' to 1
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
BROADNESS_Visualizer(BROADNESS,Options)

%Options.Labels = {'Younger adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_young_local, Options);
%Options.Labels = {'Older adults - Local effect'}; BROADNESS_Visualizer(BROADNESS_older_local, Options); 
%Options.Labels = {'Younger adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_young_global, Options);
%Options.Labels = {'Older adults - Global effect'}; BROADNESS_Visualizer(BROADNESS_older_global, Options);


%% 2c) - Convert MNI coords into brain labels for positively and negatively contributing voxels independently

% For each network ActivationTable_*.xlsx:
%   1) split voxels into positive vs negative contributions
%   2) map each voxel's MNI (mm) coordinate to an atlas label
%   3) count voxels per label and output sorted lists
%   4) ALSO: if a voxel lands on background (0) or unlabeled, assign it the
%            nearest label within a growing voxel-radius (0..maxSearchRadius)
%   5) Output per-area counts split by the radius used (exact, 1-away, 2-away, ...)
%
% IMPORTANT NOTES / LIMITATIONS (stress-test):
% - This forces a label by snapping to the nearest non-zero atlas label within a voxel radius.
% - Distance here is in *atlas voxel units* using a Chebyshev “ring radius”
%   (max(|dx|,|dy|,|dz|)). This matches cubic neighborhood shells:
%   r=1 -> 26-neighborhood shell, r=2 -> 5x5x5 shell boundary, etc.
% - If many voxels snap at large radii, interpret with caution: you’re labeling
%   “nearest parcel”, not necessarily the true anatomical assignment.
%
% FIXES INCLUDED:
% - Robust LUT parser for AAL-style .txt files with lines like: "1 Precentral_L 1"
%   (handles arbitrary whitespace and even missing last numeric id)
% - Fixed operator precedence bug in guessMNIColumns()

clear; clc;

% ----------------------- USER SETTINGS -----------------------
networkFiles = {
    '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/ActivationTable_BrainNetwork_1.xlsx'
    '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/ActivationTable_BrainNetwork_2.xlsx'
    '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/ActivationTable_BrainNetwork_3.xlsx'
};

% Choose an atlas in MNI space (NIfTI) + a label lookup (CSV/TXT).
% The AAL Atlas (version 3.1) can be downloaded here: https://www.oxcns.org/aal3.html
atlasNiftiPath  = '/main_path/AAL3/AAL3v1.nii';
atlasLabelsPath = '/main_path/AAL3/AAL3v1.nii.txt';

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

    % Save results
    [~, baseName, ~] = fileparts(xlsxPath);
    outXlsx = fullfile(outDir, sprintf('%s_voxelCountsByArea.xlsx', baseName));

    writetable(posTableByR, outXlsx, 'Sheet', 'Positive_ByRadius', 'WriteMode', 'overwritesheet');
    writetable(negTableByR, outXlsx, 'Sheet', 'Negative_ByRadius', 'WriteMode', 'overwritesheet');

    writetable(posRadiusSummary, outXlsx, 'Sheet', 'Positive_RadiusSummary', 'WriteMode', 'overwritesheet');
    writetable(negRadiusSummary, outXlsx, 'Sheet', 'Negative_RadiusSummary', 'WriteMode', 'overwritesheet');

    % Optional: keep simple totals as additional sheets
    writetable(posTableSimple, outXlsx, 'Sheet', 'Positive_SimpleTotals', 'WriteMode', 'overwritesheet');
    writetable(negTableSimple, outXlsx, 'Sheet', 'Negative_SimpleTotals', 'WriteMode', 'overwritesheet');

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

fprintf('\nDone. Outputs saved to: %s\n', outDir);

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
    %
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


%% 2d) TIMESERIES ANALYSIS: Statistics — Timepoint-wise mixed ANOVA via CLUSTER-BASED permutation (over time)
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
    YL = squeeze(BROADNESS_young_local. TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_younger x T]
    YG = squeeze(BROADNESS_young_global.TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_younger x T]
    OL = squeeze(BROADNESS_older_local. TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_older   x T]
    OG = squeeze(BROADNESS_older_global.TimeSeries_BrainNetworks(t_idx, pc, 1, :))';  % [N_older   x T]

    % ---- Defensive shape checks (fail early if misaligned) ----
    if any([size(YL,2), size(YG,2), size(OL,2), size(OG,2)] ~= Twin)
        error('Time dimension mismatch after slicing for PC %d.', pc);
    end
    if size(YL,1) ~= size(YG,1)
        error('Younger group Local/Global subject count mismatch for PC %d.', pc);
    end
    if size(OL,1) ~= size(OG,1)
        error('Older group Local/Global subject count mismatch for PC %d.', pc);
    end

    % ---- Compute OBSERVED effects (each is a 1 x T vector) ----
    % Group main (between-subject): average Local+Global within each subject, then Young minus Old.
    Y_mean = (YL + YG) ./ 2;  % [N_younger x T]
    O_mean = (OL + OG) ./ 2;  % [N_older   x T]
    obs_G  = mean(Y_mean, 1, 'omitnan') - mean(O_mean, 1, 'omitnan');  % 1 x T

    % Condition main (within-subject): (Global - Local) pooled across both groups.
    D_y   = YG - YL;                  % [N_younger x T] within-subject difference in Young
    D_o   = OG - OL;                  % [N_older   x T] within-subject difference in Old
    D_all = [D_y; D_o];               % [N_total x T] stack both groups
    obs_C = mean(D_all, 1, 'omitnan');% 1 x T pooled within-subject effect

    % Interaction: "difference of differences" = (Younger (G-L)) - (Older (G-L)).
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

        plot(ax2, t_win, mean(Y_mean,1,'omitnan'), 'LineWidth',1.5, 'Color',[0.2 0.2 0.8]); % Younger mean
        plot(ax2, t_win, mean(O_mean,1,'omitnan'), 'LineWidth',1.5, 'Color',[0.8 0.2 0.2]); % Older mean

        yl2 = ylim(ax2);
        shade_clusters_union(ax2, t_win, sig_G, yl2, [0.7 0.7 0.7], 0.15);

        title(ax2, sprintf('PC%d — Group means (Young vs Old) with Group-effect clusters (gray)', pc), 'Interpreter','none');
        xlabel(ax2,'Time (s)'); ylabel(ax2,'Activation (a.u.)');
        legend(ax2, {'Young (mean of Local+Global)','Older (mean of Local+Global)'}, 'Location','northwest');

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
        dY = mean(D_y,1,'omitnan');  % Younger (G-L)
        dO = mean(D_o,1,'omitnan');  % Older   (G-L)
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

    % Observed difference (Younger-Older or D_y - D_o depending on call site)
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

%% 2e) TIMESERIES ANALYSIS: Plotting — 2×3 PC-WISE SUMMARY PLOT using ANOVA results from 2d)

% Layout:
%   [ PC1 | PC2 | PC3
%   | LEGEND ]

% Uses:
%   - BROADNESS_young_local/global, BROADNESS_older_local/global
%   - stats_outdir and RES files created in 6a)
%   - 'time' vector in workspace

% Output:
%   - PNG saved in stats_outdir

stats_outdir = '/main_path/Output/ANOVA_Stats';

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

combo_colors = [ ...
    0.1 0.2 0.5;  % Young-Local   (blue)
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
ymin = -1100; ymax = 1000;

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
        lo = mu - se;  hi = mu + se;       % SE band
    
        % Ribbon (shaded SE)
        patch(ax, [t_plot(:); flipud(t_plot(:))], ...
                  [lo;        flipud(hi)], ...
                  combo_colors(j,:), ...
                  'FaceAlpha', ribbon_alpha, ...
                  'EdgeColor', 'none', ...
                  'HandleVisibility', 'off', ...
                  'Clipping', 'on');
    
        % Mean line on top
        plot(ax, t_plot, mu, 'LineWidth', 1.2, 'Color', combo_colors(j,:));
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
    h_leg(j) = plot(axL, NaN, NaN, '-', 'LineWidth', 3, 'Color', combo_colors(j,:));
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

%% 2f) Investigating ERFs: Latency, order and amplitude

% We compute peaks on the MEAN time series (mu) returned by extract_mu_se_proper_7
% for each PC and each plot (Young/Older x Local/Global).

% Local:
%   P50  : most positive in [0.05, 0.10]
%   MMN  : most negative in [0.10, 0.20]
%   P300 : most positive in [0.15, 0.30]
%   RON  : most negative in [0.25, 0.40]

% Global:
%   MMN  : most negative in [0.10, 0.20]
%   P300 : most positive in [0.10, 0.25]

% Helper: get peak (latency, amplitude) in a window from (time, signal)
% polarity: 'pos' (max) or 'neg' (min)
get_peak = @(tt, xx, t1, t2, polarity) local_get_peak(tt, xx, t1, t2, polarity);

% Windows definition
local_windows = { ...
    'P50',  0.05, 0.10, 'pos'; ...
    'MMN',  0.10, 0.20, 'neg'; ...
    'P300', 0.15, 0.30, 'pos'; ...
    'RON',  0.25, 0.40, 'neg' ...
    };

global_windows = { ...
    'MMN',  0.10, 0.20, 'neg'; ...
    'P300', 0.15, 0.25, 'pos'; ...
    'RON',  0.25, 0.40, 'neg' ...
    };

% --------- NEW: store detected peaks for later plotting ---------
PEAKS = struct();

fprintf('\n================ PEAK SUMMARY (mu) ================\n');
for it = 1:numel(tags)
    tag = tags{it};

    isLocal  = contains(tag, 'Local',  'IgnoreCase', true);
    isGlobal = contains(tag, 'Global', 'IgnoreCase', true);

    if ~(isLocal || isGlobal)
        continue;
    end

    fprintf('\n---- %s ----\n', tag);

    for ip = 1:numel(pc_list)
        pc = pc_list(ip);

        [mu, ~] = get_mu_se(tag, pc);
        if isempty(mu) || all(isnan(mu))
            fprintf('PC%d: mu empty/NaN -> skipping\n', pc);
            continue;
        end
        mu = mu(:);

        fprintf('PC%d\n', pc);

        if isLocal
            W = local_windows;
        else
            W = global_windows;
        end

        for iw = 1:size(W,1)
            name = W{iw,1};
            t1   = W{iw,2};
            t2   = W{iw,3};
            pol  = W{iw,4};

            % --------- NEW: skip RON for (Young_Local, PC2) and (Old_Local, PC2) ---------
            if strcmpi(name, 'RON') && pc == 2 && (strcmpi(tag, 'Young_Local') || strcmpi(tag, 'Old_Local'))
                continue;
            end

            % --------- NEW: skip P300 for (Young_Global, PC2) and (Old_Global, PC2) ---------
            if strcmpi(name, 'P300') && pc == 2 && (strcmpi(tag, 'Young_Global') || strcmpi(tag, 'Old_Global'))
                continue;
            end

            % --------- NEW: Global RON only for PC1 ---------
            if strcmpi(name, 'RON') && isGlobal && pc ~= 1
                continue;
            end

            [lat, amp] = get_peak(t_plot(:), mu, t1, t2, pol);

            fprintf('%s  %.6f  %.6f\n', name, lat, amp);

            % --------- NEW: keep for plotting (per tag, PC, component) ---------
            if ~isfield(PEAKS, tag)
                PEAKS.(tag) = struct();
            end
            pc_field = sprintf('PC%d', pc);
            if ~isfield(PEAKS.(tag), pc_field)
                PEAKS.(tag).(pc_field) = struct();
            end
            PEAKS.(tag).(pc_field).(name) = [lat, amp];
        end
    end
end
fprintf('===================================================\n');

% --------- EXPORT PEAKS TO EXCEL (4 sheets) ---------
out_xlsx = fullfile(pwd, 'Peak_Summary.xlsx');

sheet_tags = {'Young_Local','Old_Local','Young_Global','Old_Global'};

for is = 1:numel(sheet_tags)
    tag = sheet_tags{is};

    isLocal  = contains(tag, 'Local',  'IgnoreCase', true);
    isGlobal = contains(tag, 'Global', 'IgnoreCase', true);

    if isLocal
        erf_list = {'P50','MMN','P300','RON'};
    elseif isGlobal
        erf_list = {'MMN','P300','RON'};
    else
        erf_list = {};
    end

    nPC = numel(pc_list);
    varNames = cell(1, 1 + 3*nPC);
    varNames{1} = 'ERF';
    c = 2;
    for k = 1:nPC
        if k == 1
            pcName = 'First_PC';
        elseif k == 2
            pcName = 'Second_PC';
        elseif k == 3
            pcName = 'Third_PC';
        else
            pcName = sprintf('%dth_PC', k);
        end
        varNames{c} = pcName; c = c + 1;
        varNames{c} = sprintf('Latency_%d', k); c = c + 1;
        varNames{c} = sprintf('Amplitude_%d', k); c = c + 1;
    end

    nERF = numel(erf_list);
    rows = cell(nERF, 1 + 3*nPC);

    for ie = 1:nERF
        erf = erf_list{ie};

        lat_all = nan(nPC,1);
        amp_all = nan(nPC,1);
        pc_all  = nan(nPC,1);

        for ip = 1:nPC
            pc = pc_list(ip);
            pc_field = sprintf('PC%d', pc);

            lat = NaN; amp = NaN;

            if isfield(PEAKS, tag) && isfield(PEAKS.(tag), pc_field) && isfield(PEAKS.(tag).(pc_field), erf)
                v = PEAKS.(tag).(pc_field).(erf);
                if numel(v) >= 2
                    lat = v(1);
                    amp = v(2);
                end
            end

            lat_all(ip) = lat;
            amp_all(ip) = amp;
            pc_all(ip)  = pc;
        end

        valid = ~isnan(lat_all);
        pc_v  = pc_all(valid);
        lat_v = lat_all(valid);
        amp_v = amp_all(valid);

        if ~isempty(lat_v)
            M = [lat_v(:), -amp_v(:)];
            [~, idx] = sortrows(M, [1 2]);
            pc_sorted  = pc_v(idx);
            lat_sorted = lat_v(idx);
            amp_sorted = amp_v(idx);
        else
            pc_sorted  = [];
            lat_sorted = [];
            amp_sorted = [];
        end

        row = cell(1, 1 + 3*nPC);
        row{1} = erf;

        col = 2;
        for k = 1:numel(pc_sorted)
            row{col} = sprintf('PC%d', pc_sorted(k)); col = col + 1;
            row{col} = lat_sorted(k); col = col + 1;
            row{col} = amp_sorted(k); col = col + 1;
        end

        while col <= (1 + 3*nPC)
            row{col} = ''; col = col + 1;
            if col <= (1 + 3*nPC); row{col} = NaN; col = col + 1; end
            if col <= (1 + 3*nPC); row{col} = NaN; col = col + 1; end
        end

        rows(ie,:) = row;
    end

    T = cell2table(rows, 'VariableNames', varNames);
    writetable(T, out_xlsx, 'Sheet', tag);
end


% --------- local function  ---------
function [lat, amp] = local_get_peak(tt, xx, t1, t2, polarity)
    idx = find(tt >= t1 & tt <= t2);
    if isempty(idx) || isempty(xx)
        lat = NaN; amp = NaN;
        return;
    end
    xw = xx(idx);

    if all(isnan(xw))
        lat = NaN; amp = NaN;
        return;
    end

    switch lower(polarity)
        case 'pos'
            [amp, k] = max(xw);
        case 'neg'
            [amp, k] = min(xw);
        otherwise
            error('local_get_peak: polarity must be ''pos'' or ''neg''.');
    end

    lat = tt(idx(k));
end


stats_outdir = '/main_path/Output/ANOVA_Stats';

if ~exist('time','var') || isempty(time)
    error('7: time vector not found in workspace.');
end

tags   = {'Young_Local','Old_Local','Young_Global','Old_Global'};
labels = {'Younger - Local','Older - Local','Younger - Global','Older - Global'};

combo_colors = [ ...
    0.1 0.2 0.5;
    0.5 0.1 0.1;
    0.4 0.6 0.85;
    0.8 0.2 0.2];

pc_list   = [1 2 3];
pc_styles = {'-','--',':'};
pc_names  = {'PC #1','PC #2','PC #3'};

ymin = -1100; ymax = 1000;

xlim_target = [-0.1 0.8];
vlines = [0.0];

t_plot_idx = find(time >= xlim_target(1) & time <= xlim_target(2));
t_plot     = time(t_plot_idx);
if isempty(t_plot)
    error('7: No samples within [%g, %g] in the time vector.', xlim_target(1), xlim_target(2));
end

BMAP = containers.Map( ...
    tags, ...
    {BROADNESS_young_local, BROADNESS_older_local, BROADNESS_young_global, BROADNESS_older_global} ...
);

get_mu_se = @(tag, pc) extract_mu_se_proper_7(BMAP(tag), t_plot_idx, pc);

marker_map = containers.Map( ...
    {'P50','MMN','P300','RON'}, ...
    {'x','o','^','s'} ...
);

h  = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1400 1000]);

tlo = tiledlayout(2,2, 'Padding','compact','TileSpacing','compact');

axPos = cell(1,4);
for k = 1:4
    axPos{k} = nexttile(tlo, k);
end
for k = 1:4
    cla(axPos{k}); hold(axPos{k},'on');
end

left_col_x   = 0.07;
right_col_x  = 0.39;
col_w        = 0.24;
top_row_y    = 0.56;
bot_row_y    = 0.18;
row_h        = 0.30;

set(axPos{1}, 'Position', [left_col_x,  top_row_y, col_w, row_h]);
set(axPos{2}, 'Position', [right_col_x, top_row_y, col_w, row_h]);
set(axPos{3}, 'Position', [left_col_x,  bot_row_y, col_w, row_h]);
set(axPos{4}, 'Position', [right_col_x, bot_row_y, col_w, row_h]);

ribbon_alpha = 0.20;

lgd_lw   = 3.0;
lgd_fs   = 14;
lgd_loc  = 'southeast';

for it = 1:numel(tags)
    tag = tags{it};
    ax  = axPos{it};
    hold(ax,'on');

    isLocal  = contains(tag, 'Local',  'IgnoreCase', true);
    isGlobal = contains(tag, 'Global', 'IgnoreCase', true);

    for ip = 1:numel(pc_list)
        pc = pc_list(ip);
        [mu, se] = get_mu_se(tag, pc);
        if isempty(mu) || all(isnan(mu)), continue; end
        mu = mu(:);
        se = se(:);

        lo = mu - se;
        hi = mu + se;

        patch(ax, [t_plot(:); flipud(t_plot(:))], ...
                  [lo;        flipud(hi)], ...
                  combo_colors(it,:), ...
                  'FaceAlpha', ribbon_alpha, ...
                  'EdgeColor', 'none', ...
                  'HandleVisibility', 'off', ...
                  'Clipping', 'on');

        plot(ax, t_plot, mu, pc_styles{ip}, ...
            'LineWidth', 1.2, ...
            'Color', combo_colors(it,:));

        if isLocal
            W = local_windows;
        else
            W = global_windows;
        end

        pc_field = sprintf('PC%d', pc);
        if isfield(PEAKS, tag) && isfield(PEAKS.(tag), pc_field)
            for iw = 1:size(W,1)
                nm = W{iw,1};

                if strcmpi(nm, 'RON') && pc == 2 && (strcmpi(tag, 'Young_Local') || strcmpi(tag, 'Old_Local'))
                    continue;
                end

                if strcmpi(nm, 'P300') && pc == 2 && (strcmpi(tag, 'Young_Global') || strcmpi(tag, 'Old_Global'))
                    continue;
                end

                if isfield(PEAKS.(tag).(pc_field), nm)
                    pa = PEAKS.(tag).(pc_field).(nm);
                    lat = pa(1);
                    amp = pa(2);
                    if ~isnan(lat) && ~isnan(amp)
                        mk = 'x';
                        if isKey(marker_map, nm)
                            mk = marker_map(nm);
                        end

                        if strcmp(mk, 'x')
                            plot(ax, lat, amp, mk, ...
                                'Color', [0 0 0], ...
                                'LineWidth', 1.2, ...
                                'MarkerSize', 7, ...
                                'HandleVisibility', 'off', ...
                                'Clipping', 'on');
                        else
                            plot(ax, lat, amp, mk, ...
                                'Color', [0 0 0], ...
                                'MarkerFaceColor', 'none', ...
                                'LineWidth', 1.2, ...
                                'MarkerSize', 7, ...
                                'HandleVisibility', 'off', ...
                                'Clipping', 'on');
                        end
                    end
                end
            end
        end
    end

    ylim(ax, [ymin ymax]);
    xlim(ax, xlim_target);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Activation (a.u.)');
    title(ax, labels{it}, 'Interpreter','none');

    for vv = 1:numel(vlines)
        if vlines(vv) >= xlim_target(1) && vlines(vv) <= xlim_target(2)
            xline(ax, vlines(vv), ':', 'Color', [0 0 0], 'LineWidth', 0.75, 'HandleVisibility','off');
        end
    end

    h1 = plot(ax, nan, nan, '-',  'LineWidth', lgd_lw, 'Color', combo_colors(it,:));
    h2 = plot(ax, nan, nan, '--', 'LineWidth', lgd_lw, 'Color', combo_colors(it,:));
    h3 = plot(ax, nan, nan, ':',  'LineWidth', lgd_lw, 'Color', combo_colors(it,:));

    hP50  = plot(ax, nan, nan, 'x', 'Color', [0 0 0], 'LineWidth', 1.2, 'MarkerSize', 7);
    hMMN  = plot(ax, nan, nan, 'o', 'Color', [0 0 0], 'LineWidth', 1.2, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
    hP300 = plot(ax, nan, nan, '^', 'Color', [0 0 0], 'LineWidth', 1.2, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
    hRON  = plot(ax, nan, nan, 's', 'Color', [0 0 0], 'LineWidth', 1.2, 'MarkerSize', 7, 'MarkerFaceColor', 'none');

    % Bottom row (plots 3 and 4): remove "P50" from legend
    if it == 3 || it == 4
        lgd = legend(ax, [h1 h2 h3 hMMN hP300 hRON], ...
            {'PC #1','PC #2','PC #3','MMN','P300','RON'}, 'Location', lgd_loc);
    else
        lgd = legend(ax, [h1 h2 h3 hP50 hMMN hP300 hRON], ...
            {'PC #1','PC #2','PC #3','P50','MMN','P300','RON'}, 'Location', lgd_loc);
    end
    % ========================================================

    set(lgd, 'FontName', 'Helvetica', 'FontSize', lgd_fs);
end

sgtitle(tlo, 'Brain Network Time Series — PC1/PC2/PC3 as line styles (no stats)');

out_png = fullfile(stats_outdir, 'Summary_2x2_ComboWise_PCstyles_NoStats.png');
print(h, out_png, '-dpng','-r200');
close(h);
fprintf('7 saved: %s\n', out_png);


%% 3a) PHASE SPACE EMBEDDING AND RECURRENCE QUANTIFICATION ANALYSIS 

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply use the structure outputted by the BROADNESS_NetworkEstimation function
% Additional optional inputs can be provided, as described in the function. 
% For 3D visualization: 'principalcomps',[1:3] below
% For 2D visualization: 'principalcomps',[1:2] below
% For switching between individual x, y, z limits and limits standard across all plots, edit 224-225 (2D) or line
% 257-259 (3D) in BROADNESS_PhaseSpace_RQA.m: Insert commented lines for
% individual limits, and defined limits for standard limits across all plots.

%%% ------------------ COMPUTATION --------------------- %%%

% Shift here from 3D to 2D illustrations, by commenting out the unwanted
% section of the dimensionality


% 3D - Overall
RQA_BROADNESS = BROADNESS_PhaseSpace_RQA(BROADNESS,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'All data');

% 3D - Group/condition specific
RQA_BROADNESS_young_local = BROADNESS_PhaseSpace_RQA(BROADNESS_young_local,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yl');
RQA_BROADNESS_older_local = BROADNESS_PhaseSpace_RQA(BROADNESS_older_local,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Ol');
RQA_BROADNESS_young_global = BROADNESS_PhaseSpace_RQA(BROADNESS_young_global,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yg');
RQA_BROADNESS_older_global = BROADNESS_PhaseSpace_RQA(BROADNESS_older_global,'principalcomps',[1:3],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Og');


% 2D - Overall
% RQA_BROADNESS = BROADNESS_PhaseSpace_RQA(BROADNESS,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'All data');

% 2D - Group/condition specific
% RQA_BROADNESS_young_local = BROADNESS_PhaseSpace_RQA(BROADNESS_young_local,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yl');
% RQA_BROADNESS_older_local = BROADNESS_PhaseSpace_RQA(BROADNESS_older_local,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Ol');
% RQA_BROADNESS_young_global = BROADNESS_PhaseSpace_RQA(BROADNESS_young_global,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Yg');
% RQA_BROADNESS_older_global = BROADNESS_PhaseSpace_RQA(BROADNESS_older_global,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on', 'plotlabel', 'Og');


%% 3b) Export average phase space values

% After running 3D plots in 3a):

% Put items in a list with filenames
items = { ...
    RQA_BROADNESS_young_local, 'young_local_phase_space.csv'; ...
    RQA_BROADNESS_older_local, 'older_local_phase_space.csv'; ...
    RQA_BROADNESS_young_global, 'young_global_phase_space.csv'; ...
    RQA_BROADNESS_older_global, 'older_global_phase_space.csv'  ...
};

for k = 1:size(items,1)
    RQA  = items{k,1};
    fout = items{k,2};

    t = RQA.PhaseSpace.time(:)';          % 1 x T
    XYZ = RQA.PhaseSpace.coords{1};       % [T x 3] (assumes condition = 1)

    % Basic sanity checks (avoid silent wrong output)
    if size(XYZ,2) ~= 3
        error('Expected 3 PCs for %s, but got %d columns.', fout, size(XYZ,2));
    end
    if size(XYZ,1) ~= numel(t)
        error('Time length mismatch for %s (time=%d, coords=%d).', fout, numel(t), size(XYZ,1));
    end

    out = [ ...
        t; ...
        XYZ(:,1)'; ...
        XYZ(:,2)'; ...
        XYZ(:,3)' ...
    ];  % 4 x T

    writematrix(out, fout);
end



%% 3c) Plot differences in 3D phase-space coordinates
%  1) LOCAL:  Young_local - Older_local
%  2) GLOBAL: Young_global - Older_global
%  3) LOCAL vs GLOBAL (group-averaged): mean(Young_local,Older_local) - mean(Young_global,Older_global)
%
% Within each plot, x/y/z differences are overlaid:
%   - x: solid
%   - y: dashed
%   - z: dotted
%
% Segment coloring (UPDATED):
%   Colors no longer encode "young/local vs old/global".
%   Instead, each time point uses the same jet-based time coloring as the
%   phase space video: color progresses with time index (early->late).
%
% X-axis forced to [-0.1 0.8] seconds (data clipped to what exists).
% Manually set y-axis limits via ylims_user (set [] for auto).

% ------------------- USER SETTINGS -------------------
RQA_y_local  = RQA_BROADNESS_young_local;
RQA_o_local  = RQA_BROADNESS_older_local;

RQA_y_global = RQA_BROADNESS_young_global;
RQA_o_global = RQA_BROADNESS_older_global;

cond_to_plot = 1;            % which condition to use (1..nCond)
xlims = [-0.1 0.8];          % requested time window (s)
ylims_user = [-800 800];     % <<< MANUAL Y LIMITS (set [] to auto)
lw = 2;                      % line width

% ------------------- GLOBAL FIGURE FONT SETTINGS -------------------
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');
set(groot, 'defaultLegendFontName', 'Helvetica');

% ------------------- FIGURE 1: LOCAL (Young - Older) -------------------
figure('Color','w');
plot_one_context(RQA_y_local, RQA_o_local, cond_to_plot, xlims, ylims_user, ...
    lw, 'LOCAL (Young - Older)');

% ------------------- FIGURE 2: GLOBAL (Young - Older) ------------------
figure('Color','w');
plot_one_context(RQA_y_global, RQA_o_global, cond_to_plot, xlims, ylims_user, ...
    lw, 'GLOBAL (Young - Older)');

% ------------------- FIGURE 3: LOCAL vs GLOBAL (avg(Y,O)) --------------
figure('Color','w');
plot_local_vs_global_avg(RQA_y_local, RQA_o_local, RQA_y_global, RQA_o_global, ...
    cond_to_plot, xlims, ylims_user, lw);

% =================== LOCAL FUNCTIONS ===================

function plot_one_context(RQA_y, RQA_o, cond_to_plot, xlims, ylims_user, lw, plotTitle)

    % ----- Validate -----
    require_phasespace(RQA_y, 'RQA_y');
    require_phasespace(RQA_o, 'RQA_o');

    coords_y = RQA_y.PhaseSpace.coords;
    coords_o = RQA_o.PhaseSpace.coords;

    if cond_to_plot > numel(coords_y) || cond_to_plot > numel(coords_o)
        error('cond_to_plot exceeds available conditions.');
    end

    PS_y = coords_y{cond_to_plot};  % [T x 3]
    PS_o = coords_o{cond_to_plot};  % [T x 3]

    if size(PS_y,2) < 3 || size(PS_o,2) < 3
        error('Phase-space coordinates do not have 3 columns. Did you run with ''principalcomps'',[1:3]?');
    end

    t_y = RQA_y.PhaseSpace.time(:);
    t_o = RQA_o.PhaseSpace.time(:);

    % ----- Align time bases -----
    [t, PS_o_aligned, PS_y_aligned] = align_to_young_time(t_y, PS_y, t_o, PS_o);

    % ----- FULL colormap (EXACTLY like phase space video: jet(full_length)) -----
    cmapFull = jet(numel(t));

    % ----- Clip to requested time window -----
    % Compare magnitudes by treating numbers as positive
    % i.e., diff = abs(Young) - abs(Older)
    [tW, diffPS, cmapW] = clip_and_diff(t, abs(PS_y_aligned), abs(PS_o_aligned), xlims, plotTitle, cmapFull);

    % ----- Plot -----
    hold on; grid on; box on;

    plot_time_colored_segments(tW, diffPS(:,1), '-',  lw, cmapW); % PC1 solid
    plot_time_colored_segments(tW, diffPS(:,2), '--', lw, cmapW); % PC2 dashed
    plot_time_colored_segments(tW, diffPS(:,3), ':',  lw, cmapW); % PC3 dotted

    xlim(xlims);
    apply_ylims(ylims_user);

    xlabel('Time (s)');
    ylabel('Difference');
    title(sprintf('%s (cond %d)', plotTitle, cond_to_plot), 'Interpreter','none');

    % ----- Legend (bottom-right, large, Helvetica) -----
    h1 = plot(nan, nan, '-',  'LineWidth', lw, 'Color', [0 0 0]);
    h2 = plot(nan, nan, '--', 'LineWidth', lw, 'Color', [0 0 0]);
    h3 = plot(nan, nan, ':',  'LineWidth', lw, 'Color', [0 0 0]);
    lgd = legend([h1 h2 h3], {'PC #1', 'PC #2', 'PC #3'}, 'Location', 'southeast');
    set(lgd, 'FontName', 'Helvetica', 'FontSize', 14);

    drawnow;
end

function plot_local_vs_global_avg(RQA_y_local, RQA_o_local, RQA_y_global, RQA_o_global, ...
    cond_to_plot, xlims, ylims_user, lw)

    require_phasespace(RQA_y_local,  'RQA_y_local');
    require_phasespace(RQA_o_local,  'RQA_o_local');
    require_phasespace(RQA_y_global, 'RQA_y_global');
    require_phasespace(RQA_o_global, 'RQA_o_global');

    % Grab phase space coords for the chosen condition
    PS_yL = RQA_y_local.PhaseSpace.coords{cond_to_plot};
    PS_oL = RQA_o_local.PhaseSpace.coords{cond_to_plot};
    t_yL  = RQA_y_local.PhaseSpace.time(:);
    t_oL  = RQA_o_local.PhaseSpace.time(:);

    PS_yG = RQA_y_global.PhaseSpace.coords{cond_to_plot};
    PS_oG = RQA_o_global.PhaseSpace.coords{cond_to_plot};
    t_yG  = RQA_y_global.PhaseSpace.time(:);
    t_oG  = RQA_o_global.PhaseSpace.time(:);

    if size(PS_yL,2) < 3 || size(PS_oL,2) < 3 || size(PS_yG,2) < 3 || size(PS_oG,2) < 3
        error('Local/Global phase-space coordinates do not have 3 columns. Did you run with ''principalcomps'',[1:3]?');
    end

    % Build LOCAL average on LOCAL young timebase
    [tL, PS_oL_aligned, PS_yL_aligned] = align_to_young_time(t_yL, PS_yL, t_oL, PS_oL);
    local_avg = 0.5*(PS_yL_aligned + PS_oL_aligned);

    % Build GLOBAL average, interpolated onto LOCAL timebase (so subtraction is aligned)
    [tG, PS_oG_aligned, PS_yG_aligned] = align_to_young_time(t_yG, PS_yG, t_oG, PS_oG);
    global_avg_Gtime = 0.5*(PS_yG_aligned + PS_oG_aligned);

    global_avg = zeros(size(local_avg));
    global_avg(:,1) = interp1(tG, global_avg_Gtime(:,1), tL, 'linear', 'extrap');
    global_avg(:,2) = interp1(tG, global_avg_Gtime(:,2), tL, 'linear', 'extrap');
    global_avg(:,3) = interp1(tG, global_avg_Gtime(:,3), tL, 'linear', 'extrap');

    % ----- FULL colormap based on the LOCAL timebase length (matches video indexing) -----
    cmapFull = jet(numel(tL));

    % Clip + diff
    % Compare magnitudes by treating numbers as positive
    % i.e., diff = abs(Local_avg) - abs(Global_avg)
    plotTitle = 'LOCAL vs GLOBAL (avg(Y,O))';
    [tW, diffPS, cmapW] = clip_and_diff(tL, abs(local_avg), abs(global_avg), xlims, plotTitle, cmapFull); % Local - Global (magnitudes)

    % Plot
    hold on; grid on; box on;

    plot_time_colored_segments(tW, diffPS(:,1), '-',  lw, cmapW); % PC1 solid
    plot_time_colored_segments(tW, diffPS(:,2), '--', lw, cmapW); % PC2 dashed
    plot_time_colored_segments(tW, diffPS(:,3), ':',  lw, cmapW); % PC3 dotted

    xlim(xlims);
    apply_ylims(ylims_user);

    xlabel('Time (s)');
    ylabel('Difference');
    title(sprintf('%s (cond %d)', plotTitle, cond_to_plot), 'Interpreter','none');

    % ----- Legend (bottom-right, large, Helvetica) -----
    h1 = plot(nan, nan, '-',  'LineWidth', lw, 'Color', [0 0 0]);
    h2 = plot(nan, nan, '--', 'LineWidth', lw, 'Color', [0 0 0]);
    h3 = plot(nan, nan, ':',  'LineWidth', lw, 'Color', [0 0 0]);
    lgd = legend([h1 h2 h3], {'PC #1', 'PC #2', 'PC #3'}, 'Location', 'southeast');
    set(lgd, 'FontName', 'Helvetica', 'FontSize', 14);

    drawnow;
end

function require_phasespace(RQA, nameStr)
    if ~isfield(RQA,'PhaseSpace') || ~isfield(RQA.PhaseSpace,'coords') || ~isfield(RQA.PhaseSpace,'time')
        error('%s is missing PhaseSpace.coords or PhaseSpace.time. Run the modified BROADNESS_PhaseSpace_RQA first.', nameStr);
    end
end

function [t, PS_o_aligned, PS_y_aligned] = align_to_young_time(t_y, PS_y, t_o, PS_o)
    % Align older onto young time base (or keep if identical)
    same_time = (numel(t_y) == numel(t_o)) && all(abs(t_y - t_o) < 1e-12);
    if same_time
        t = t_y;
        PS_o_aligned = PS_o;
        PS_y_aligned = PS_y;
    else
        t = t_y;
        PS_y_aligned = PS_y;
        PS_o_aligned = zeros(size(PS_y,1), 3);
        PS_o_aligned(:,1) = interp1(t_o, PS_o(:,1), t, 'linear', 'extrap');
        PS_o_aligned(:,2) = interp1(t_o, PS_o(:,2), t, 'linear', 'extrap');
        PS_o_aligned(:,3) = interp1(t_o, PS_o(:,3), t, 'linear', 'extrap');
    end

    % Enforce same length safely
    n = min([numel(t), size(PS_y_aligned,1), size(PS_o_aligned,1)]);
    t = t(1:n);
    PS_y_aligned = PS_y_aligned(1:n,:);
    PS_o_aligned = PS_o_aligned(1:n,:);
end

function [tW, diffPS, cmapW] = clip_and_diff(t, A, B, xlims, contextLabel, cmapFull)
    % Clip to xlims, then compute A-B.
    % ALSO return cmapW as the subset of a FULL-length jet colormap,
    % so colors match the phase-space video exactly per timepoint index.

    if nargin < 6 || isempty(cmapFull)
        cmapFull = jet(numel(t));
    end

    if xlims(1) < min(t) || xlims(2) > max(t)
        warning('%s: requested x-limits [%.3f %.3f] exceed available data range [%.3f %.3f]. Data will only appear where available.', ...
            contextLabel, xlims(1), xlims(2), min(t), max(t));
    end
    idx = (t >= xlims(1)) & (t <= xlims(2));
    if ~any(idx)
        error('%s: no data points fall within xlims = [%.3f %.3f].', contextLabel, xlims(1), xlims(2));
    end

    tW = t(idx);
    diffPS = A(idx,:) - B(idx,:);
    cmapW = cmapFull(idx,:);
end

function apply_ylims(ylims_user)
    if ~isempty(ylims_user)
        if ~isnumeric(ylims_user) || numel(ylims_user) ~= 2 || any(~isfinite(ylims_user)) || ylims_user(1) >= ylims_user(2)
            error('ylims_user must be [] or a finite 1x2 vector with ylims_user(1) < ylims_user(2).');
        end
        ylim(ylims_user);
    end
end

function plot_time_colored_segments(x, y, lineStyle, lw, cmap)
    % Draw line segment-by-segment, color decided by time index (like phase space video).
    % If there is only a single sample in the clipped window, plotting segments would
    % draw nothing (because there are 0 segments). In that case, draw a visible point.

    if isempty(x) || isempty(y)
        return;
    end
    hold on;

    % Ensure cmap length matches number of points
    if size(cmap,1) < numel(x)
        cmap = jet(numel(x));
    end

    if numel(x) == 1
        % No segments possible -> draw a point so the figure isn't empty
        plot(x, y, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 18, ...
            'Color', cmap(1,:));
        return;
    end

    nSeg = numel(x) - 1;
    for i = 1:nSeg
        % Use the color of the starting point (consistent with per-timepoint coloring)
        col = cmap(i,:);
        plot(x(i:i+1), y(i:i+1), 'LineStyle', lineStyle, 'LineWidth', lw, 'Color', col);
    end
end


%% 3d) STATS - RQA ANOVA TEST

% Requires:
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

% Container for results (now with q-values)
Res = table('Size',[numel(metricNames) 7], ...
            'VariableTypes',{'string','double','double','double','double','double','double'}, ...
            'VariableNames',{'Metric','p_Group','p_Condition','p_Interaction', ...
                             'q_Group','q_Condition','q_Interaction'});

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

    % Mixed effects model: fixed Group, Condition, Interaction effects. Random intercept per subject
    lme = fitlme(Tlong, 'Value ~ Group*Condition + (1|Subject)');

    % F-tests for fixed effects (Satterthwaite df)
    A = anova(lme,'DFMethod','Satterthwaite');
    pG  = A.pValue(strcmp(A.Term,'Group'));
    pC  = A.pValue(strcmp(A.Term,'Condition'));
    pIx = A.pValue(strcmp(A.Term,'Group:Condition'));

    Res.Metric(m)        = string(met);
    Res.p_Group(m)       = pG;
    Res.p_Condition(m)   = pC;
    Res.p_Interaction(m) = pIx;

    fprintf('%-8s  Group p=%.4g | Condition p=%.4g | Interaction p=%.4g\n', ...
            met, pG, pC, pIx);
end

%% Joint FDR correction (BH) across all tests (Group/Condition/Interaction × metrics)

pG  = Res.p_Group;
pC  = Res.p_Condition;
pIx = Res.p_Interaction;

% Stack all p-values
pAll = [pG; pC; pIx];

% Benjamini–Hochberg (custom)
qAll = bh_fdr(pAll);

M = numel(metricNames);

Res.q_Group       = qAll(1:M);
Res.q_Condition   = qAll(M+1:2*M);
Res.q_Interaction = qAll(2*M+1:3*M);

disp(Res);

%% FDR correction helper

function q = bh_fdr(p)
    % Benjamini–Hochberg FDR correction (vectorized, supports NaNs)

    q = nan(size(p));
    [ps, idx] = sort(p(:), 'ascend', 'MissingPlacement','last');

    m = sum(~isnan(ps));
    ranks = (1:m)';

    adj = ps(1:m) .* m ./ ranks;
    adj = cummin(flipud(adj));
    adj = flipud(adj);
    adj(adj>1) = 1;

    q(idx(1:m)) = adj;
end

%% 3e) Plotting RQA metrics — ONE FIGURE (2x4 grid), VIOLINS per metric (using draw_violin_mean_SD)

plot_outdir = '/main_path/Output/RQA';
if ~exist(plot_outdir,'dir')
    mkdir(plot_outdir);
end

% ---------- USER-ADJUSTABLE PARAMETER ----------
% Proportion of the data range added as extra empty space above the highest
% data point.
EXTRA_TOP_PAD = 0.30;   % e.g. 0.2 (little extra), 0.3 (default), 0.5 (lots)

% Colors per combo (consistent with Piece 1)
% Order here matches conditions:
% 1: Young-Local, 2: Old-Local, 3: Young-Global, 4: Old-Global

comb_colors = [ ...
    0.1 0.2 0.5;  % Young-Local   (blue)
    0.4 0.6 0.85; % Young-Global  (light blue)
    0.5 0.1 0.1;  % Old-Local     (red)
    0.8 0.2 0.2]; % Old-Global    (light red)

% X-axis / legend order inside each panel:
condLabelsOrder = {'Young-Local','Old-Local','Young-Global','Old-Global'};

% -------------------------------------------------------------------------
% 1) Get basic info and metric names
% -------------------------------------------------------------------------
RQA_YL = RQA_BROADNESS_young_local.RQA_metrics;   % young, local (cell array)
RQA_OL = RQA_BROADNESS_older_local.RQA_metrics;   % old, local   (cell array)
RQA_YG = RQA_BROADNESS_young_global.RQA_metrics;  % young, global
RQA_OG = RQA_BROADNESS_older_global.RQA_metrics;  % old, global

nY = numel(RQA_YL);
nO = numel(RQA_OL);

% Sanity checks (same as before)
if numel(RQA_YG) ~= nY
    error('Young-Global and Young-Local have different number of subjects.');
end
if numel(RQA_OG) ~= nO
    error('Old-Global and Old-Local have different number of subjects.');
end

% Extract metric names from first non-empty table
AllCells = [RQA_YL(:); RQA_OL(:); RQA_YG(:); RQA_OG(:)];
idx1 = find(~cellfun(@isempty, AllCells), 1, 'first');
if isempty(idx1)
    error('No non-empty RQA tables found.');
end
metricNames = AllCells{idx1}.Properties.VariableNames;  % e.g. RR, L, DET, ...

% --- CUSTOM METRIC PLOTTING ORDER ---
desiredOrder = {'RR','L','TT','V_max','ENTR','DET','LAM','DIV'}; 

% Keep only metrics that actually exist, and append any leftovers at the end
desiredOrder = desiredOrder(ismember(desiredOrder, metricNames));
leftovers    = metricNames(~ismember(metricNames, desiredOrder));
metricNames  = [desiredOrder leftovers];

% -------------------------------------------------------------------------
% 2) Prepare figure: 2x4 tiles (for 8 metrics)
% -------------------------------------------------------------------------
nMetrics = numel(metricNames);  % should be 8

fig = figure('Color','w','Units','pixels','Position',[30 30 2600 1700]);
tl  = tiledlayout(fig, 2, 4, 'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% 3) Loop over metrics and create violins per metric using draw_violin_mean_SD
% -------------------------------------------------------------------------
for m = 1:nMetrics
    metName = metricNames{m};
    ax = nexttile(tl, m); 
    hold(ax,'on');

    % -------------------------------------------------------------
    % Collect values for this metric across subjects and conditions
    % -------------------------------------------------------------
    YL_vals = nan(nY,1);
    YG_vals = nan(nY,1);
    OL_vals = nan(nO,1);
    OG_vals = nan(nO,1);

    % Young
    for s = 1:nY
        Tloc  = RQA_YL{s};
        Tglob = RQA_YG{s};

        if isempty(Tloc) || ~istable(Tloc) || ...
           isempty(Tglob) || ~istable(Tglob) || ...
           ~ismember(metName, Tloc.Properties.VariableNames) || ...
           ~ismember(metName, Tglob.Properties.VariableNames)
            continue;
        end

        YL_vals(s) = Tloc{1, metName};
        YG_vals(s) = Tglob{1, metName};
    end

    % Old
    for s = 1:nO
        Tloc  = RQA_OL{s};
        Tglob = RQA_OG{s};

        if isempty(Tloc) || ~istable(Tloc) || ...
           isempty(Tglob) || ~istable(Tglob) || ...
           ~ismember(metName, Tloc.Properties.VariableNames) || ...
           ~ismember(metName, Tglob.Properties.VariableNames)
            continue;
        end

        OL_vals(s) = Tloc{1, metName};
        OG_vals(s) = Tglob{1, metName};
    end

    % Group into a cell array for easier looping (index meaning fixed):
    % 1: Young-Local, 2: Young-Global, 3: Old-Local, 4: Old-Global
    groupVals = {
        YL_vals(~isnan(YL_vals));  % 1: Young-Local
        YG_vals(~isnan(YG_vals));  % 2: Young-Global
        OL_vals(~isnan(OL_vals));  % 3: Old-Local
        OG_vals(~isnan(OG_vals));  % 4: Old-Global
    };

    % Check if there is any data at all
    if all(cellfun(@isempty, groupVals))
        title(ax, metName, 'Interpreter','none');
        ylabel(ax, 'Value');
        text(ax, 0.5, 0.5, 'No data', 'HorizontalAlignment','center');
        box(ax, 'on');
        continue;
    end

    % -------------------------------------------------------------
    % Draw violins using the same helper as in Piece 1
    % Order on x-axis: Young-Local, Old-Local, Young-Global, Old-Global
    % -------------------------------------------------------------
    xpositions = 1:4;
    plot_order = [1 3 2 4];  % indices into groupVals / colors

    for ip = 1:numel(plot_order)
        g = plot_order(ip);             % which condition to plot
        thisVals = groupVals{g};
        if isempty(thisVals)
            continue;
        end
        % draw_violin_mean_SD should be on the path (same as Piece 1)
        draw_violin_mean_SD(ax, thisVals, xpositions(ip), comb_colors(g,:));
    end

    % -------------------------------------------------------------
    % Get data-driven y-range (robustly), same logic as Piece 1
    % -------------------------------------------------------------
    kids = findall(ax);      % all children under this axis
    Y = [];

    for h = reshape(kids,1,[])
        if isprop(h,'YData')
            y = get(h,'YData');
            if ~isempty(y)
                y = y(:);
                if isnumeric(y)
                    Y = [Y; y]; %#ok<AGROW>
                end
            end
        end
    end

    finiteY = Y(isfinite(Y));

    if isempty(finiteY)
        ymin = 0;
        ymax = 1;
    else
        ymin = min(finiteY);
        ymax = max(finiteY);
    end

    if ymax == ymin
        ymax = ymin + 1;
    end

    yr = ymax - ymin;

    % ---------- Apply padding ----------
    pad_lower = 0.10 * yr;           % fixed bottom padding
    pad_upper = EXTRA_TOP_PAD * yr;  % user-controlled top padding

    ylim(ax, [ymin - pad_lower, ymax + pad_upper]);
    xlim(ax, [0.5 4.5]);

    % ---------- Cosmetics (axes) ----------
    set(ax, 'XTick', xpositions, 'XTickLabel', condLabelsOrder, ...
            'TickLabelInterpreter','none');
    xtickangle(ax, 30);
    ylabel(ax, metName, 'Interpreter','none');
    title(ax, metName, 'FontWeight','bold','Interpreter','none');
    grid(ax,'on'); 
    box(ax,'on');
end

% -------------------------------------------------------------------------
% Legend tile (same style as Piece 1)
% -------------------------------------------------------------------------
axL = nexttile(tl, 9);
cla(axL);
hold(axL,'on');
axis(axL,'off');

legend_items  = gobjects(0);
legend_labels = strings(0);
for k = 1:4
    legend_items(end+1) = plot(axL, NaN, NaN, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', comb_colors(k,:), 'MarkerEdgeColor', 'k'); %#ok<AGROW>
    legend_labels(end+1) = condLabelsOrder{k}; %#ok<AGROW>
end

lg = legend(axL, legend_items, legend_labels, 'Location', 'northwest');
set(lg, 'Interpreter','none', 'FontSize',13, 'Box','off');

% -------------------------------------------------------------------------
% Title, fonts, save specification
% -------------------------------------------------------------------------
sgtitle(tl, 'RQA Metrics: Violin plots per condition', 'FontWeight','bold');

% Set global font
set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica');

out_png = fullfile(plot_outdir, 'RQA_violin_like_2x4_metrics.png');
out_fig = fullfile(plot_outdir, 'RQA_violin_like_2x4_metrics.fig');

set(fig, 'PaperPositionMode','auto');
print(fig, out_png, '-dpng', '-r300');
saveas(fig, out_fig);

fprintf('RQA violin-like plot saved to:\n  %s\n  %s\n', out_png, out_fig);



%% 4a) SPATIAL GRADIENTS EMBEDDING AND CLUSTERING ANALYSIS

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply use the structure outputted by the BROADNESS_NetworkEstimation function
% Additional optional inputs can be provided, as described in the function. 
% Delete 'scatterplots','all' to only visualize optimal k in the cluster plot.

load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;

%%% ------------------ COMPUTATION --------------------- %%%

% Run 1 of these 5 lines in seperate

% Overall
%SPATIAL_GRADIENTS_BROADNESS = BROADNESS_SpatialGradients(BROADNESS,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords,'outpath', '/main_path/Output');

% Group/condition specific
%SPATIAL_GRADIENTS_BROADNESS_young_local = BROADNESS_SpatialGradients(BROADNESS_young_local,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);
%SPATIAL_GRADIENTS_BROADNESS_older_local = BROADNESS_SpatialGradients(BROADNESS_older_local,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);
%SPATIAL_GRADIENTS_BROADNESS_young_global = BROADNESS_SpatialGradients(BROADNESS_young_global,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);
%SPATIAL_GRADIENTS_BROADNESS_older_global = BROADNESS_SpatialGradients(BROADNESS_older_global,'principalcomps',[1:3],'evalclusters',1,'mni_coords', Options.MNI_coords);

%% 4b) Finding the optimal cluster solution with 5000 iterations for k = 2:40

load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;

S_GRAD_Final = BROADNESS_SpatialGradients(BROADNESS,'principalcomps',[1:3],'evalclusters', 5000, 'nclusters', 2:40, 'mni_coords', Options.MNI_coords,'outpath', '/main_path/Output','scatterplots', 'all');

%% 4c) Plot Sillouhete plot of optimal k frequency count for all iterations

freqTbl = S_GRAD_Final.optimalK_Frequency;

figure;
b = bar(freqTbl.K, freqTbl.Count, 'FaceColor', [0 0.45 0]);
xlabel('Optimal number of clusters (K)');
ylabel('Frequency across evalclusters runs');
title('Distribution of silhouette-optimal K');
grid on;

%% 4d) Plot 3 stacked jittered data points for PC1/PC2/PC3 values across all clusters
% Reads CSVs: OptimalK_14_Cluster_XX_PC_points.csv
% Folder: /main_path/Output/BROADNESS_Output/ClusterPCcoords

clear; close all; clc;

dataDir   = '/main_path/Output/BROADNESS_Output/ClusterPCcoords';
optimalK  = 16; % Inserted manually after obtaining results from 4b and 4c

% Limits
xLimits = [-0.1 0.1];

% Visual settings
markerSize  = 18;    % dot size
jitterAmp   = 0.08;  % vertical jitter within each box (0 = no jitter)
boxHeight   = 0.55;  % height of each box in y-units (visual only)
lineCenters = [3, 2, 1];  % top-to-bottom: PC1, PC2, PC3 (stacked)

% Color map must match scatter: jet(k)*0.9
cmap = jet(optimalK) * 0.9;

% Which PCs to plot (must match the CSV column names)
pcNames = {'PC1','PC2','PC3'};

% ---------------- LOAD DATA ----------------
clusterVals = cell(optimalK, 3);  % {cluster, pcIndex} -> vector of values

for cl = 1:optimalK
    fn = fullfile(dataDir, sprintf('OptimalK_%d_Cluster_%02d_PC_points.csv', optimalK, cl));
    if ~exist(fn,'file')
        error('Missing file: %s', fn);
    end

    T = readtable(fn);

    % Validate required columns exist
    for p = 1:3
        if ~ismember(pcNames{p}, T.Properties.VariableNames)
            error('File %s does not contain column "%s". Columns are: %s', ...
                  fn, pcNames{p}, strjoin(T.Properties.VariableNames, ', '));
        end
    end

    clusterVals{cl,1} = T.(pcNames{1});
    clusterVals{cl,2} = T.(pcNames{2});
    clusterVals{cl,3} = T.(pcNames{3});
end

%% One image per cluster: 3 stacked 1D point-rugs
clear; close all; clc;

dataDir  = '/main_path/Output/BROADNESS_Output/ClusterPCcoords';
optimalK = 16; % Inserted after running 4b and 4c to identify optimal number of clusters

% Match limits
lims = [-0.1 0.1];

% Match cluster colors
cmap = jet(optimalK) * 0.9;

% Output folder
outDir = fullfile(dataDir, 'Cluster1DPointRugs');
if ~exist(outDir,'dir'); mkdir(outDir); end

% -------------------- USER-EDITABLE LAYOUT VARIABLES --------------------
lineSep   = 0.0060;             % baseline-to-baseline spacing
laneY     = [2 1 0] * lineSep; % BN1 top, BN2 middle, BN3 bottom

leftLabels = {'BN1','BN2','BN3'};
labelX     = lims(1) - 0.0040;
labelFS    = 12;

yPadTop    = 0.020;
yPadBottom = 0.0040;

vTickHalfHeight = 0.0018;   % <<< half-height of vertical ticks on bars
% -----------------------------------------------------------------------

% Plot settings
yJitter = 0.0017;
ptSize  = 42;

xTickPositions = [lims(1), 0, lims(2)];

for cl = 1:optimalK
    fn = fullfile(dataDir, sprintf('OptimalK_%d_Cluster_%02d_PC_points.csv', optimalK, cl));
    if ~exist(fn,'file')
        warning('Missing file: %s (skipping)', fn);
        continue;
    end

    T = readtable(fn);

    needed = {'PC1','PC2','PC3'};
    if ~all(ismember(needed, T.Properties.VariableNames))
        warning('File %s missing PC columns.', fn);
        continue;
    end

    X = T.PC1(:);
    Y = T.PC2(:);
    Z = T.PC3(:);

    keep = (X>=lims(1) & X<=lims(2)) & ...
           (Y>=lims(1) & Y<=lims(2)) & ...
           (Z>=lims(1) & Z<=lims(2));

    Xp = X(keep);
    Yp = Y(keep);
    Zp = Z(keep);

    % --- Figure ---
    fig = figure('Color','w', 'Visible','off');
    ax  = axes(fig); hold(ax,'on');

    % Baseline lanes
    for k = 1:3
        plot(ax, lims, [laneY(k) laneY(k)], 'k-', 'LineWidth', 1.0);
    end

    % Vertical ticks on each baseline
    for k = 1:3
        for xt = xTickPositions
            plot(ax, [xt xt], ...
                laneY(k) + [-vTickHalfHeight vTickHalfHeight], ...
                'k-', 'LineWidth', 1.0);
        end
    end

    % Left labels (Helvetica, not bold)
    for k = 1:3
        text(ax, labelX, laneY(k), leftLabels{k}, ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','middle', ...
            'FontName','Helvetica', ...
            'FontSize', labelFS, ...
            'FontWeight','normal');
    end

    % Point rugs
    if ~isempty(Xp)
        scatter(ax, Xp, laneY(1) + (rand(size(Xp))-0.5)*2*yJitter, ...
            ptSize, 'MarkerFaceColor', cmap(cl,:), 'MarkerEdgeColor','none');
    end
    if ~isempty(Yp)
        scatter(ax, Yp, laneY(2) + (rand(size(Yp))-0.5)*2*yJitter, ...
            ptSize, 'MarkerFaceColor', cmap(cl,:), 'MarkerEdgeColor','none');
    end
    if ~isempty(Zp)
        scatter(ax, Zp, laneY(3) + (rand(size(Zp))-0.5)*2*yJitter, ...
            ptSize, 'MarkerFaceColor', cmap(cl,:), 'MarkerEdgeColor','none');
    end

    % Axes formatting
    xlim(ax, lims);
    ylim(ax, [min(laneY)-yPadBottom, max(laneY)+yPadTop]);

    xticks(ax, xTickPositions);
    box(ax,'off'); grid(ax,'off');

    % Hide x-axis line but keep numbers
    set(ax, 'XColor', 'none');     % hides axis line & default ticks
    set(ax, 'XTickLabelMode','auto');  % keep numeric labels

    set(ax,'YColor','none','FontName','Helvetica');

    title(ax, sprintf('Cluster %02d (OptimalK=%d, n=%d)', ...
        cl, optimalK, numel(Xp)), 'FontWeight','bold');

    % Save
    outPng = fullfile(outDir, ...
        sprintf('OptimalK_%d_Cluster_%02d_point_rugs.png', optimalK, cl));
    exportgraphics(fig, outPng, 'Resolution', 300);
    close(fig);

    fprintf('Saved: %s\n', outPng);
end

%% 4e) Export MNI coordinates of all clusters

% Drop-in block (not a function). Requires MATLAB's niftiinfo/niftiread.
% Exports: Index, X, Y, Z, Activation
%
% NOTES:
% - Uses the NIfTI affine: info.Transform.T
% - Activation is the voxel value in the NIfTI (for your current pipeline, likely 0/1)
% - If Transform is wrong/missing, coordinates will be wrong.

% ---- USER SETTINGS (edit these) ----
optimalK = 16; % Inserted after running 4b and 4c to identify optimal number of clusters
niiDir  = '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/Clusters_k=16_Final_With_Niftis/' % where SpatialGradients_*.nii live
outDir  = '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/Clusters_k=16_Final_With_Niftis/' % where tables will be saved
pattern = sprintf('SpatialGradients_OptimalK_%d_Cluster_*.nii*', optimalK);

% ---- SAFETY CHECKS ----
if ~(exist('niftiinfo','file')==2 && exist('niftiread','file')==2)
    error('This block requires niftiinfo/niftiread (MATLAB built-in). Update MATLAB or use your NIfTI toolbox + affine parsing.');
end
if ~exist(niiDir,'dir')
    error('NIfTI folder not found: %s', niiDir);
end
if ~exist(outDir,'dir'); mkdir(outDir); end

files = dir(fullfile(niiDir, pattern));
if isempty(files)
    warning('No files found for pattern: %s', fullfile(niiDir, pattern));
else
    fprintf('Found %d cluster NIfTIs. Exporting tables to: %s\n', numel(files), outDir);
end

for f = 1:numel(files)
    fn = fullfile(files(f).folder, files(f).name);

    % Read header + volume
    info = niftiinfo(fn);
    vol  = double(niftiread(info));

    % Handle potential 4th dimension of size 1
    if ndims(vol) == 4 && size(vol,4) == 1
        vol = vol(:,:,:,1);
    end

    % Nonzero voxels = "active"
    linIdx = find(vol ~= 0);

    if isempty(linIdx)
        warning('No nonzero voxels in: %s', files(f).name);
        continue;
    end

    % Convert linear indices -> voxel subscripts (i,j,k)
    [i,j,k] = ind2sub(size(vol), linIdx);
    act = vol(linIdx);

    % Convert voxel subscripts -> MNI/world using affine
    % MATLAB uses 1-based subscripts; Transform.T is consistent for niftiinfo.
    T = info.Transform.T;                 % 4x4
    vox = [i(:) j(:) k(:) ones(numel(i),1)];
    xyz = vox * T;                        % Nx4
    X = xyz(:,1); Y = xyz(:,2); Z = xyz(:,3);

    % Build table
    Index = (1:numel(act))';
    TBL = table(Index, X, Y, Z, act, ...
        'VariableNames', {'Index','X','Y','Z','Activation'});

    % Write outputs
    [~,base,~] = fileparts(files(f).name);
    outCsv  = fullfile(outDir, [base '_MNIcoords.csv']);
    writetable(TBL, outCsv);

    % Try XLSX too (may fail on some Linux setups)
    outXlsx = fullfile(outDir, [base '_MNIcoords.xlsx']);
    try
        writetable(TBL, outXlsx);
    catch
        % ignore; CSV is the reliable output
    end

    fprintf('Exported %s (%d voxels)\n', outCsv, height(TBL));

    % ---- Optional sanity checks ----
    % For an 8mm MNI grid ~ multiples of 2 or 8 depending on template.
    % If these look wild, the header transform is likely not correct.
    % disp(TBL(1:min(5,height(TBL)),:));
end

%% 4f) Converts cluster MNI coords to AAL labels: Outputs per labels cluster sheet in one Excel workbook.
%
% Takes all cluster coordinate tables exported in Code 1 (*_MNIcoords.csv or .xlsx),
% assigns AAL (or other integer-coded atlas) labels voxelwise, and writes ONE
% Excel file with ONE sheet per cluster containing:
% Area | Dist0 | Dist1 | ... | DistR | Total
%
% where DistR = number of voxels whose nearest-label snapping used radius R
% (Chebyshev ring radius in ATLAS VOXELS; r=0 = exact atlas hit).

addpath('/home/mathiasha/MATLAB_Add-Ons/Collections/spm_25.01.02/spm/')
clear; clc;

% ----------------------- USER SETTINGS -----------------------
optimalK = 16;

% Folder containing the exported MNI tables from Code 1
% (e.g., SpatialGradients_OptimalK_16_Cluster_1_MNIcoords.csv)
clusterTableDir = '/main_path/Output/BROADNESS_Output/BROADNESS_nifti/Clusters_k=16_Final_With_Niftis/';

% Atlas in MNI space + label lookup table
atlasNiftiPath  = '/main_path/AAL3/AAL3v1.nii';
atlasLabelsPath = '/main_path/AAL3/AAL3v1.nii.txt';

unknownLabelName = "Unknown/Background";

% Max snapping distance in ATLAS VOXELS (0 = exact only)
maxSearchRadius = 3;

% Output Excel (one workbook; one sheet per cluster)
outXlsx = fullfile(clusterTableDir, sprintf('Clusters_k=%d_AAL_labels_by_cluster.xlsx', optimalK));

% ---- Dependencies / safety checks ----
if exist('spm_vol','file')~=2 || exist('spm_read_vols','file')~=2
    error('This script requires SPM on the MATLAB path (spm_vol, spm_read_vols).');
end
assert(exist(clusterTableDir,'dir')==7, 'Cluster table folder not found: %s', clusterTableDir);
assert(exist(atlasNiftiPath,'file')==2, 'Atlas NIfTI not found: %s', atlasNiftiPath);
assert(exist(atlasLabelsPath,'file')==2, 'Atlas labels file not found: %s', atlasLabelsPath);

% ---- Load atlas volume + LUT ----
V = spm_vol(atlasNiftiPath);
A = spm_read_vols(V);
lut = loadAtlasLUT(atlasLabelsPath, unknownLabelName);

% ---- Find cluster coordinate tables ----
% Prefer CSV if both exist, otherwise use XLSX.
csvFiles  = dir(fullfile(clusterTableDir, sprintf('SpatialGradients_OptimalK_%d_Cluster_*_MNIcoords.csv',  optimalK)));
xlsxFiles = dir(fullfile(clusterTableDir, sprintf('SpatialGradients_OptimalK_%d_Cluster_*_MNIcoords.xlsx', optimalK)));

if isempty(csvFiles) && isempty(xlsxFiles)
    error('No cluster MNI tables found in: %s', clusterTableDir);
end

% Build map: clusterNumber -> filepath, preferring CSV
clusterMap = containers.Map('KeyType','double','ValueType','char');

for f = 1:numel(xlsxFiles)
    cnum = parseClusterNumber(xlsxFiles(f).name);
    if ~isnan(cnum)
        clusterMap(cnum) = fullfile(xlsxFiles(f).folder, xlsxFiles(f).name);
    end
end
for f = 1:numel(csvFiles)
    cnum = parseClusterNumber(csvFiles(f).name);
    if ~isnan(cnum)
        % CSV overwrites XLSX preference
        clusterMap(cnum) = fullfile(csvFiles(f).folder, csvFiles(f).name);
    end
end

% Sort cluster numbers
clusterNums = sort(cell2mat(keys(clusterMap)));

fprintf('Found %d cluster tables.\n', numel(clusterNums));
fprintf('Writing output workbook:\n  %s\n', outXlsx);

% ---- Process each cluster and write one sheet per cluster ----
for iC = 1:numel(clusterNums)
    cnum = clusterNums(iC);
    inPath = clusterMap(cnum);

    % Read table
    T = readtable(inPath);

    % Find MNI columns (X,Y,Z) robustly
    [xCol, yCol, zCol] = guessMNIColumns(T);

    XYZ = [T.(xCol), T.(yCol), T.(zCol)];
    if ~isnumeric(XYZ) || size(XYZ,2)~=3
        error('Cluster %d: MNI columns are not numeric 3D coords. Detected: %s %s %s', cnum, xCol, yCol, zCol);
    end

    % Map MNI -> atlas label with snapping  (FIXED: function returns 3 outputs)
    [labels, usedRadius, labelIDs] = mniToAtlasLabelWithSnap2(XYZ, V, A, lut, unknownLabelName, maxSearchRadius);

    % Count labels by distance radius used (0..R)
    countsByR = countLabelsByRadius(labels, usedRadius, maxSearchRadius);

    % Convert to the exact kind of table Code 2 produced (Area + Dist0..DistR + Total)
    outTable = countsByRadiusToTable(countsByR, maxSearchRadius, "Area");

    % Optional: raw per-voxel mapping too
    voxelTable = table( ...
        (1:size(XYZ,1))', XYZ(:,1), XYZ(:,2), XYZ(:,3), labelIDs(:), string(labels(:)), usedRadius(:), ...
        'VariableNames', {'Index','X','Y','Z','AtlasID','Area','Distance'} );

    % Excel sheet name
    sheetName = sprintf('Cluster_%d', cnum);

    % Write the summary counts table to the cluster sheet
    writetable(outTable, outXlsx, 'Sheet', sheetName, 'WriteMode', 'overwritesheet');

    % Also write the per-voxel assignments below the summary (same sheet),
    % leaving a blank row between (Excel-friendly).
    try
        % Write voxelTable starting at (row = height(outTable)+3, col = 1)
        startRow = height(outTable) + 3;
        writetable(voxelTable, outXlsx, 'Sheet', sheetName, 'Range', sprintf('A%d', startRow));
    catch ME
        % If MATLAB/Excel writer on Linux is finicky, at least keep the summary.
        warning('Cluster %d: Could not append per-voxel table; kept summary only. (%s) %s', cnum, sheetName, ME.message);
    end

    % Print quick diagnostics
    fracSnapped = mean(usedRadius(:) > 0);
    fprintf('Cluster %d: %d voxels | %.1f%% snapped (radius>0) | top area: %s\n', ...
        cnum, size(XYZ,1), 100*fracSnapped, string(outTable.Area(1)));
end

fprintf('\nDone.\nWorkbook saved:\n  %s\n', outXlsx);

% ========================= LOCAL FUNCTIONS =========================

function cnum = parseClusterNumber(filename)
% Extract cluster number from names like:
%   SpatialGradients_OptimalK_16_Cluster_3_MNIcoords.csv
    cnum = NaN;
    tok = regexp(filename, 'Cluster_(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        cnum = str2double(tok{1});
    end
end

function [labels, usedRadius, labelIDs] = mniToAtlasLabelWithSnap2(XYZmm, V, A, lut, unknownLabelName, maxR)
% Returns:
%  labels: assigned region name for each coordinate (string)
%  usedRadius:
%    0..maxR = snapped using that voxel shell radius (0 means exact voxel non-zero label)
%    -1      = could not find any non-zero label within maxR or out-of-bounds
%  labelIDs:
%    integer atlas parcel ID used (0 = unknown/background)
%
% NOTE: "radius" is Chebyshev shells in ATLAS VOXELS, not mm.

    n = size(XYZmm,1);
    labels     = strings(n,1);
    usedRadius = -1 * ones(n,1);
    labelIDs   = zeros(n,1);

    M  = V.mat;       %#ok<NASGU> % voxel->mm (kept for clarity)
    iM = inv(V.mat);  % mm->voxel

    volSize = size(A);

    for iPt = 1:n
        mm  = [XYZmm(iPt,:), 1]';
        vox = iM * mm;
        ijk0 = round(vox(1:3))';

        % Out of bounds -> Unknown
        if any(isnan(ijk0)) || any(ijk0 < 1) || ...
                ijk0(1) > volSize(1) || ijk0(2) > volSize(2) || ijk0(3) > volSize(3)
            labels(iPt)     = unknownLabelName;
            usedRadius(iPt) = -1;
            labelIDs(iPt)   = 0;
            continue;
        end

        val0 = A(ijk0(1), ijk0(2), ijk0(3));
        if isnan(val0); val0 = 0; end
        key0 = double(round(val0));

        % Exact hit must be non-zero AND in LUT
        if key0 ~= 0 && isKey(lut, key0)
            labels(iPt)     = string(lut(key0));
            usedRadius(iPt) = 0;
            labelIDs(iPt)   = key0;
            continue;
        end

        % Otherwise, search shells rad=1..maxR for nearest non-zero labels
        found = false;
        for rad = 1:maxR
            keysInShell = collectNonzeroKeysInShell(A, ijk0, rad);

            if ~isempty(keysInShell)
                % Choose the most frequent label key in this shell (mode; tiebreak smallest)
                chosenKey = modeTiebreakSmallest(keysInShell);

                if isKey(lut, chosenKey)
                    labels(iPt)   = string(lut(chosenKey));
                    labelIDs(iPt) = chosenKey;
                else
                    % LUT missing this key (unexpected) -> mark unknown but keep ID for debugging
                    labels(iPt)   = unknownLabelName;
                    labelIDs(iPt) = chosenKey;
                end

                usedRadius(iPt) = rad;
                found = true;
                break;
            end
        end

        if ~found
            labels(iPt)     = unknownLabelName;
            usedRadius(iPt) = -1;
            labelIDs(iPt)   = 0;
        end
    end
end

function col = pickFirstMatching(vars, candidates)
    col = "";
    lowerVars = lower(vars);
    for c = candidates
        idx = find(lowerVars == lower(string(c)), 1, 'first');
        if ~isempty(idx)
            col = vars(idx);
            return;
        end
    end
end

function [id, rUsed] = lookupWithSnap(A, ijk, dims, maxR)
% Returns atlas id at ijk if valid and nonzero, else nearest nonzero within Chebyshev radius.
% Chebyshev shell radius r: max(|dx|,|dy|,|dz|)==r

    id = 0;
    rUsed = 0;

    if isInside(ijk, dims)
        v = A(ijk(1), ijk(2), ijk(3));
        if ~isnan(v) && v ~= 0
            id = round(v);
            rUsed = 0;
            return;
        end
    end

    % Snap search
    for r = 1:maxR
        bestId = 0;

        % Iterate only over the shell (ring) at radius r
        for dx = -r:r
            for dy = -r:r
                for dz = -r:r
                    if max(abs([dx,dy,dz])) ~= r
                        continue; % not on the shell
                    end
                    p = ijk + [dx dy dz];
                    if ~isInside(p, dims)
                        continue;
                    end
                    v = A(p(1), p(2), p(3));
                    if ~isnan(v) && v ~= 0
                        bestId = round(v);
                        break;
                    end
                end
                if bestId ~= 0, break; end
            end
            if bestId ~= 0, break; end
        end

        if bestId ~= 0
            id = bestId;
            rUsed = r;
            return;
        end
    end

    % if still 0 -> unknown
    id = 0;
    rUsed = maxR + 1; % indicates "not found within maxR"
end

function tf = isInside(ijk, dims)
    tf = all(ijk >= 1) && ijk(1) <= dims(1) && ijk(2) <= dims(2) && ijk(3) <= dims(3);
end


%% 5a) BRAIN NETWORK MODULARITY CALCULATION - PCA APPLIED TO INDIVIDUAL DATA

% RELOAD THINGS TO AVOID CONFUSION
clear 
close all
clc

% Setup directories
path_home = '/main_path/BROADNESS_MEG_AuditoryRecognition-main/BROADNESS_Toolbox';
addpath(path_home)
BROADNESS_Startup(path_home);
addpath(fullfile(matlabroot,'toolbox','stats','stats'),'-begin') % makes sure the pca function is the standard function in matlab

%%% For each participant we calculate how many networks (PCs) it takes to
%%% to reach a certain threshold of explained varience in the brain and use
%%% this as a measure of brain network modularity/partitioning.

load('/main_path/Data/Mathias_MMN/MMNSubtracted_Average_SignFixed.mat');
DATA = dum;

% Remove first participant and define time vector

% 100 ms baseline
time = -0.100:0.004:0.8;
DATA = DATA(:,101:326,:,2:78);

DATA_global = DATA(:,:,1,:);
DATA_local = DATA(:,:,2,:);

size(DATA)

% Load groups and adjust subid after removing first participant
load('/main_path/Data/Mathias_MMN/groups.mat'); % older, young
older = older - 1;
young = young - 1;

older_subj = older(:);   % row-wise versions
young_subj = young(:);

%%% ------------------ BROADNESS --------------------- %%%

% Run BROADNESS network estimation (default parameters)
numSubjects = size(DATA_global, 4);

BROADNESS_global = cell(numSubjects,1);
BROADNESS_local = cell(numSubjects,1);

for s = 1:numSubjects
    % Extract subject data (keeping dimensionality intact)
    sub_global = DATA_global(:,:,:,s);
    sub_local  = DATA_local(:,:,:,s);

    % Run BROADNESS for this subject
    BROADNESS_global{s} = BROADNESS_NetworkEstimation(sub_global, time);
    BROADNESS_local{s}  = BROADNESS_NetworkEstimation(sub_local, time);

    fprintf('Subject %d/%d completed\n', s, numSubjects);
end


%% 5b) Calculate effective dimensionality (ED) on PCA data applied to individuals and plot per participant

%  This function estimates effective dimensionality (ED) and quadratic
%  Rényi entropy (H2) across frequencies, given the eigenspectrum of a
%  covariance matrix. It uses the entropy-based index derived by Pirk et al.
%  (2012) to estimate the effective number of uncorrelated measurements, as
%  described in Del Giudice (2020).
%
%  It can be applied to a 1-D vector of eigenvalues, as well as directly to
%  the FREQ.evals output produced by FREQNESS_NetworkEstimation().
%  If FREQ.evals is given as an input, this can consist of either a 2D or
%  3D matrix, depending on whether FREQNESS was run on individual
%  participants or at the group level.
%
%  When FREQ.evals is provided as input, the ED output will reflect the
%  effective number of uncorrelated measurements across the frequency
%  spectrum. If FREQ.evals contains participants as a 3rd dimension, then ED
%  will also include an additional dimension for multiple participants.
%  This will allow to visualize the grand-average entropy landscape and
%  eventually carry out statistical testing.

% Assumes one already has:
%   BROADNESS_global : nSubs x 1 cell
%   BROADNESS_local  : nSubs x 1 cell
%
% Each cell contains a struct with field:
%   .Variance_BrainNetworks : 225 x 1 double

% ---- Settings ----
fieldName = 'Variance_BrainNetworks';

% ---- Compute ED vectors ----
ED_global = computeED_fromBROADNESS(BROADNESS_global, fieldName);
ED_local  = computeED_fromBROADNESS(BROADNESS_local,  fieldName);

% ---- Quick sanity checks / summary ----
fprintf('GLOBAL: computed ED for %d/%d participants.\n', sum(~isnan(ED_global)), numel(ED_global));
fprintf('LOCAL : computed ED for %d/%d participants.\n', sum(~isnan(ED_local)),  numel(ED_local));

fprintf('GLOBAL ED: mean = %.3f, sd = %.3f\n', mean(ED_global,'omitnan'), std(ED_global,'omitnan'));
fprintf('LOCAL  ED: mean = %.3f, sd = %.3f\n', mean(ED_local,'omitnan'),  std(ED_local,'omitnan'));

% Optional: paired difference
ED_diff = ED_global - ED_local;
fprintf('DIFF (global-local): mean = %.3f, sd = %.3f\n', mean(ED_diff,'omitnan'), std(ED_diff,'omitnan'));

% ===== User-adjustable legend parameters =====
legendFontSize  = 11;        % controls legend text size
legendTokenSize = [24 12];   % controls line/marker size in legend [length height]
% ============================================

figure; 

plot(ED_global, 'o-', ...
    'Color', [0.75 0.6 0.9], ...   % light purple
    'LineWidth', 1.5); 
hold on;

plot(ED_local, 'o-', ...
    'Color', [0.4 0.1 0.6], ...    % dark purple
    'LineWidth', 1.5);

xlabel('Participant'); 
ylabel('Effective Dimensionality (ED)');
title('ED per participant');

grid on;

lgd = legend({'Global','Local'});
lgd.Location   = 'south';
lgd.Box        = 'off';
lgd.FontSize   = legendFontSize;
lgd.ItemTokenSize = legendTokenSize;

%% -------- Local function(s) --------
function ED = computeED_fromBROADNESS(BROADNESS_cell, fieldName)
    % Returns ED as [nSubs x 1] double. Participants that fail checks -> NaN.

    nSubs = numel(BROADNESS_cell);
    ED    = nan(nSubs, 1);

    for s = 1:nSubs
        x = BROADNESS_cell{s};

        % Basic structure checks
        if isempty(x) || ~isstruct(x) || ~isfield(x, fieldName)
            warning('Sub %d: missing struct or field "%s". Setting ED=NaN.', s, fieldName);
            continue;
        end

        lam = x.(fieldName);

        % Ensure numeric column vector
        if isempty(lam) || ~isnumeric(lam)
            warning('Sub %d: "%s" is empty/non-numeric. Setting ED=NaN.', s, fieldName);
            continue;
        end
        lam = lam(:); % force column

        % Clean invalid values
        lam(~isfinite(lam)) = 0;

        % ED formula expects nonnegative spectrum; clamp tiny negatives if present
        % (If seeing large negatives, something upstream is wrong.)
        if any(lam < -1e-12)
            warning('Sub %d: large negative values in "%s". Check upstream. Clamping negatives to 0.', s, fieldName);
        end
        lam(lam < 0) = 0;

        % Compute ED
        num   = (sum(lam))^2;
        denom = sum(lam.^2);

        % Guard against division by zero
        denom = max(denom, realmin);

        ED(s) = num / denom;
    end
end


%% 5c) Plotting ED for each YL, OL, YG, OG - VIOLIN PLOT

% Requires variables in workspace:
%   ED_local   [nSubs x 1]
%   ED_global  [nSubs x 1]
%   young      [1 x 37] participant indices (1..nSubs)
%   older      [1 x 40] participant indices (1..nSubs)
%
% Also requires the helper function on path:
%   draw_violin_mean_SD(ax, values, xposition, color)

% ---------------- User-adjustable ----------------
EXTRA_TOP_PAD = 0.30;   % extra headroom above max datapoint (for markers/text)
plot_outdir   = '';     % set to a folder to save; leave '' to not save

% Colors (same as RQA code)
% 1: Young-Local, 2: Old-Local, 3: Young-Global, 4: Old-Global
comb_colors = [ ...
    0.1 0.2 0.7;  ... Young-Local, blue
    0.4 0.6 0.85; ... Young-Global, light blue
    0.5 0.1 0.1;  ... Old-Local, red
    0.8 0.2 0.2]; ... Old-Global, light red

condLabelsOrder = {'Younger - Local','Older - Local','Younger - Global','Older - Global'};

% ---------------- Collect data ----------------
young = young(:);
older = older(:);

% Sanity checks: if these fail, the "young/older" are not indices into ED vectors
if any(young < 1) || any(older < 1) || any(young > numel(ED_local)) || any(older > numel(ED_local))
    error('Entries in "young" or "older" are out of bounds for ED vectors. Check participant numbering (IDs vs indices).');
end

YL_vals = ED_local(young);
OL_vals = ED_local(older);
YG_vals = ED_global(young);
OG_vals = ED_global(older);

% Group into a cell array for easier looping (index meaning fixed):
% 1: Young-Local, 2: Young-Global, 3: Old-Local, 4: Old-Global
groupVals = {
    YL_vals(~isnan(YL_vals) & isfinite(YL_vals));
    YG_vals(~isnan(YG_vals) & isfinite(YG_vals));
    OL_vals(~isnan(OL_vals) & isfinite(OL_vals));
    OG_vals(~isnan(OG_vals) & isfinite(OG_vals));
};

if all(cellfun(@isempty, groupVals))
    error('All ED groups are empty after removing NaN/Inf. Fix upstream ED computation or indexing.');
end

% ---------------- Plot (violin-style) ----------------
fig = figure('Color','w','Units','pixels','Position',[100 100 1400 700]);
tl  = tiledlayout(fig, 1, 2, 'TileSpacing','compact','Padding','compact');

ax = nexttile(tl, 1);
hold(ax,'on');

xpositions = 1:4;

% We want x-axis order: Young-Local, Old-Local, Young-Global, Old-Global
% But groupVals order is:   1=YL, 2=YG, 3=OL, 4=OG
plot_order = [1 3 2 4];  % indices into groupVals / colors
for ip = 1:numel(plot_order)
    g = plot_order(ip);
    thisVals = groupVals{g};
    if isempty(thisVals)
        continue;
    end
    draw_violin_mean_SD(ax, thisVals, xpositions(ip), comb_colors(g,:));
end

% ---------------- Data-driven y-limits (same logic as the RQA code) ----------------
kids = findall(ax);
Y = [];
for h = reshape(kids,1,[])
    if isprop(h,'YData')
        y = get(h,'YData');
        if ~isempty(y) && isnumeric(y)
            Y = [Y; y(:)]; %#ok<AGROW>
        end
    end
end

finiteY = Y(isfinite(Y));
if isempty(finiteY)
    ymin = 0; ymax = 1;
else
    ymin = min(finiteY);
    ymax = max(finiteY);
end
if ymax == ymin
    ymax = ymin + 1;
end

yr = ymax - ymin;
pad_lower = 0.10 * yr;
pad_upper = EXTRA_TOP_PAD * yr;

ylim(ax, [ymin - pad_lower, ymax + pad_upper]);
xlim(ax, [0.5 4.5]);

% ---------------- Cosmetics ----------------
set(ax, 'XTick', xpositions, 'XTickLabel', condLabelsOrder, ...
        'TickLabelInterpreter','none');
xtickangle(ax, 30);
ylabel(ax, 'Effective Dimensionality (ED)', 'Interpreter','none');
title(ax, 'ED by Group and Condition', 'FontWeight','bold', 'Interpreter','none');
grid(ax,'on');
box(ax,'on');

% ---------------- Legend tile (matching the RQA style) ----------------
axL = nexttile(tl, 2);
cla(axL);
hold(axL,'on');
axis(axL,'off');

legend_items  = gobjects(0);
legend_labels = strings(0);

% Legend order should match the x-axis labels (not plot_order indices)
% x-axis: YL, OL, YG, OG -> color rows: 1,2,3,4 respectively
for k = 1:4
    legend_items(end+1) = plot(axL, NaN, NaN, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', comb_colors(k,:), 'MarkerEdgeColor', 'k'); %#ok<AGROW>
    legend_labels(end+1) = condLabelsOrder{k}; %#ok<AGROW>
end

lg = legend(axL, legend_items, legend_labels, 'Location', 'northwest');
set(lg, 'Interpreter','none', 'FontSize',13, 'Box','off');

% Global font
set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica');

% Optional: overall title
sgtitle(tl, 'Effective Dimensionality (ED): Violin plots per condition', 'FontWeight','bold');

% ---------------- Optional save ----------------
if ~isempty(plot_outdir)
    if ~exist(plot_outdir,'dir'); mkdir(plot_outdir); end
    out_png = fullfile(plot_outdir, 'ED_violin_4groups.png');
    out_fig = fullfile(plot_outdir, 'ED_violin_4groups.fig');
    set(fig, 'PaperPositionMode','auto');
    print(fig, out_png, '-dpng', '-r300');
    saveas(fig, out_fig);
    fprintf('ED violin plot saved to:\n  %s\n  %s\n', out_png, out_fig);
end

%% 5d) ED statistics: Normality checks + 2x2 mixed ANOVA (Group x Condition)
% Assumes already in workspace:
%   ED_local   [numSubjects x 1]
%   ED_global  [numSubjects x 1]
%   young      [nY x 1 or 1 x nY] indices into ED vectors
%   older      [nO x 1 or 1 x nO] indices into ED vectors

% -----------------------------
% 1) Split ED into groups
% -----------------------------
numSubjects = numel(ED_local);

young_idx = young(:);
older_idx = older(:);

% Sanity check: indices valid
if any(young_idx < 1) || any(older_idx < 1) || ...
   any(young_idx > numSubjects) || any(older_idx > numSubjects)
    error('young/older contain indices outside 1..numSubjects. If these are IDs, map IDs -> indices first.');
end

Y_local_ED  = ED_local(young_idx);
O_local_ED  = ED_local(older_idx);
Y_global_ED = ED_global(young_idx);
O_global_ED = ED_global(older_idx);

fprintf('\nED means:\n');
fprintf(' Local  - Young mean = %.3f, Older mean = %.3f\n', ...
        mean(Y_local_ED,'omitnan'), mean(O_local_ED,'omitnan'));
fprintf(' Global - Young mean = %.3f, Older mean = %.3f\n', ...
        mean(Y_global_ED,'omitnan'), mean(O_global_ED,'omitnan'));

% Optional: check amount of missing values (important for rmANOVA)
fprintf('Missing ED values (NaN): Local=%d, Global=%d\n', ...
        sum(isnan(ED_local)), sum(isnan(ED_global)));

% -----------------------------
% 2) Normality checks (Anderson-Darling)
% -----------------------------
fprintf('\n--- Normality Testing for ED (Anderson-Darling) ---\n');

[h1,p1] = adtest(Y_local_ED(~isnan(Y_local_ED)));
fprintf('Local Young:  h=%d (1=non-normal), p=%.4f\n', h1, p1);

[h2,p2] = adtest(O_local_ED(~isnan(O_local_ED)));
fprintf('Local Older:  h=%d (1=non-normal), p=%.4f\n', h2, p2);

[h3,p3] = adtest(Y_global_ED(~isnan(Y_global_ED)));
fprintf('Global Young: h=%d (1=non-normal), p=%.4f\n', h3, p3);

[h4,p4] = adtest(O_global_ED(~isnan(O_global_ED)));
fprintf('Global Older: h=%d (1=non-normal), p=%.4f\n', h4, p4);

% -----------------------------
% 3) 2x2 mixed ANOVA: Group (Young/Older) x Condition (Local/Global)
% -----------------------------
% IMPORTANT limitation: fitrm/ranova uses listwise deletion if any repeated
% measure is missing. So we explicitly keep only subjects with BOTH Local
% and Global ED present, otherwise effective N can drop silently.

valid = isfinite(ED_local) & isfinite(ED_global);

if ~all(valid)
    fprintf('\nNote: %d/%d subjects have both Local and Global ED. Using only these for mixed ANOVA.\n', ...
            sum(valid), numSubjects);
end

ED_local_valid  = ED_local(valid);
ED_global_valid = ED_global(valid);

% Build group factor for valid subjects
group = strings(sum(valid),1);

% Map original indices -> valid-subset indices
valid_idx = find(valid);

% Mark group labels within the valid subset
group(ismember(valid_idx, young_idx)) = "Young";
group(ismember(valid_idx, older_idx)) = "Older";

% If someone is in neither group, that's a problem
if any(group == "")
    missingLabelSubs = valid_idx(group=="");
    error('Some valid subjects were not labeled Young/Older. Check young/older indices. Example subject indices: %s', mat2str(missingLabelSubs(1:min(10,end))'));
end

group = categorical(group);

% Table for repeated-measures model: two columns = within-subject factor
T = table(group, ED_local_valid, ED_global_valid, ...
          'VariableNames', {'Group','Local','Global'});

% Within-subject design table
Within = table({'Local'; 'Global'}, 'VariableNames', {'Condition'});

% Repeated-measures model: Group = between factor, Condition = within factor
rm = fitrm(T, 'Local,Global ~ Group', 'WithinDesign', Within);

% Mixed ANOVA table (includes Group, Condition, Group:Condition)
ranovatbl = ranova(rm, 'WithinModel', 'Condition');

% -------- Extract p-values from ranovatbl --------
% Row names usually include:
% '(Intercept)', 'Group', 'Error', '(Intercept):Condition', 'Group:Condition', 'Error(Condition)'
p_group       = ranovatbl{'Group',                 'pValue'};   % Group main effect
p_condition   = ranovatbl{'(Intercept):Condition', 'pValue'};   % Condition main effect
p_interaction = ranovatbl{'Group:Condition',       'pValue'};   % Interaction

% Bonferroni correction for 3 tests (Group, Condition, Interaction)
pvals = [p_group, p_condition, p_interaction];
pvals_corr = min(pvals * 3, 1);

fprintf('\n--- Mixed ANOVA (ED) ---\n');
fprintf('Group main effect:      p = %.4f (Bonferroni-corr = %.4f)\n', ...
        p_group, pvals_corr(1));
fprintf('Condition main effect:  p = %.4f (Bonferroni-corr = %.4f)\n', ...
        p_condition, pvals_corr(2));
fprintf('Group x Condition int.: p = %.4f (Bonferroni-corr = %.4f)\n', ...
        p_interaction, pvals_corr(3));

% Optional: display the full ANOVA table for transparency/debugging
disp(ranovatbl);

%% 5e) Cumulative varience explained per component for each group in each condition - calculation

%%% --------------------------------------------------------------- %%%
%%% Cumulative variance explained (PC1–PC100) for local & global
%%% --------------------------------------------------------------- %%%

% How many PCs to include on the x-axis
maxPC = 30;

% Sanity: make sure we don't exceed available PCs
nPC_local  = numel(BROADNESS_local{1}.Variance_BrainNetworks);
nPC_global = numel(BROADNESS_global{1}.Variance_BrainNetworks);
maxPC = min([maxPC, nPC_local, nPC_global]);  % in case vectors are shorter

pcRange = 1:maxPC;

% Preallocate: rows = PCs, columns = subjects
cumVar_local  = nan(maxPC, numSubjects);
cumVar_global = nan(maxPC, numSubjects);

for s = 1:numSubjects
    % Local variance explained (per PC) for subject s
    v_loc = BROADNESS_local{s}.Variance_BrainNetworks(:);
    v_glob = BROADNESS_global{s}.Variance_BrainNetworks(:);
    
    % Defensive: limit to maxPC in case length differs slightly
    nLoc  = min(maxPC, numel(v_loc));
    nGlob = min(maxPC, numel(v_glob));
    
    cumVar_local(1:nLoc,  s) = cumsum(v_loc(1:nLoc));
    cumVar_global(1:nGlob, s) = cumsum(v_glob(1:nGlob));
end

% Make sure group indices are column vectors
young_idx = young_subj(:);
older_idx = older_subj(:);

% (Optional) sanity checks 
% assert(all(young_idx >= 1 & young_idx <= numSubjects), 'young_subj has invalid indices');
% assert(all(older_idx >= 1 & older_idx <= numSubjects), 'older_subj has invalid indices');

% Split into groups: matrices [PC x subjects_in_group]
Y_local = cumVar_local(:, young_idx);
O_local = cumVar_local(:, older_idx);

Y_global = cumVar_global(:, young_idx);
O_global = cumVar_global(:, older_idx);

% Group sizes
nY = size(Y_local, 2);
nO = size(O_local, 2);

% Mean and SEM across subjects (per PC)
meanY_loc = mean(Y_local, 2, 'omitnan');
meanO_loc = mean(O_local, 2, 'omitnan');
semY_loc  = std(Y_local, 0, 2, 'omitnan') ./ sqrt(nY);
semO_loc  = std(O_local, 0, 2, 'omitnan') ./ sqrt(nO);

meanY_glob = mean(Y_global, 2, 'omitnan');
meanO_glob = mean(O_global, 2, 'omitnan');
semY_glob  = std(Y_global, 0, 2, 'omitnan') ./ sqrt(nY);
semO_glob  = std(O_global, 0, 2, 'omitnan') ./ sqrt(nO);

% Colors for plotting (MATLAB default-ish)
colYoung = [0 0.4470 0.7410];   % blue
colOlder = [0.8500 0.3250 0.0980]; % orange


%% 5f) Explained cumulative variance plot: Local vs Global & Young vs Older 

set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultTextFontName', 'Helvetica');

% ===================== USER PARAMETERS =====================
legendFontSize   = 12;     % change legend text size here
legendTokenSize  = [28 14];% change legend symbol size here: [lineLength height] (points)
legendLineWidth  = 2.5;    % change legend line thickness here (applies to plotted lines too)
% ===========================================================

% ED means:
% Local  - Young mean = 3.924, Older mean = 3.949
% Global - Young mean = 5.615, Older mean = 5.728

comb_colors = [ ...
    0.1 0.2 0.7;   % Young-Local   (blue)
    0.5 0.1 0.1;   % Old-Local     (red)
    0.4 0.6 0.85;  % Young-Global  (light blue)
    0.8 0.2 0.2];  % Old-Global    (light red)

colYL = comb_colors(1,:); % Young Local
colOL = comb_colors(2,:); % Older Local
colYG = comb_colors(3,:); % Young Global
colOG = comb_colors(4,:); % Older Global

figure; hold on;

% --- LOCAL: Young ---
fill([pcRange, fliplr(pcRange)], ...
     [(meanY_loc - semY_loc)', fliplr((meanY_loc + semY_loc)')], ...
     colYL, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
hYL = plot(pcRange, meanY_loc, 'Color', colYL, 'LineWidth', legendLineWidth, ...
           'LineStyle', '-');

% --- LOCAL: Older ---
fill([pcRange, fliplr(pcRange)], ...
     [(meanO_loc - semO_loc)', fliplr((meanO_loc + semO_loc)')], ...
     colOL, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
hOL = plot(pcRange, meanO_loc, 'Color', colOL, 'LineWidth', legendLineWidth, ...
           'LineStyle', '-');

% --- GLOBAL: Young ---
fill([pcRange, fliplr(pcRange)], ...
     [(meanY_glob - semY_glob)', fliplr((meanY_glob + semY_glob)')], ...
     colYG, 'FaceAlpha', 0.10, 'EdgeColor', 'none');
hYG = plot(pcRange, meanY_glob, 'Color', colYG, 'LineWidth', legendLineWidth, ...
           'LineStyle', '-');

% --- GLOBAL: Older ---
fill([pcRange, fliplr(pcRange)], ...
     [(meanO_glob - semO_glob)', fliplr((meanO_glob + semO_glob)')], ...
     colOG, 'FaceAlpha', 0.10, 'EdgeColor', 'none');
hOG = plot(pcRange, meanO_glob, 'Color', colOG, 'LineWidth', legendLineWidth, ...
           'LineStyle', '-');

xlabel('Principal component (PC)');
ylabel('Cumulative variance explained');
title('Cumulative variance explained: Local vs Global, Young vs Older');
grid on; box on;

% Axis limits
xlim([1 30]);
yMaxAll = max([meanY_loc + semY_loc; meanO_loc + semO_loc; ...
               meanY_glob + semY_glob; meanO_glob + semO_glob]);
ylim([30 yMaxAll * 1.05]);

% --- Add ED dotted lines (vertical) + legend entries ---

% Effective dimensionality (x-positions)
ED_YL = 3.924;  % Younger - Local
ED_OL = 3.949;  % Older   - Local
ED_YG = 5.615;  % Younger - Global
ED_OG = 5.728;  % Older   - Global

% Helper: y-value on curve at arbitrary x (linear interp)
y_at = @(x, yvec) interp1(pcRange, yvec(:), x, 'linear', 'extrap');

% Bottom of plot in data units (x-axis baseline of the plot)
yl = ylim;
y0 = yl(1);

% Vertical dotted lines (same color as corresponding curves)
hED_YL = plot([ED_YL ED_YL], [y0 y_at(ED_YL, meanY_loc)], ':', ...
    'Color', colYL, 'LineWidth', legendLineWidth);
hED_OL = plot([ED_OL ED_OL], [y0 y_at(ED_OL, meanO_loc)], ':', ...
    'Color', colOL, 'LineWidth', legendLineWidth);
hED_YG = plot([ED_YG ED_YG], [y0 y_at(ED_YG, meanY_glob)], ':', ...
    'Color', colYG, 'LineWidth', legendLineWidth);
hED_OG = plot([ED_OG ED_OG], [y0 y_at(ED_OG, meanO_glob)], ':', ...
    'Color', colOG, 'LineWidth', legendLineWidth);

% Update legend to include dotted lines
lgd = legend([hYL, hED_YL, hOL, hED_OL, hYG, hED_YG, hOG, hED_OG], ...
       {'Younger - Local ± SEM', ...
        'Effective Dimensionality', ...
        'Older - Local ± SEM', ...
        'Effective Dimensionality', ...
        'Younger - Global ± SEM', ...
        'Effective Dimensionality', ...
        'Older - Global ± SEM', ...
        'Effective Dimensionality'}, ...
        'Location', 'southeast');

% Apply legend sizing parameters
lgd.FontSize = legendFontSize;          % legend text size
lgd.ItemTokenSize = legendTokenSize;    % legend symbol size [lineLength height]


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
% Leonardo Bonetti:         leonardo.bonetti@clin.au.dk
%                           leonardo.bonetti@psych.ox.ac.uk
% Mattia Rosso:             mattia.rosso@clin.au.dk
% Chiara Malvaso:           chiara.malvaso@studio.unibo.it
% Mathias Houe Andersen     mathiasha@drcmr.dk
%
%  Please cite the first BROADNESS paper if using the BROADNESS toolbox:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science.
%  BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
%  https://doi.org/10.1002/advs.202507878
