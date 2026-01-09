%% Create .func files

% --- Config ---
paths = {'/yourpath/folder_with_nifti_files'};

% Mid thickness
% Download these files also from the GitHub repository
Lsurf = '/yourpath/folder_with_nifti_files/ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii';
Rsurf = '/yourpath/folder_with_nifti_files/ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii';

% Absolute, verified path to wb_command (change if needed)
wb = '/yourpath//workbench-rh_linux64-v2.1.0/workbench/exe_rh_linux64/wb_command';

assert(exist(wb,'file')==2, 'wb_command not found at: %s', wb);

% Process ALL defined paths instead of only 20
idx = 1:numel(paths);   % was: idx = 20;

madeL = 0; madeR = 0;

for ii = idx
    p = paths{ii};
    if ~isfolder(p); warning('Not a folder: %s', p); continue; end

    list = [dir(fullfile(p,'*.nii')); dir(fullfile(p,'*.nii.gz'))];
    fprintf('Path %d: %s — %d files found\n', ii, p, numel(list));
    if isempty(list); continue; end

    for jj = 1:numel(list)
        inima = fullfile(list(jj).folder, list(jj).name);

        % Robust stem
        nm = list(jj).name;
        if endsWith(nm, '.nii.gz'); stem = erase(nm,'.nii.gz');
        elseif endsWith(nm, '.nii'); stem = erase(nm,'.nii');
        else; warning('Skipping non-NIfTI: %s', inima); continue; end

        outL = fullfile(list(jj).folder, [stem '_L.func.gii']);
        outR = fullfile(list(jj).folder, [stem '_R.func.gii']);

        cmdL = sprintf('"%s" -volume-to-surface-mapping "%s" "%s" "%s" -trilinear', wb, inima, Lsurf, outL);
        cmdR = sprintf('"%s" -volume-to-surface-mapping "%s" "%s" "%s" -trilinear', wb, inima, Rsurf, outR);

        [stL, msgL] = system(cmdL); if stL~=0, error('wb_command Left failed:\n%s', msgL); end
        [stR, msgR] = system(cmdR); if stR~=0, error('wb_command Right failed:\n%s', msgR); end

        madeL = madeL + exist(outL,'file')==2;
        madeR = madeR + exist(outR,'file')==2;
        fprintf('OK: %s\n', stem);
    end

    % sanity check per path
    nL = numel(dir(fullfile(p,'*_L.func.gii')));
    nR = numel(dir(fullfile(p,'*_R.func.gii')));
    fprintf('Outputs in %s -> L:%d  R:%d\n', p, nL, nR);
end

fprintf('TOTAL created -> L:%d  R:%d\n', madeL, madeR);


%% Viewing in workbench

% wb_view
% load surface: 100307.L.midthickness.32k_fs_LR.surf.gii
% load functonal: PCA_ActivationPattern_BrainNetwork_01_L.func.gii