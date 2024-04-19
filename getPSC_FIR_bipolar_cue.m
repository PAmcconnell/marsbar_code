% getPSC_FIR_bipolar_cue.m
%
% This MATLAB script calculates the percent signal change (PSC) and extracts Finite Impulse Response (FIR)
% time courses for regions of interest (ROIs) across multiple visits and conditions using the MarsBaR toolbox.
% It processes fMRI data for a set list of subjects, reading their corresponding SPM.mat files and converts
% .nii ROI files to .mat format if necessary before computing PSC and FIR for predefined conditions.
%
% Requirements:
%   - MarsBaR must be running (https://marsbar-toolbox.github.io/faq.html)
%   - MATLAB with SPM installed for file operations
%   - readCsv.m & writeCsv.m functions for reading and writing CSV files
%
% Inputs:
%   - Subject information from 'E:\NIAAAr01\fmri_data\subjList.csv'
%   - Corresponding SPM.mat and ROI.mat files from a structured directory. Automatically converts .nii ROI files
%     to .mat if necessary using MarsBaR conversion functions.
%
% Outputs:
%   - CSV files containing percent signal change data and FIR time courses for each subject, visit, and ROI.
%     Files are saved in 'E:\NIAAAr01\Results\cue_analysis\' and named according to the subject ID, visit, and ROI name.
%
% Process:
%   1. Initialize FIR parameters, define output directory, and read the subject list.
%   2. Check and convert .nii ROI files to .mat format if necessary.
%   3. Process each subject by reading their SPM.mat file for each visit.
%   4. Load and process each ROI.mat file to extract both PSC and FIR time courses for predefined conditions.
%   5. Output results to CSV files.
%
% Usage Notes:
%   - Ensure MarsBaR and SPM are active and properly configured before running this script.
%   - Verify the directory structure and file existence to prevent runtime errors.
%
% Example:
%   Execute this script in MATLAB after configuring MarsBaR and organizing the data according to the described structure.
%
% Author:
%   Adapted by Patrick A. McConnell and Taylor Salo
%   Original MarsBaR code by Matthew Brett
%
% Last updated:
%   [20240419]
%
% See also:
%   maroi, mardo, get_marsy, estimate, event_specs, event_types_named, writeCsv, event_fitted_fir, spm_vol

% Clear all variables from workspace
clear all;

% Check that MARSBAR is running (i.e., maroi function is available)
if ~exist('maroi', 'file')
    error('MARSBAR must be running');
end

% Define output directory for results
outDir = 'E:\NIAAAr01\Results\cue_analysis\';
if ~exist(outDir, 'dir')
    mkdir(outDir);  % Create directory if it does not exist
end

% Attempt to read subject list
subjectListPath = 'E:\NIAAAr01\fmri_data\subjList.csv';
if ~exist(subjectListPath, 'file')
    error(['Subject list file does not exist: ', subjectListPath]);
end
subjectsList = readCsv(subjectListPath);
subjects = removeEmptyCells(subjectsList{1,1}.col);

% Define visits
visits = {'V1', 'V2', 'V3'};

% Define conditions - NOTE output (i.e., PSC for each design matrix vector) is NOT listed as in the design matrix, 
% but sorted ALPHABETICALLY
conditions = {'Alc', 'Bev'};

% Initialize PSC and FIR parameters
tr = 2;  % Repetition time (TR) in seconds
bin_size = tr;  % Duration of each FIR bin in seconds (i.e., TR)
dur = 24;  % Duration of the event in seconds
fir_length = dur;  % Total duration of the FIR analysis in seconds
n_bins = fir_length / bin_size;  % Number of FIR bins

% Process each subject, visit, and ROI
for iSubj = 1:length(subjects)
    for iVisit = 1:length(visits)
        spmDir = fullfile('E:\NIAAAr01\fMRI_data', subjects{iSubj}, visits{iVisit}, 'optcens', '01_censored');
        if ~exist(spmDir, 'dir')
            warning(['SPM directory does not exist: ', spmDir]);
            continue;
        end
        spmFile = fullfile(spmDir, 'SPM.mat');
        if ~exist(spmFile, 'file')
            warning(['SPM file does not exist: ', spmFile]);
            continue;
        end
        load(spmFile);

        % Here roiDir is labeled 'NSrois', the script is agnostic to roi filenames
        roiDir = fullfile('E:\NIAAAr01\fMRI_data', subjects{iSubj}, visits{iVisit}, 'NSrois');
        if ~exist(roiDir, 'dir')
            warning(['ROI directory does not exist: ', roiDir]);
            continue;
        end

        % Check and convert .nii ROI files to _roi.mat format if necessary
        niiFiles = dir(fullfile(roiDir, '*.nii'));
        for niiFile = niiFiles'
            matFilename = [strtok(niiFile.name, '.'), '_roi.mat'];
            matPath = fullfile(roiDir, matFilename);
            if ~exist(matPath, 'file')  % If the .mat file does not exist, convert the .nii file
                img = fullfile(roiDir, niiFile.name);
                o = maroi_image(struct('vol', spm_vol(img), 'binarize', 0, 'func', 'img'));
                o = maroi_matrix(o);
                saveroi(o, matPath);
                fprintf('Converted %s to %s\n', niiFile.name, matFilename);
            end
        end

        roiFiles = dir(fullfile(roiDir, '*_roi.mat'));
        
        for iRoi = 1:length(roiFiles)
            roiPath = fullfile(roiFiles(iRoi).folder, roiFiles(iRoi).name);
            if ~exist(roiPath, 'file')
                warning(['ROI file does not exist: ', roiPath]);
                continue;
            end
            R = maroi(roiPath);
            D = mardo(spmFile);
            Y = get_marsy(R, D, 'mean');
            E = estimate(D, Y);

            [e_specs, e_names] = event_specs(E);

            % Initialize output structure for both PSC and FIR
            outStruct = cell(1, 2 + length(conditions) * (n_bins + 1) + length(conditions));
            outStruct{1} = {'Subject', subjects{iSubj}};
            outStruct{2} = {'Visit', visits{iVisit}};
            colCounter = 3;

            % Extract data for PSC and FIR
            for jCond = 1:length(conditions)
                idx = find(strcmp(e_names, conditions{jCond}));
                if isempty(idx)
                    for k = 1:n_bins
                        outStruct{colCounter} = {sprintf('%s_%d', conditions{jCond}, k), NaN};
                        colCounter = colCounter + 1;
                    end
                    outStruct{colCounter} = {sprintf('%s_mean', conditions{jCond}), NaN};
                    colCounter = colCounter + 1;
                    outStruct{colCounter} = {conditions{jCond}, NaN};
                    colCounter = colCounter + 1;
                else
                    opts = struct('single', 1, 'percent', 1);
                    fir_tc = event_fitted_fir(E, e_specs(:, idx), bin_size, n_bins, opts);
                    for k = 1:n_bins
                        outStruct{colCounter} = {sprintf('%s_%d', conditions{jCond}, k), fir_tc(k)};
                        colCounter = colCounter + 1;
                    end
                    mean_fir = mean(fir_tc);
                    outStruct{colCounter} = {sprintf('%s_mean', conditions{jCond}), mean_fir};
                    colCounter = colCounter + 1;
                    psc = event_signal(E, e_specs(:, idx), dur);
                    outStruct{colCounter} = {conditions{jCond}, psc};
                    colCounter = colCounter + 1;
                end
            end

            % Output filename for FIR and PSC
            [~, roiName, ~] = fileparts(roiPath);
            outputFile = fullfile(outDir, sprintf('%s_%s_%s_FIR_PSC.csv', subjects{iSubj}, visits{iVisit}, roiName));
            writeCsv(outStruct, outputFile);
        end
    end
end
