% getPSC_bipolar_cue.m
%
% This script calculates the percent signal change (PSC) for regions of interest (ROIs) 
% from fMRI data for subjects across multiple visits and conditions using the MarsBaR 
% toolbox for ROI analysis of SPM data.
%
% Requirements:
% - MarsBaR must be running (https://marsbar-toolbox.github.io/faq.html)
% - ROIs must be saved in .mat format (not .nii format).
%
% Inputs:
% - No direct function inputs, but the script reads subject information and file paths
%   from 'E:\NIAAAr01\fmri_data\subjList.csv' and processes corresponding SPM.mat files 
%   and ROI.mat files located in structured directories.
%
% Outputs:
% - CSV files containing percent signal change data for each subject, visit, and ROI.
%   These files are saved in the directory 'E:\NIAAAr01\Results\cue_analysis\'.
%   Each file is named according to the subject ID, visit, and ROI name.
%
% Process Overview:
% 1. Loads subject list and initializes variables for visits and conditions.
% 2. Loops through each subject, visit, and ROI file:
%    - Loads the SPM.mat file for the specific visit of the subject.
%    - Loads each ROI file from the visit-specific directory.
%    - Calculates the mean PSC for specified conditions ('Alc' and 'Bev').
% 3. Outputs results to a CSV file for each ROI, including the subject ID and visit.
%
% Notes:
% - The script assumes a specific directory structure and file naming convention.
% - It's recommended to check paths and file availability before execution.
% - Consider adapting this script into a function for more versatile use cases.
%
% Written by Matthew Brett and adapted by Patrick A. McConnell and Taylor Salo.
% Last updated: [20240419]
%
% See also maroi, mardo, get_marsy, estimate, event_specs, event_types_named, writeCsv
% https://marsbar-toolbox.github.io/faq.html
%
% Example:
% Assuming the necessary files are in place and MarsBaR is active, running this script
% in the MATLAB command window will process all data as configured and output CSV files
% to the specified directory.

% Clear all variables from workspace
clear all

% Check that MARSBAR is running (i.e., maroi function is available) % REVIEW
if ~exist('maroi', 'file')
    error('MARSBAR must be running');
end

% Define output directory for percent signal change (PSC) results
outDir = 'E:\NIAAAr01\Results\cue_analysis\';
if ~exist(outDir, 'dir')
    mkdir(outDir); % Create directory if it does not exist
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

% Process each subject, visit, and ROI
for iSubj = 1:length(subjects)
    for iVisit = 1:length(visits)
        % Define the directory for SPM files
        spmDir = fullfile('E:\NIAAAr01\fMRI_data', subjects{iSubj}, visits{iVisit}, 'optcens', '01_censored');
        if ~exist(spmDir, 'dir')
            warning(['SPM directory does not exist: ', spmDir]);
            continue; % Skip to next iteration
        end
        spmFile = fullfile(spmDir, 'SPM.mat');
        if ~exist(spmFile, 'file')
            warning(['SPM file does not exist: ', spmFile]);
            continue; % Skip to next iteration
        end
        load(spmFile);

        % ROI path for the current visit
        roiDir = fullfile('E:\NIAAAr01\fMRI_data', subjects{iSubj}, visits{iVisit}, 'NSrois');
        if ~exist(roiDir, 'dir')
            warning(['ROI directory does not exist: ', roiDir]);
            continue; % Skip to next iteration
        end
        roiFiles = dir(fullfile(roiDir, '*_roi.mat'));
        
        for iRoi = 1:length(roiFiles)
            roiPath = fullfile(roiFiles(iRoi).folder, roiFiles(iRoi).name);
            if ~exist(roiPath, 'file')
                warning(['ROI file does not exist: ', roiPath]);
                continue; % Skip to next iteration
            end
            R = maroi(roiPath);
            
            % Initialize output structure
            outStruct = cell(1, 2 + length(conditions));
            outStruct{1} = {'Subject', subjects{iSubj}};
            outStruct{2} = {'Visit', visits{iVisit}};
            colCounter = 3;
            
            % Marsbar specific code
            D = mardo(spmFile);
            Y = get_marsy(R, D, 'mean');
            E = estimate(D, Y);

            [e_specs, e_names] = event_specs(E);

            % Extract data for relevant conditions
            for jCond = 1:length(conditions)
                idx = find(strcmp(e_names, conditions{jCond}));
                dur = 24; % Define condition duration in seconds
                if isempty(idx)
                    outStruct{colCounter} = {conditions{jCond}, NaN};
                else
                    psc = event_signal(E, e_specs(:, idx), dur);
                    outStruct{colCounter} = {conditions{jCond}, psc};
                end
                colCounter = colCounter + 1;
            end
            
            % Output filename
            [~, roiName, ~] = fileparts(roiPath);
            outputFile = fullfile(outDir, [subjects{iSubj} '_' visits{iVisit} '_' roiName '_cue.csv']);
            if ~writeCsv(outStruct, outputFile)
                warning(['Failed to write output file: ', outputFile]);
            end
        end
    end
end
