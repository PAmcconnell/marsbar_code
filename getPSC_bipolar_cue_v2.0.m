
% NOTE: MARSBAR must be running
% Written by Matthew Brett (MARSBAR toolbox: https://marsbar-toolbox.github.io/faq.html)
% Adapted by Patrick McConnell and Taylor Salo
% TODO: Consider changing this to a function. The only inputs are spmFile, roiFiles, and firLength.

% Clear all variables from workspace
clear all

% Check that MARSBAR is running (i.e., maroi function is available) % REVIEW
if ~exist('maroi', 'file')
    error('MARSBAR must be running');
end

% Define output directory for percent signal change (PSC) results
outDir = 'E:\NIAAAr01\Results\cue_analysis\';

% Read subject list
subjectsList = readCsv('E:\NIAAAr01\fmri_data\subjList.csv');
subjects = removeEmptyCells(subjectsList{1,1}.col);
visits = {'V1', 'V2', 'V3'};
conditions = {'Alc', 'Bev'};

% Process each subject, visit, and ROI
for iSubj = 1:length(subjects)
    for iVisit = 1:length(visits)
        % Define the directory for SPM files
        spmDir = fullfile('E:\NIAAAr01\fMRI_data', subjects{iSubj}, visits{iVisit}, 'optcens', '01_censored');
        spmFile = fullfile(spmDir, 'SPM.mat');
        load(spmFile);

        % ROI path for the current visit
        roiDir = fullfile('E:\NIAAAr01\fMRI_data', subjects{iSubj}, visits{iVisit}, 'NSrois');
        roiFiles = dir(fullfile(roiDir, '*_roi.mat'));
        
        for iRoi = 1:length(roiFiles)
            roiPath = fullfile(roiFiles(iRoi).folder, roiFiles(iRoi).name);
            R = maroi(roiPath);
            
            % Initialize output structure
            outStruct = cell(1, 2 + length(conditions));
            outStruct{1} = {'Subject', subjects{iSubj}};
            outStruct{2} = {'Visit', visits{iVisit}};
            colCounter = 3;
            
            % Marsbar specific code from https://marsbar-toolbox.github.io/faq.html
            D = mardo(spmFile);
            Y = get_marsy(R, D, 'mean');
            E = estimate(D, Y);

            [e_specs, e_names] = event_specs(E);
            ets = event_types_named(E);

            % Extract data for relevant conditions
            for jCond = 1:length(conditions)
                idx = find(strcmp(e_names, conditions{jCond}));
                dur = 24 % Define condition duration in seconds
                if ~isempty(idx)
                    psc = event_signal(E, e_specs(:, idx), dur);
                    outStruct{colCounter} = {conditions{jCond}, psc};
                else
                    outStruct{colCounter} = {conditions{jCond}, NaN};
                end
                colCounter = colCounter + 1;
            end
            
            % Output filename
            [~, roiName, ~] = fileparts(roiPath);
            outputFile = fullfile(outDir, [subjects{iSubj} '_' visits{iVisit} '_' roiName '_cue.csv']);
            writeCsv(outStruct, outputFile);
        end
    end
end
