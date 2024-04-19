% To extract mean PSC from 1st level fMRI models using SPM/Marsbar toolbox.

% The output of this script (i.e., PSC for each design matrix vector) is NOT listed as in the design matrix,
% but sorted ALPHABETICALLY


% MARSBAR must be running for script to execute

% Written by Matthew Brett (MARSBAR toolbox)
% Adapted by Patrick McConnell and Taylor Salo

clear all
outDir = '/TRAIN_Data/fmri/smoke_emote/tasks/er_task/z_project/';

subjectsList = readCsv('/TRAIN_Data/fmri/smoke_emote/tasks/er_task/z_project/smoke_emote.csv');
allSubjects = {};

subjects = removeEmptyCells(subjectsList{1,1}.col);
version = removeEmptyCells(subjectsList{1,2}.col);
subjectsDir = removeEmptyCells(subjectsList{1,3}.col);
status = removeEmptyCells(subjectsList{1,4}.col);
rootDir = '/TRAIN_Data/fmri/';
task = '/tasks/er_task/subjects/';

% These conditions must be listed below according to the alphabetical order defined in the 1LM to correspond
% with the output (not in the order listed in the design matrix)

conditions= {'NR' 'NV' 'neutral' 'PR' 'PV'};

% ROI must be in .mat format, not .nii

roiFiles = {
'/TRAIN_Data/fmri/smoke_emote/tasks/er_task/z_models/081518/mcconnell/negativeView_OST_ns/R_AMY_nm_roi.mat'};

% roiFiles = {'/TRAIN_Data/rois/NAcc/RVS_sphere_3-12_8_-8_roi.mat'
% '/TRAIN_Data/rois/NAcc/LVS_sphere_3--12_8_-8_roi.mat'
% '/TRAIN_Data/rois/NAcc/Left_pCore_roi.mat'            
% '/TRAIN_Data/rois/NAcc/Left_pShell_roi.mat'
% '/TRAIN_Data/rois/NAcc/Right_pCore_roi.mat'          
% '/TRAIN_Data/rois/NAcc/Right_pShell_roi.mat'
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_AStr_L_MNI_roi.mat'
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_AStr_R_MNI_roi.mat'
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_CM_L_MNI_roi.mat'  
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_CM_R_MNI_roi.mat'  
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_LB_L_MNI_roi.mat'  
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_LB_R_MNI_roi.mat'  
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_SF_L_MNI_roi.mat'  
% '/TRAIN_Data/rois/amygdala/ROI_Amygdala_SF_R_MNI_roi.mat'};


%%

% Begin ROI loop
for iRoi = 1:length(roiFiles)
    % set up outStruct
    outStruct{1}.header{1} = 'Subject';
    outStruct{2}.header{1} = 'Status';
    outStruct{3}.header{1} = 'Study';
    colCounter = 4;
    for jCond = 1:length(conditions)
        outStruct{colCounter}.header{1} = [conditions{jCond}];
        colCounter = colCounter + 1;
    end
    
    % Make marsbar ROI object
    R = maroi(roiFiles{iRoi});
    
    % Begin Subject Loop
    for jSubj = 1:length(subjects)
        outStruct{1}.col{jSubj} = subjects{jSubj};
        spmFile = fullfile(rootDir, subjectsDir{jSubj}, task, subjects{jSubj}, version{jSubj}, '/spm_er_neuromorph_noSmooth/SPM.mat');
        load(spmFile);
        
        % Make marsbar design object
        D = mardo(SPM);
        % Fetch Data into marsbar data object
        Y = get_marsy(R, D, 'mean');
        % Estimate design on ROI data
        E = estimate(D, Y);
        % Get definitions of all events in model
        [e_specs, e_names] = event_specs(E);
        n_events = size(e_specs, 2);
        % Get contrasts from original design
        xCon = get_contrasts(D);
        % Put contrasts from original design back into design object
        E = set_contrasts(E, xCon);
        % get design betas
        b = betas(E);
        % Get compound event types structure
        ets = event_types_named(E);
        n_event_types = length(ets);
        % get stats and stuff for all contrasts into statistics structure
        marsS = compute_contrasts(E, 1:length(xCon));
        colCounter = 4;
        
        % Begin condition loop
        
        for e_s = 1:n_event_types
            %  define block duration or zero for event-related
            dur = 10;
            
            % calculate PSC for trial across conditions
            
            pct_ev(:,e_s) = event_signal(E, ets(e_s).e_spec, dur);
         
            % Print results to outStruct
         
            outStruct{colCounter}.col{jSubj} = pct_ev(e_s);
            outStruct{2}.col{jSubj} = status{jSubj};
            outStruct{3}.col{jSubj} = subjectsDir{jSubj};
            colCounter = colCounter + 1;
        end
    end

    %%
    
    % Print outStruct to csv file
    [~, roiName, ~] = fileparts(roiFiles{iRoi});
    
    % IMPORTANT: directory needs to have been created in advance or there will
    % be an 'error using fprintf'
    writeCsv(outStruct, [outDir '/' roiName '_meanPSC_noSmooth_10sec.csv']);
end

