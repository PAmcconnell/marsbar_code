
% The output of this script (i.e., PSC for each design matrix vector) is NOT listed as in the design matrix, 
% but sorted ALPHABETICALLY
% MARSBAR must be running
% Written by Matthew Brett (MARSBAR toolbox)
% Adapted by Patrick McConnell and Taylor Salo
% TODO
% Consider changing this to a function. The only inputs are spmFile, roiFiles, and firLength.

clear all
outDir = '/TRAIN_Data/fmri/bipolar_alcohol/z_project/cue/';

subjectsList = readCsv('/TRAIN_Data/fmri/bipolar_alcohol/z_project/cue/bipolar_alcohol_subjList.csv');
allSubjects = {}; 

subjects = removeEmptyCells(subjectsList{1,1}.col);
version = removeEmptyCells(subjectsList{1,2}.col);
subjectsDir = removeEmptyCells(subjectsList{1,3}.col);
status = removeEmptyCells(subjectsList{1,4}.col);
rootDir = '/TRAIN_Data/fmri/';
task = 'functional/cue/spm_cue_noST/';  
conditions= {'Blur' 'Alcohol' 'Beverage' 'Rest'};

% ROI must be in .mat format, not .nii
roiFiles = {'/TRAIN_Data/rois/studies/bipolar_alcohol/gross_BA_PF_mask_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/IBASPM71_L_Nacc_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/IBASPM71_R_Nacc_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/LDS_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/L_IFG_op_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/L_INS_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/L_MFG_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/LVS_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/mPFC_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/RDS_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/R_IFG_op_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/R_INS_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/R_MFG_5mm_roi.mat'
'/TRAIN_Data/rois/studies/bipolar_alcohol/RVS_5mm_roi.mat'};

for iRoi = 1:length(roiFiles)
% set up outStruct
   outStruct{1}.header{1} = 'Subject';
   outStruct{2}.header{1} = 'Status';
   outStruct{3}.header{1} = 'Study';
   colCounter = 4;
   for jCond = 1:length(conditions)
      outStruct{colCounter}.header{1} = [conditions{jCond} '_mean'];
      colCounter = colCounter + 1;
   end

   R = maroi(roiFiles{iRoi});

   for jSubj = 1:length(subjects)
      outStruct{1}.col{jSubj} = subjects{jSubj};
      spmFile = fullfile(rootDir, subjectsDir{jSubj}, version{jSubj}, subjects{jSubj}, task, 'SPM.mat');
      load(spmFile);

      D = mardo(SPM); 
      Y = get_marsy(R, D, 'mean');
      E = estimate(D, Y);

      [e_specs, e_names] = event_specs(E);
      n_events = size(e_specs, 2);
      ets = event_types_named(E);
      colCounter = 4;
      
      for kCond = 1:length(conditions)
            dur = 24;
            e_s = ets(1,kCond).e_spec(2);
            psc = event_signal(E,  e_specs(:,e_s), dur);
            outStruct{colCounter}.col{jSubj} = psc;
            outStruct{2}.col{jSubj} = status{jSubj};
            outStruct{3}.col{jSubj} = subjectsDir{jSubj};
            colCounter = colCounter + 1;   
      end
   end
%        clear SPM
% end  
   [~, roiName, ~] = fileparts(roiFiles{iRoi});
   writeCsv(outStruct, [outDir '/' roiName '_cue.csv']);
%    clear outStruct
end
