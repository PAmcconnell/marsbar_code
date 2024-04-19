% The output of this script (i.e., PSC for each design matrix vector) is NOT listed as in the design matrix, 
% but sorted ALPHABETICALLY
% MARSBAR must be running
% Written by Matthew Brett (MARSBAR toolbox)
% Adapted by Patrick McConnell and Taylor Salo
% TODO
% Consider changing this to a function. The only inputs are spmFile, roiFiles, and firLength.

clear all
outDir = '/TRAIN_Data/fmri/smoke_emote/tasks/er_task/z_project/';

subjectsList = readCsv('/TRAIN_Data/fmri/smoke_emote/tasks/er_task/z_project/smoke_emote.csv');
allSubjects = {}; 

subjects = removeEmptyCells(subjectsList{1,1}.col);
version = removeEmptyCells(subjectsList{1,2}.col);
subjectsDir = removeEmptyCells(subjectsList{1,3}.col);
status = removeEmptyCells(subjectsList{1,4}.col);
rootDir = '/TRAIN_Data/fmri/';
task = '/er_task/spm_er_neuromorph_noSmooth/';  
conditions = {'FreqGo_Inc' 'InfGo_Correct' 'InfGo_Inc' 'NG_Correct' 'NG_Inc'};

% ROI must be in .mat format, not .nii
roiFiles = {};

nBins = 12;

for iRoi = 1:length(roiFiles)
   % set up outStruct
   outStruct{1}.header{1} = 'Subject';
   outStruct{2}.header{1} = 'Status';
   outStruct{3}.header{1} = 'Study';
   colCounter = 4;
   for jCond = 1:length(conditions)
      for kBin = 1:nBins
         outStruct{colCounter}.header{1} = [conditions{jCond} '_' num2str(kBin)];
         colCounter = colCounter + 1;
      end
      outStruct{colCounter}.header{1} = [conditions{jCond} '_mean'];
      colCounter = colCounter + 1;
   end

   R = maroi(roiFiles{iRoi});

   for jSubj = 1:length(subjects)
      outStruct{1}.col{jSubj} = subjects{jSubj};
%       outStruct{2}.col{jSubj} = status{jSubj};
%       outStruct{3}.col{jSubj} = subjectsDir{jSubj};
      spmFile = fullfile(rootDir, subjectsDir{jSubj}, version{jSubj}, subjects{jSubj}, task, 'SPM.mat');
      load(spmFile);
      
      D = mardo(SPM); 
      Y = get_marsy(R, D, 'mean');
      E = estimate(D, Y);

      ets = event_types_named(E);
      nConditionsFound = length(ets);

      % binSize = TR
      binSize = SPM.xY.RT;

      options = struct('single', 1, 'percent', 1);
   
      colCounter = 4;
      for kCond = 1:length(conditions)
         conditionIndex = [];
         for mCond = 1:length(ets)
            conditionName = ets(mCond).name;
            if strcmp(conditionName, conditions{kCond})
               conditionIndex = mCond;
            end
         end
         if ~isempty(conditionIndex)
            fir_tc = event_fitted_fir(E, ets(conditionIndex).e_spec, binSize, nBins, options);
            for mBin = 1:nBins
               outStruct{colCounter}.col{jSubj} = fir_tc(mBin);
               outStruct{2}.col{jSubj} = status{jSubj};
               outStruct{3}.col{jSubj} = subjectsDir{jSubj};
               colCounter = colCounter + 1;
            end
            outStruct{colCounter}.col{jSubj} = mean(fir_tc);
            colCounter = colCounter + 1;
         else
            for mBin = 1:nBins
               outStruct{colCounter}.col{jSubj} = '';
               colCounter = colCounter + 1;
            end
            outStruct{colCounter}.col{jSubj} = '';
            colCounter = colCounter + 1;
         end
      end
      clear SPM
   end
   [~, roiName, ~] = fileparts(roiFiles{iRoi});
   writeCsv(outStruct, [outDir '/' roiName '_FIR.csv']);
   clear outStruct
   end