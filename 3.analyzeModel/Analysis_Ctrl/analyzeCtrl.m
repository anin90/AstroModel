%% Load Data
tic;
% Load iPS_Models
FileName   = 'iAstro_iPS_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
% Load Primary_Models
FileName   = 'iAstro_Primary_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to run 'compareCbModel'
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')

%% iAstro_iPS_TP
%% iAstro_Primary_TP vs iAstro_iPS_Ctrl_TP
model_a = iAstro_Primary_TP;
model_b = iAstro_iPS_Ctrl_TP;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
Table_m = table(m);
writetable(Table_m,'FSR_iAstro_Primary_TP_vs_iAstro_iPS_Ctrl_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
n = setdiff(dat.fluxSpanTable.Rxn, m);
Table_n = table(n);
writetable(Table_n,'UnChanged_iAstro_Primary_TP_vs_iAstro_iPS_Ctrl_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_Primary_TP_vs_iPS_Ctrl_TP = FluxDiff;
clear model_a model_b FluxDiff dat m n Table_m Table_n

%% iAstro_iPS_T
%% iAstro_Primary_T vs iAstro_iPS_Ctrl_T
model_a = iAstro_Primary_T;
model_b = iAstro_iPS_Ctrl_T;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
Table_m = table(m);
writetable(Table_m,'FSR_iAstro_Primary_T_vs_iAstro_iPS_Ctrl_T.csv','WriteVariableNames', true, 'Delimiter','\t');
n = setdiff(dat.fluxSpanTable.Rxn, m);
Table_n = table(n);
writetable(Table_n,'UnChanged_iAstro_Primary_T_vs_iAstro_iPS_Ctrl_T.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_Primary_T_vs_iPS_Ctrl_T = FluxDiff;
clear model_a model_b FluxDiff dat m n Table_m Table_n

%%
clearvars -except FluxDiff_Primary_TP_vs_iPS_Ctrl_TP FluxDiff_Primary_T_vs_iPS_Ctrl_T

%% save data to mat file
save('iAstro_FluxDiff_Primary_vs_iPSCtrl.mat'); 

%%
toc;
