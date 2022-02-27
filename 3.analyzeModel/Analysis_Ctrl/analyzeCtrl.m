%% Load Data
tic;
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to run 'compareCbModel'
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')

%% iAstro_iPS_TP_abs
%% iAstro_iPS_Ctrl_TP_abs vs iAstro_Primary_TP_abs
model_a = iAstro_iPS_Ctrl_TP_abs;
model_b = iAstro_Primary_TP_abs;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
Table_m = table(m);
writetable(Table_m,'FSR_iAstro_iPS_Ctrl_TP_abs_vs_iAstro_Primary_TP_abs.csv','WriteVariableNames', true, 'Delimiter','\t');
n = setdiff(dat.fluxSpanTable.Rxn, m);
Table_n = table(n);
writetable(Table_n,'UnChanged_iAstro_iPS_Ctrl_TP_abs_vs_iAstro_Primary_TP_abs.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_abs_vs_Primary_TP_abs = FluxDiff;
clear model_a model_b FluxDiff dat m n Table_m Table_n

%% iAstro_iPS_TP_norm
%% iAstro_iPS_Ctrl_TP_norm vs iAstro_Primary_TP_norm
model_a = iAstro_iPS_Ctrl_TP_norm;
model_b = iAstro_Primary_TP_norm;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
Table_m = table(m);
writetable(Table_m,'FSR_iAstro_iPS_Ctrl_TP_norm_vs_iAstro_Primary_TP_norm.csv','WriteVariableNames', true, 'Delimiter','\t');
n = setdiff(dat.fluxSpanTable.Rxn, m);
Table_n = table(n);
writetable(Table_n,'UnChanged_iAstro_iPS_Ctrl_TP_norm_vs_iAstro_Primary_TP_norm.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_norm_vs_Primary_TP_norm = FluxDiff;
clear model_a model_b FluxDiff dat m n Table_m Table_n

%%
clearvars -except FluxDiff_iPS_Ctrl_TP_abs_vs_Primary_TP_abs FluxDiff_iPS_Ctrl_TP_norm_vs_Primary_TP_norm

%% save data to mat file
save('iAstro_FluxDiff_iPSCtrl_vs_Primary.mat'); 

%%
toc;
