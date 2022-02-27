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
%% iAstro_iPS_Ctrl_TP vs iAstro_iPS_BD_TP
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_BD_TP;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP vs iAstro_iPS_BD_R_TP
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_BD_R_TP;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_R_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_R_TP = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP vs iAstro_iPS_BD_NR_TP
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_BD_NR_TP;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_NR_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_NR_TP = FluxDiff;
clear model_a model_b dat m T FluxDiff


%% iAstro_iPS_T
%% iAstro_iPS_Ctrl_T vs iAstro_iPS_BD_T
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_BD_T;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_T_vs_iAstro_iPS_BD_T.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_T_vs_iPS_BD_T = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_T vs iAstro_iPS_BD_R_T
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_BD_R_T;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_T_vs_iAstro_iPS_BD_R_T.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_T_vs_iPS_BD_R_T = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_T vs iAstro_iPS_BD_NR_T
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_BD_NR_T;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_T_vs_iAstro_iPS_BD_NR_T.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_T_vs_iPS_BD_NR_T = FluxDiff;
clear model_a model_b dat m T FluxDiff

%%
clearvars -except FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_R_TP FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_NR_TP FluxDiff_iPS_Ctrl_T_vs_iPS_BD_T FluxDiff_iPS_Ctrl_T_vs_iPS_BD_R_T FluxDiff_iPS_Ctrl_T_vs_iPS_BD_NR_T

%% save data to mat file
save('iAstro_FluxDiff_BD.mat'); 

%%
toc;
