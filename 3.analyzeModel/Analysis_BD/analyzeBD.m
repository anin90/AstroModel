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
%% iAstro_iPS_Ctrl_TP_abs vs iAstro_iPS_BD_TP_abs
model_a = iAstro_iPS_Ctrl_TP_abs;
model_b = iAstro_iPS_BD_TP_abs;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_abs_vs_iAstro_iPS_BD_TP_abs.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_TP_abs = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP_abs vs iAstro_iPS_BD_R_TP_abs
model_a = iAstro_iPS_Ctrl_TP_abs;
model_b = iAstro_iPS_BD_R_TP_abs;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_abs_vs_iAstro_iPS_BD_R_TP_abs.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_R_TP_abs = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP_abs vs iAstro_iPS_BD_NR_TP_abs
model_a = iAstro_iPS_Ctrl_TP_abs;
model_b = iAstro_iPS_BD_NR_TP_abs;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_abs_vs_iAstro_iPS_BD_NR_TP_abs.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_NR_TP_abs = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_TP_norm_t1
%% iAstro_iPS_Ctrl_TP_norm_t1 vs iAstro_iPS_BD_TP_norm_t1
model_a = iAstro_iPS_Ctrl_TP_norm_t1;
model_b = iAstro_iPS_BD_TP_norm_t1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_norm_t1_vs_iAstro_iPS_BD_TP_norm_t1.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_TP_norm_t1 = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP_norm_t1 vs iAstro_iPS_BD_R_TP_norm_t1
model_a = iAstro_iPS_Ctrl_TP_norm_t1;
model_b = iAstro_iPS_BD_R_TP_norm_t1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_norm_t1_vs_iAstro_iPS_BD_R_TP_norm_t1.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_R_TP_norm_t1 = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP_norm_t1 vs iAstro_iPS_BD_NR_TP_norm_t1
model_a = iAstro_iPS_Ctrl_TP_norm_t1;
model_b = iAstro_iPS_BD_NR_TP_norm_t1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_norm_t1_vs_iAstro_iPS_BD_NR_TP_norm_t1.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_NR_TP_norm_t1 = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_TP_norm_t2
%% iAstro_iPS_Ctrl_TP_norm_t2 vs iAstro_iPS_BD_TP_norm_t2
model_a = iAstro_iPS_Ctrl_TP_norm_t2;
model_b = iAstro_iPS_BD_TP_norm_t2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_TP_norm_t2.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_norm_t2_vs_iPS_BD_TP_norm_t2 = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP_norm_t2 vs iAstro_iPS_BD_R_TP_norm_t2
model_a = iAstro_iPS_Ctrl_TP_norm_t2;
model_b = iAstro_iPS_BD_R_TP_norm_t2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_R_TP_norm_t2.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_norm_t2_vs_iPS_BD_R_TP_norm_t2 = FluxDiff;
clear model_a model_b dat m T FluxDiff

%% iAstro_iPS_Ctrl_TP_norm_t2 vs iAstro_iPS_BD_NR_TP_norm_t2
model_a = iAstro_iPS_Ctrl_TP_norm_t2;
model_b = iAstro_iPS_BD_NR_TP_norm_t2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_NR_TP_norm_t2.csv','WriteVariableNames', true, 'Delimiter','\t');
FluxDiff_iPS_Ctrl_TP_norm_t2_vs_iPS_BD_NR_TP_norm_t2 = FluxDiff;
clear model_a model_b dat m T FluxDiff

%%
clearvars -except FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_TP_abs FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_R_TP_abs FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_NR_TP_abs FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_TP_norm_t1 FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_R_TP_norm_t1 FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_NR_TP_norm_t1 FluxDiff_iPS_Ctrl_TP_norm_t2_vs_iPS_BD_TP_norm_t2 FluxDiff_iPS_Ctrl_TP_norm_t2_vs_iPS_BD_R_TP_norm_t2 FluxDiff_iPS_Ctrl_TP_norm_t2_vs_iPS_BD_NR_TP_norm_t2

%% save data to mat file
save('iAstro_FluxDiff_BD.mat'); 

%%
toc;
