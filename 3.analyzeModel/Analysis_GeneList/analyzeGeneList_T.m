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

%% Retain Primary and iPS_Ctrl models
clearvars -except iAstro_Primary_T iAstro_iPS_Ctrl_T

%% add path to run 'compareCbModel', 'downregulateGene'
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')

%% iAstro_Primary_T
%% CKMT2 (1160)
[~, GeneToRxnList, iAstro_Primary_T_CKMT2] = downregulateGene(iAstro_Primary_T, '1160');
model_a = iAstro_Primary_T;
model_b = iAstro_Primary_T_CKMT2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_T_CKMT2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_T_CKMT2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_T_CKMT2 

%% INPP5A (3632)
[~, GeneToRxnList, iAstro_Primary_T_INPP5A] = downregulateGene(iAstro_Primary_T, '3632');
model_a = iAstro_Primary_T;
model_b = iAstro_Primary_T_INPP5A;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_T_INPP5A.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_T_INPP5A = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_T_INPP5A

%% PLA2G6 (8398)
[~, GeneToRxnList, iAstro_Primary_T_PLA2G6] = downregulateGene(iAstro_Primary_T, '8398');
model_a = iAstro_Primary_T;
model_b = iAstro_Primary_T_PLA2G6;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_T_PLA2G6.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_T_PLA2G6 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_T_PLA2G6

%% GBA (2629)
[~, GeneToRxnList, iAstro_Primary_T_GBA] = downregulateGene(iAstro_Primary_T, '2629');
model_a = iAstro_Primary_T;
model_b = iAstro_Primary_T_GBA;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_T_GBA.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_T_GBA = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_T_GBA

%% PLCB4 (5332)
[~, GeneToRxnList, iAstro_Primary_T_PLCB4] = downregulateGene(iAstro_Primary_T, '5332');
model_a = iAstro_Primary_T;
model_b = iAstro_Primary_T_PLCB4;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_T_PLCB4.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_T_PLCB4 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_T_PLCB4

%% iAstro_iPS_Ctrl_T
%% CKMT2 (1160)
[~, GeneToRxnList, iAstro_iPS_Ctrl_T_CKMT2] = downregulateGene(iAstro_iPS_Ctrl_T, '1160');
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_Ctrl_T_CKMT2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_T_CKMT2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_T_CKMT2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_T_CKMT2

%% INPP5A (3632)
[~, GeneToRxnList, iAstro_iPS_Ctrl_T_INPP5A] = downregulateGene(iAstro_iPS_Ctrl_T, '3632');
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_Ctrl_T_INPP5A;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_T_INPP5A.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_T_INPP5A = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_T_INPP5A

%% PLA2G6 (8398)
[~, GeneToRxnList, iAstro_iPS_Ctrl_T_PLA2G6] = downregulateGene(iAstro_iPS_Ctrl_T, '8398');
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_Ctrl_T_PLA2G6;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_T_PLA2G6.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_T_PLA2G6 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_T_PLA2G6

%% GBA (2629)
[~, GeneToRxnList, iAstro_iPS_Ctrl_T_GBA] = downregulateGene(iAstro_iPS_Ctrl_T, '2629');
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_Ctrl_T_GBA;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_T_GBA.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_T_GBA = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_T_GBA

%% PLCB4 (5332)
[~, GeneToRxnList, iAstro_iPS_Ctrl_T_PLCB4] = downregulateGene(iAstro_iPS_Ctrl_T, '5332');
model_a = iAstro_iPS_Ctrl_T;
model_b = iAstro_iPS_Ctrl_T_PLCB4;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_T_PLCB4.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_T_PLCB4 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_T_PLCB4

%% save data to mat file
save('iAstro_FluxDiff_GeneList_T.mat'); 

%%
toc;
