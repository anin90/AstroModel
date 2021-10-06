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
clearvars -except iAstro_Primary_TP iAstro_iPS_Ctrl_TP

%% add path to run 'compareCbModel', 'downregulateGene'
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')

%% iAstro_Primary_TP
%% CKMT2 (1160)
[~, GeneToRxnList, iAstro_Primary_TP_CKMT2] = downregulateGene(iAstro_Primary_TP, '1160');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_CKMT2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_CKMT2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_CKMT2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_CKMT2 

%% INPP5A (3632)
[~, GeneToRxnList, iAstro_Primary_TP_INPP5A] = downregulateGene(iAstro_Primary_TP, '3632');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_INPP5A;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_INPP5A.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_INPP5A = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_INPP5A

%% PLA2G6 (8398)
[~, GeneToRxnList, iAstro_Primary_TP_PLA2G6] = downregulateGene(iAstro_Primary_TP, '8398');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_PLA2G6;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_PLA2G6.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_PLA2G6 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_PLA2G6

%% GBA (2629)
[~, GeneToRxnList, iAstro_Primary_TP_GBA] = downregulateGene(iAstro_Primary_TP, '2629');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_GBA;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_GBA.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_GBA = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_GBA

%% PLCB4 (5332)
[~, GeneToRxnList, iAstro_Primary_TP_PLCB4] = downregulateGene(iAstro_Primary_TP, '5332');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_PLCB4;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_PLCB4.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_PLCB4 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_PLCB4

%% iAstro_iPS_Ctrl_TP
%% CKMT2 (1160)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_CKMT2] = downregulateGene(iAstro_iPS_Ctrl_TP, '1160');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_CKMT2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_CKMT2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_CKMT2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_CKMT2

%% INPP5A (3632)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_INPP5A] = downregulateGene(iAstro_iPS_Ctrl_TP, '3632');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_INPP5A;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_INPP5A.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_INPP5A = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_INPP5A

%% PLA2G6 (8398)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_PLA2G6] = downregulateGene(iAstro_iPS_Ctrl_TP, '8398');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_PLA2G6;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_PLA2G6.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_PLA2G6 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_PLA2G6

%% GBA (2629)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_GBA] = downregulateGene(iAstro_iPS_Ctrl_TP, '2629');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_GBA;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_GBA.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_GBA = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_GBA

%% PLCB4 (5332)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_PLCB4] = downregulateGene(iAstro_iPS_Ctrl_TP, '5332');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_PLCB4;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_PLCB4.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_PLCB4 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_PLCB4

%% save data to mat file
save('iAstro_FluxDiff_GeneList_TP.mat'); 

%%
toc;