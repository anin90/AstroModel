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

%% NMNAT1 (64802)
[~, GeneToRxnList, iAstro_Primary_TP_NMNAT1] = downregulateGene(iAstro_Primary_TP, '64802');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_NMNAT1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_NMNAT1.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_NMNAT1 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_NMNAT1

%% ASPA (443)
[~, GeneToRxnList, iAstro_Primary_TP_ASPA] = downregulateGene(iAstro_Primary_TP, '443');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_ASPA;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_ASPA.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_ASPA = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_ASPA

%% SLC5A8 (160728)
[~, GeneToRxnList, iAstro_Primary_TP_SLC5A8] = downregulateGene(iAstro_Primary_TP, '160728');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_SLC5A8;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_SLC5A8.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_SLC5A8 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_SLC5A8

%% PLA2G12A (81579)
[~, GeneToRxnList, iAstro_Primary_TP_PLA2G12A] = downregulateGene(iAstro_Primary_TP, '81579');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_PLA2G12A;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_PLA2G12A.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_PLA2G12A = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_PLA2G12A

%% SLC29A2 (3177)
[~, GeneToRxnList, iAstro_Primary_TP_SLC29A2] = downregulateGene(iAstro_Primary_TP, '3177');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_SLC29A2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_SLC29A2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_SLC29A2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_SLC29A2

%% COASY (80347)
[~, GeneToRxnList, iAstro_Primary_TP_COASY] = downregulateGene(iAstro_Primary_TP, '80347');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_COASY;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_COASY.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_COASY = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_COASY

%% SYNJ2 (8871)
[~, GeneToRxnList, iAstro_Primary_TP_SYNJ2] = downregulateGene(iAstro_Primary_TP, '8871');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_SYNJ2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_SYNJ2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_SYNJ2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_SYNJ2

%% PYGL (5836)
[~, GeneToRxnList, iAstro_Primary_TP_PYGL] = downregulateGene(iAstro_Primary_TP, '5836');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_PYGL;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_PYGL.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_PYGL = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_PYGL

%% NEU1 (4758)
[~, GeneToRxnList, iAstro_Primary_TP_NEU1] = downregulateGene(iAstro_Primary_TP, '4758');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_NEU1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_NEU1.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_NEU1 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_NEU1

%% UGT2A1 (10941)
[~, GeneToRxnList, iAstro_Primary_TP_UGT2A1] = downregulateGene(iAstro_Primary_TP, '10941');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_UGT2A1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_UGT2A1.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_UGT2A1 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_UGT2A1

%% SLC22A9 (114571)
[~, GeneToRxnList, iAstro_Primary_TP_SLC22A9] = downregulateGene(iAstro_Primary_TP, '114571');
model_a = iAstro_Primary_TP;
model_b = iAstro_Primary_TP_SLC22A9;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_Primary_TP_SLC22A9.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_Primary_TP_SLC22A9 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_SLC22A9

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

%% NMNAT1 (64802)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_NMNAT1] = downregulateGene(iAstro_iPS_Ctrl_TP, '64802');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_NMNAT1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_NMNAT1.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_NMNAT1 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_NMNAT1

%% ASPA (443)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_ASPA] = downregulateGene(iAstro_iPS_Ctrl_TP, '443');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_ASPA;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_ASPA.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_ASPA = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_ASPA

%% SLC5A8 (160728)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_SLC5A8] = downregulateGene(iAstro_iPS_Ctrl_TP, '160728');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_SLC5A8;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_SLC5A8.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_SLC5A8 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_SLC5A8

%% PLA2G12A (81579)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_PLA2G12A] = downregulateGene(iAstro_iPS_Ctrl_TP, '81579');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_PLA2G12A;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_PLA2G12A.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_PLA2G12A = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_PLA2G12A

%% SLC29A2 (3177)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_SLC29A2] = downregulateGene(iAstro_iPS_Ctrl_TP, '3177');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_SLC29A2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_SLC29A2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_SLC29A2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_SLC29A2

%% COASY (80347)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_COASY] = downregulateGene(iAstro_iPS_Ctrl_TP, '80347');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_COASY;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_COASY.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_COASY = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_COASY

%% SYNJ2 (8871)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_SYNJ2] = downregulateGene(iAstro_iPS_Ctrl_TP, '8871');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_SYNJ2;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_SYNJ2.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_SYNJ2 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_SYNJ2

%% PYGL (5836)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_PYGL] = downregulateGene(iAstro_iPS_Ctrl_TP, '5836');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_PYGL;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_PYGL.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_PYGL = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_PYGL

%% NEU1 (4758)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_NEU1] = downregulateGene(iAstro_iPS_Ctrl_TP, '4758');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_NEU1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_NEU1.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_NEU1 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_NEU1

%% UGT2A1 (10941)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_UGT2A1] = downregulateGene(iAstro_iPS_Ctrl_TP, '10941');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_UGT2A1;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_UGT2A1.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_UGT2A1 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_UGT2A1

%% SLC22A9 (114571)
[~, GeneToRxnList, iAstro_iPS_Ctrl_TP_SLC22A9] = downregulateGene(iAstro_iPS_Ctrl_TP, '114571');
model_a = iAstro_iPS_Ctrl_TP;
model_b = iAstro_iPS_Ctrl_TP_SLC22A9;
[FluxDiff] = compareCbModel(model_a, model_b);
dat = FluxDiff;
m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T = table(m);
writetable(T,'FSR_iPS_Ctrl_TP_SLC22A9.csv','WriteVariableNames',0);
FluxDiff.GeneToRxnList = GeneToRxnList;
FluxDiff_iPS_Ctrl_TP_SLC22A9 = FluxDiff;
clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_iPS_Ctrl_TP_SLC22A9

%% save data to mat file
save('iAstro_FluxDiff_GeneList_TP.mat'); 

%%
toc;
