%% Load Data
tic;
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% Retain Primary and iPS_Ctrl models
clearvars -except iAstro_Primary_TP_abs iAstro_iPS_Ctrl_TP_abs iAstro_iPS_SCZ_Ctrl_TP_abs 

%% add path to run 'compareCbModel', 'downregulateGeneStatic', 'downregulateGeneDynamic'
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')

%% iAstro_Primary_TP
%% CKMT2 (1160)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_CKMT2] = downregulateGeneStatic(iAstro_Primary_TP_abs, '1160');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_CKMT2;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_CKMT2.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_CKMT2 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_CKMT2 

%% INPP5A (3632)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_INPP5A] = downregulateGeneStatic(iAstro_Primary_TP_abs, '3632');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_INPP5A;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_INPP5A.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_INPP5A = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_INPP5A

%% PLA2G6 (8398)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_PLA2G6] = downregulateGeneStatic(iAstro_Primary_TP_abs, '8398');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_PLA2G6;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_PLA2G6.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_PLA2G6 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_PLA2G6

%% GBA (2629)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_GBA] = downregulateGeneStatic(iAstro_Primary_TP_abs, '2629');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_GBA;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_GBA.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_GBA = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_GBA

%% PLCB4 (5332)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_PLCB4] = downregulateGeneStatic(iAstro_Primary_TP_abs, '5332');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_PLCB4;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_PLCB4.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_PLCB4 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_PLCB4

%% NMNAT1 (64802)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_NMNAT1] = downregulateGeneStatic(iAstro_Primary_TP_abs, '64802');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_NMNAT1;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_NMNAT1.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_NMNAT1 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_NMNAT1

%% ASPA (443)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_ASPA] = downregulateGeneStatic(iAstro_Primary_TP_abs, '443');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_ASPA;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_ASPA.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_ASPA = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_ASPA

%% SLC5A8 (160728)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_SLC5A8] = downregulateGeneStatic(iAstro_Primary_TP_abs, '160728');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_SLC5A8;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_SLC5A8.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_SLC5A8 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_SLC5A8

%% PLA2G12A (81579)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_PLA2G12A] = downregulateGeneStatic(iAstro_Primary_TP_abs, '81579');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_PLA2G12A;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_PLA2G12A.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_PLA2G12A = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_PLA2G12A

%% SLC29A2 (3177)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_SLC29A2] = downregulateGeneStatic(iAstro_Primary_TP_abs, '3177');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_SLC29A2;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_SLC29A2.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_SLC29A2 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_SLC29A2

%% COASY (80347)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_COASY] = downregulateGeneStatic(iAstro_Primary_TP_abs, '80347');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_COASY;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_COASY.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_COASY = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_COASY

%% SYNJ2 (8871)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_SYNJ2] = downregulateGeneStatic(iAstro_Primary_TP_abs, '8871');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_SYNJ2;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_SYNJ2.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_SYNJ2 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_SYNJ2

%% PYGL (5836)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_PYGL] = downregulateGeneStatic(iAstro_Primary_TP_abs, '5836');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_PYGL;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_PYGL.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_PYGL = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_PYGL

%% NEU1 (4758)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_NEU1] = downregulateGeneStatic(iAstro_Primary_TP_abs, '4758');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_NEU1;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_NEU1.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_NEU1 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_NEU1

%% UGT2A1 (10941)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_UGT2A1] = downregulateGeneStatic(iAstro_Primary_TP_abs, '10941');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_UGT2A1;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_UGT2A1.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_UGT2A1 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_UGT2A1

%% SLC22A9 (114571)
    [~, GeneToRxnList, iAstro_Primary_TP_abs_SLC22A9] = downregulateGeneStatic(iAstro_Primary_TP_abs, '114571');
    model_a = iAstro_Primary_TP_abs;
    model_b = iAstro_Primary_TP_abs_SLC22A9;
    [FluxDiff] = compareCbModel(model_a, model_b);
    dat = FluxDiff;
    m = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
    T = table(m);
    writetable(T,'FSR_Primary_TP_abs_SLC22A9.csv','WriteVariableNames',0);
    FluxDiff.GeneToRxnList = GeneToRxnList;
    FluxDiff_Primary_TP_abs_SLC22A9 = FluxDiff;
    clear model_a model_b dat m T FluxDiff GeneToRxnList iAstro_Primary_TP_abs_SLC22A9

%% iAstro_iPS_Ctrl_TP


%% save data to mat file
save('iAstro_FluxDiff_GeneLoFAnalysis.mat');

%%
toc;
