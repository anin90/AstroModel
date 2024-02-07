tStart = tic;

%% Load Data
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays')

%% model_abs:

    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_abs, 'Lanz_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Lanz_vs_Ctrl_abs.csv','WriteVariableNames', true, 'Delimiter','\t');
    
    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_abs, 'Rivera_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Rivera_vs_Ctrl_abs.csv','WriteVariableNames', true, 'Delimiter','\t');    
    
    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_abs, 'Akkouh_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Akkouh_vs_Ctrl_abs.csv','WriteVariableNames', true, 'Delimiter','\t');   

%% model_norm_t1:

    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_norm_t1, 'Lanz_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Lanz_vs_Ctrl_norm_t1.csv','WriteVariableNames', true, 'Delimiter','\t');
    
    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_norm_t1, 'Rivera_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Rivera_vs_Ctrl_norm_t1.csv','WriteVariableNames', true, 'Delimiter','\t');    
    
    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_norm_t1, 'Akkouh_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Akkouh_vs_Ctrl_norm_t1.csv','WriteVariableNames', true, 'Delimiter','\t');   

%% model_norm_t2:

    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_norm_t2, 'Lanz_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Lanz_vs_Ctrl_norm_t2.csv','WriteVariableNames', true, 'Delimiter','\t');
    
    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_norm_t2, 'Rivera_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Rivera_vs_Ctrl_norm_t2.csv','WriteVariableNames', true, 'Delimiter','\t');    
    
    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_norm_t2, 'Akkouh_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Akkouh_vs_Ctrl_norm_t2.csv','WriteVariableNames', true, 'Delimiter','\t');   

%%
clearvars

%%
toc;