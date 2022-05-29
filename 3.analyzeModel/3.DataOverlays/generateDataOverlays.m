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

%% Li_analysis:

    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_abs, 'Lanz_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Ctrl_vs_Lanz.csv','WriteVariableNames', true, 'Delimiter','\t');
    
    dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_abs, 'Rivera_Li.mat');
    Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
    writetable(Tbl,'AnalysisFiles/Ctrl_vs_Rivera.csv','WriteVariableNames', true, 'Delimiter','\t');    
    
%     dat = mapDiffExpToRxns(iAstro_iPS_Ctrl_TP_abs, 'Akkouh_Li.mat');
%     Tbl = table(vertcat(union(dat.rxnList_down, dat.rxnList_up)));
%     writetable(Tbl,'AnalysisFiles/Ctrl_vs_Akkouh.csv','WriteVariableNames', true, 'Delimiter','\t');   
    
%%
clear all

%%
toc;