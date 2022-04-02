tStart = tic;

%% Load Data
% Load iAstro_Models
% FileName   = 'iAstro_Models.mat';
% FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
% File = fullfile(FolderName, FileName);
% load(File);
% clear File FileName FolderName

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_MTA_TestRuns/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Omics/Transcriptomics/RNA_Seq/2.Vadodaria_2021/Raw_data/6.analysis/MTA_input_bpad_untreated/')

%% DATASET: SOURCE to TAREGT
%% Vadodaria: Ctrl_abs to BD_abs
%     run import_MTA_Ctrl_to_BD_Lumped.m                      % DEGs (SOURCE to TAREGT)
%     model = iAstro_iPS_Ctrl_TP_abs;                         % Model (SOURCE)
%     [Ctrl_to_BD] = generateMTAscores(model, diff_exprs);    % Run MTA

%% Vadodaria: BD_abs to Ctrl_abs
%     run import_MTA_BD_Lumped_to_Ctrl.m                      % DEGs (SOURCE to TAREGT)
%     model = iAstro_iPS_BD_TP_abs;                           % Model (SOURCE)
%     [BD_to_Ctrl] = generateMTAscores(model, diff_exprs);    % Run MTA
    
%% Vadodaria: BD_NR_abs to BD_R_abs

%% GSE32323: Healthy to Cancer
    load('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_MTA_TestRuns/ExampleFiles_GSE32323/model.mat');
    load('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_MTA_TestRuns/ExampleFiles_GSE32323/GSE32323.mat');
    model = defineHumanMediaDMEM(model, flag);
    exprs.source = GSE32323.healthy_tab;
    exprs.target = GSE32323.cancer_tab;
    exprs.smean = mean(exprs.source, 2);
    exprs.tmean = mean(exprs.target, 2);
    [mta_tbl_filt, del_rxnID, mta_score] = generateMTAscores_b(model, exprs);
    [del_rxnID_FUT9] = intersect(findRxnsActiveWithGenes(model, '10690.1')', del_rxnID);
    [del_rxnID_AKR7A2] = intersect(findRxnsActiveWithGenes(model, '8574.1')', del_rxnID);
    
%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except del_rxnID mta_score del_rxnID_FUT9 del_rxnID_AKR7A2
