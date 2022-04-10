tStart = tic;

%% Load Data
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_MTA_TestRuns/')

%% PREPARE DATA

    % Vadodaria_Ctrl
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
    expression_Ctrl = prepareExpressionMatrix(data);
    
    % Vadodaria_BD
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat');
    expression_BD = prepareExpressionMatrix(data);
    
%% DATASET: SOURCE to TAREGT
%% Vadodaria: Ctrl_abs to BD_abs
    model = iAstro_iPS_Ctrl_TP_abs;
    model = prepareModel(model);


%% Vadodaria: BD_abs to Ctrl_abs

%% Vadodaria: BD_NR_abs to BD_R_abs

%% GSE32323: Healthy to Cancer
    % load('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_MTA_TestRuns/ExampleFiles_GSE32323/model.mat');
    % load('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_MTA_TestRuns/ExampleFiles_GSE32323/GSE32323.mat');
    % load('/home/anirudh/Downloads/papers_to_read/mta/MSB-13-956-s009/EV2/GSE32323/model.mat')
    % load('/home/anirudh/Downloads/papers_to_read/mta/MSB-13-956-s009/EV2/GSE32323/GSE32323.mat')
    % model = defineHumanMediaDMEM(model, flag);
    % exprs.source = GSE32323.healthy_tab;
    % exprs.target = GSE32323.cancer_tab;
    % exprs.smean = mean(exprs.source, 2);
    % exprs.tmean = mean(exprs.target, 2);
    % [mta_tbl_filt, del_rxnID, mta_score] = generateMTAscores_b(model, exprs);
    % [del_rxnID_FUT9] = intersect(findRxnsActiveWithGenes(model, '10690.1')', del_rxnID);
    % [del_rxnID_AKR7A2] = intersect(findRxnsActiveWithGenes(model, '8574.1')', del_rxnID);
    % [del_rxnID_PTEN] = intersect(findRxnsActiveWithGenes(model, '5728.1')', del_rxnID);
    
%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except model expression_Ctrl expression_BD
