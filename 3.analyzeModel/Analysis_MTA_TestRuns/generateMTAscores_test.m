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
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Omics/Transcriptomics/RNA_Seq/2.Vadodaria_2021/Raw_data/6.analysis/MTA_input_bpad_untreated/')

%% SOURCE to TAREGT
%% Ctrl_abs to BD_abs
    run import_MTA_Ctrl_to_BD_Lumped.m                      % DEGs (SOURCE to TAREGT)
    model = iAstro_iPS_Ctrl_TP_abs;                         % Model (SOURCE)
    [Ctrl_to_BD] = generateMTAscores(model, diff_exprs);    % Run MTA

%% BD_abs to Ctrl_abs
    run import_MTA_BD_Lumped_to_Ctrl.m                      % DEGs (SOURCE to TAREGT)
    model = iAstro_iPS_BD_TP_abs;                           % Model (SOURCE)
    [BD_to_Ctrl] = generateMTAscores(model, diff_exprs);    % Run MTA
    
%% BD_NR_abs to BD_R_abs

%% 

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except iAstro_iPS_Ctrl_TP_abs iAstro_iPS_BD_TP_abs Ctrl_to_BD BD_to_Ctrl
