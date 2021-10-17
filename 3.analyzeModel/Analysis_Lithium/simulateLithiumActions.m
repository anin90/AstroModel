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

%% add path to run 'constrainCbModel'
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')

%% Li_analysis: TP_models 

tic; [iAstro_Primary_TP_Akkouh_Li] = constrainCbModel(iAstro_Primary_TP, 'Akkouh_Redone_Li.mat'); toc;
tic; [iAstro_Primary_TP_GSE66276_Li] = constrainCbModel(iAstro_Primary_TP, 'GSE66276_Li.mat'); toc;
tic; [iAstro_Primary_TP_GSE132397_Li] = constrainCbModel(iAstro_Primary_TP, 'GSE132397_Li.mat'); toc;

tic; [iAstro_iPS_Ctrl_TP_Akkouh_Li] = constrainCbModel(iAstro_iPS_Ctrl_TP, 'Akkouh_Redone_Li.mat'); toc;
tic; [iAstro_iPS_Ctrl_TP_GSE66276_Li] = constrainCbModel(iAstro_iPS_Ctrl_TP, 'GSE66276_Li.mat'); toc;
tic; [iAstro_iPS_Ctrl_TP_GSE132397_Li] = constrainCbModel(iAstro_iPS_Ctrl_TP, 'GSE132397_Li.mat'); toc;

tic; [iAstro_iPS_BD_TP_Akkouh_Li] = constrainCbModel(iAstro_iPS_BD_TP, 'Akkouh_Redone_Li.mat'); toc;
tic; [iAstro_iPS_BD_TP_GSE66276_Li] = constrainCbModel(iAstro_iPS_BD_TP, 'GSE66276_Li.mat'); toc;
tic; [iAstro_iPS_BD_TP_GSE132397_Li] = constrainCbModel(iAstro_iPS_BD_TP, 'GSE132397_Li.mat'); toc;

%%
tic; save('iAstro_Lithium.mat'); toc;

%%
clear all
%%
toc;