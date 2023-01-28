%%
tStart = tic;

%% Load Data
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/')

%% initialize variable 'k'
k=1;

%% Zhang et al.
%% Primary

    % Zhang_Primary_abs_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';    
    model = iAstro_Primary_TP_abs;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Zhang_Primary_norm_t1_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
    model = iAstro_Primary_TP_norm_t1;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;
    
    % Zhang_Primary_norm_t2_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/GSE73721_HMA_CTX.mat';
    model = iAstro_Primary_TP_norm_t2;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% Vadodaria et al.
%% iPS_Ctrl

    % Vadodaria_Ctrl_abs_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/Vadodaria_Control_Untreated.mat';
    model = iAstro_iPS_Ctrl_TP_abs;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_Ctrl_norm_t1_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat';
    model = iAstro_iPS_Ctrl_TP_norm_t1;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;
    
    % Vadodaria_Ctrl_norm_t2_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/Vadodaria_Control_Untreated.mat';
    model = iAstro_iPS_Ctrl_TP_norm_t2;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;
    
%% Final Table    
modelStatsMatSol = [modelID, modelStatsMatSol];
modelStatsMatSol.Properties.VariableNames{1} = 'Dataset';
modelStatsMatSol.Properties.VariableNames{2} = 'Phenotype';
modelStatsMatSol.Properties.VariableNames{3} = 'ExpThreshold';
modelStatsMatSol.Properties.VariableNames{4} = 'MEM';
writetable(modelStatsMatSol, 'modelStatsMatFilesFinal/modelStatsMatSolFinal.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except modelStatsMatSol