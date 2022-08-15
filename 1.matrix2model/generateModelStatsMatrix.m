%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/')

%% initialize variable 'k'
k=1;

%% Primary_Ctrl

    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelID{k,1} = 'v12_TP_abs_1';
    k = k+1;
    
%% Final Table    
modelStatsMat = [modelID, modelStatsMat];
modelStatsMat.Properties.VariableNames{1} = 'model_ID';

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except modelStatsMat
