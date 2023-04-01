%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/')

%% initialize variable 'k'
k=1;

%% Zhang et al.
%% Primary

    % Zhang_Primary_abs_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;
    
    % Zhang_Primary_abs_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'GIMME';
    k = k+1;    
    
    % Zhang_Primary_abs_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'MBA';
    k = k+1;      

    % Zhang_Primary_abs_FastCore
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
    model = data.FastCore_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'FastCore';
    k = k+1;    
    
    % Zhang_Primary_norm_t1_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;
    
    % Zhang_Primary_norm_t1_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'GIMME';
    k = k+1;    
    
    % Zhang_Primary_norm_t1_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'MBA';
    k = k+1;      

    % Zhang_Primary_norm_t1_FastCore
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
    model = data.FastCore_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'FastCore';
    k = k+1;    

    % Zhang_Primary_norm_t2_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/matrix2models_norm_t2_v12.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;
    
    % Zhang_Primary_norm_t2_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/matrix2models_norm_t2_v12.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'GIMME';
    k = k+1;    
    
    % Zhang_Primary_norm_t2_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/GSE73721_HMA_CTX.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/matrix2models_norm_t2_v12.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Zhang';
    modelID{k,2} = 'Primary';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'MBA';
    k = k+1;      

    % Primary_norm_t2_4
    % FastCore InFeasible
    % k = k+1;   

%% Vadodaria et al.
%% iPS_Ctrl

    % Vadodaria_Ctrl_abs_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_Ctrl_abs_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_Ctrl_abs_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_Ctrl_abs_4
    % FastCore InFeasible
    % k = k+1;  
    
    % Vadodaria_Ctrl_norm_t1_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_Ctrl_norm_t1_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_Ctrl_norm_t1_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_Ctrl_norm_t1_4
    % FastCore InFeasible
    % k = k+1;  
    
    % Vadodaria_Ctrl_norm_t2_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_Ctrl_norm_t2_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_Ctrl_norm_t2_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/Vadodaria_Control_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_Ctrl_norm_t2_4
    % FastCore InFeasible
    % k = k+1;  

%% iPS_BD

    % Vadodaria_BD_abs_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_abs_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_BD_abs_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_BD_abs_4
    % FastCore InFeasible
    % k = k+1;  

    % Vadodaria_BD_norm_t1_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_norm_t1_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_BD_norm_t1_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_BD_norm_t1_4
    % FastCore InFeasible
    % k = k+1;      

    % Vadodaria_BD_norm_t2_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_norm_t2_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_BD_norm_t2_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/Vadodaria_BD_Untreated.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'MBA';
    k = k+1;

%     % iPS_BD_norm_t2_4
%     % FastCore InFeasible
%     % k = k+1;              
    
%% iPS_BD_R

    % Vadodaria_BD_R_abs_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_R_abs_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_BD_R_abs_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_BD_R_abs_4
    % FastCore InFeasible
    % k = k+1;  

    % Vadodaria_BD_R_norm_t1_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_R_norm_t1_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_BD_R_norm_t1_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_BD_R_norm_t1_4
    % FastCore InFeasible
    % k = k+1;  
    
    % Vadodaria_BD_R_norm_t2_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_R_norm_t2_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_BD_R_norm_t2_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_Responder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_BD_R_norm_t2_4
    % FastCore InFeasible
    % k = k+1;          

%% iPS_BD_NR

    % Vadodaria_BD_NR_abs_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_NR_abs_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'GIMME';
    k = k+1;

    % Vadodaria_BD_NR_abs_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/matrix2models_abs_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'MBA';
    k = k+1;

    % iPS_BD_NR_abs_4
    % FastCore InFeasible
    % k = k+1;  

    % Vadodaria_BD_NR_norm_t1_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Vadodaria_BD_NR_norm_t1_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'GIMME'; 
    k = k+1;

    % Vadodaria_BD_NR_norm_t1_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'MBA'; 
    k = k+1;

    % iPS_BD_NR_norm_t1_4
    % FastCore InFeasible
    % k = k+1;  

    % Vadodaria_BD_NR_norm_t2_iMAT
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.iMAT_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT'; 
    k = k+1;

    % Vadodaria_BD_NR_norm_t2_GIMME
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.GIMME_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'GIMME'; 
    k = k+1;

    % Vadodaria_BD_NR_norm_t2_MBA
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_NonResponder.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
    model = data.MBA_model_TP;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Vadodaria';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'MBA';    
    k = k+1;

    % iPS_BD_NR_norm_t2_4
    % FastCore InFeasible
    % k = k+1;  
    
%% Koskuvi et al.
%% Control

    % Koskuvi_Ctrl_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/abs/Koskuvi_Control.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Koskuvi_Ctrl_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Koskuvi_Ctrl_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t2/Koskuvi_Control.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'Ctrl';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% HT

    % Koskuvi_HT_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/abs/Koskuvi_HT.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'HT';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Koskuvi_HT_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/Koskuvi_HT.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'HT';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Koskuvi_HT_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t2/Koskuvi_HT.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'HT';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% ST

    % Koskuvi_ST_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/abs/Koskuvi_ST.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'ST';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Koskuvi_ST_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/Koskuvi_ST.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'ST';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Koskuvi_ST_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t2/Koskuvi_ST.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Koskuvi';
    modelID{k,2} = 'ST';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% Biju et al.
%% Ctrl_NIH1

    % Biju_Ctrl_NIH1_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/1.Ctrl_NIH1/v1/abs/Biju_Ctrl_NIH1.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/1.Ctrl_NIH1/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'Ctrl_NIH1';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_Ctrl_NIH1_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/1.Ctrl_NIH1/v1/norm_t1/Biju_Ctrl_NIH1.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/1.Ctrl_NIH1/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'Ctrl_NIH1';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_Ctrl_NIH1_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/1.Ctrl_NIH1/v1/norm_t2/Biju_Ctrl_NIH1.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/1.Ctrl_NIH1/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'Ctrl_NIH1';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% Ctrl_Fam

    % Biju_Ctrl_Fam_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/2.Ctrl_Fam/v1/abs/Biju_Ctrl_Fam.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/2.Ctrl_Fam/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'Ctrl_Fam';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_Ctrl_Fam_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/2.Ctrl_Fam/v1/norm_t1/Biju_Ctrl_Fam.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/2.Ctrl_Fam/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'Ctrl_Fam';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_Ctrl_Fam_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/2.Ctrl_Fam/v1/norm_t2/Biju_Ctrl_Fam.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/2.Ctrl_Fam/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'Ctrl_Fam';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% BD

    % Biju_BD_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/3.BD/v1/abs/Biju_BD.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/3.BD/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_BD_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/3.BD/v1/norm_t1/Biju_BD.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/3.BD/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_BD_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/3.BD/v1/norm_t2/Biju_BD.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/3.BD/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% BD_R

    % Biju_BD_R_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/4.BD_R/v1/abs/Biju_BD_R.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/4.BD_R/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_BD_R_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/4.BD_R/v1/norm_t1/Biju_BD_R.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/4.BD_R/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_BD_R_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/4.BD_R/v1/norm_t2/Biju_BD_R.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/4.BD_R/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD_R';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;

%% BD_NR

    % Biju_BD_NR_abs
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/5.BD_NR/v1/abs/Biju_BD_NR.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/5.BD_NR/v1/abs/matrix2models_abs_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'abs';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_BD_NR_norm_t1
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/5.BD_NR/v1/norm_t1/Biju_BD_NR.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/5.BD_NR/v1/norm_t1/matrix2models_norm_t1_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t1';
    modelID{k,4} = 'iMAT';
    k = k+1;

    % Biju_BD_NR_norm_t2
    expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/5.BD_NR/v1/norm_t2/Biju_BD_NR.mat';
    data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/4.Biju/1.Veh/5.BD_NR/v1/norm_t2/matrix2models_norm_t2_v1.mat');
    model = data;
    modelStatsMat = printModelStats(model, expMat);
    modelStatsMatSol(k,:) = modelStatsMat;
    modelID{k,1} = 'Biju';
    modelID{k,2} = 'BD_NR';
    modelID{k,3} = 'norm_t2';
    modelID{k,4} = 'iMAT';
    k = k+1;
    
%% Final Table    
modelStatsMatSol = [modelID, modelStatsMatSol];
modelStatsMatSol.Properties.VariableNames{1} = 'Dataset';
modelStatsMatSol.Properties.VariableNames{2} = 'Phenotype';
modelStatsMatSol.Properties.VariableNames{3} = 'ExpThreshold';
modelStatsMatSol.Properties.VariableNames{4} = 'MEM';
writetable(modelStatsMatSol, 'modelStatsMatFiles/modelStatsMatSol.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except modelStatsMatSol
