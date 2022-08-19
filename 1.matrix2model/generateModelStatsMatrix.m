%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/')

%% initialize variable 'k'
k=1;

%% Primary

%     % Primary_abs_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_abs_1';
%     k = k+1;
%     
%     % Primary_abs_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_abs_2';
%     k = k+1;    
%     
%     % Primary_abs_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_abs_3';
%     k = k+1;      
% 
%     % Primary_abs_4
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/abs/matrix2models_abs_v12.mat');
%     model = data.FastCore_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_abs_4';
%     k = k+1;    
%     
%     % Primary_norm_t1_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_norm_t1_1';
%     k = k+1;
%     
%     % Primary_norm_t1_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_norm_t1_2';
%     k = k+1;    
%     
%     % Primary_norm_t1_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_norm_t1_3';
%     k = k+1;      
% 
%     % Primary_norm_t1_4
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/matrix2models_norm_t1_v12.mat');
%     model = data.FastCore_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_norm_t1_4';
%     k = k+1;    
% 
%     % Primary_norm_t2_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/matrix2models_norm_t2_v12.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_norm_t2_1';
%     k = k+1;
%     
%     % Primary_norm_t2_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/matrix2models_norm_t2_v12.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_norm_t2_2';
%     k = k+1;    
%     
%     % Primary_norm_t2_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/GSE73721_HMA_CTX.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t2/matrix2models_norm_t2_v12.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'Primary_norm_t2_3';
%     k = k+1;      
% 
%     % Primary_norm_t2_4
%     % FastCore InFeasible
%     % k = k+1;   

%% iPS_Ctrl

%     % iPS_Ctrl_abs_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_abs_1';
%     k = k+1;
% 
%     % iPS_Ctrl_abs_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_abs_2';
%     k = k+1;
% 
%     % iPS_Ctrl_abs_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_abs_3';
%     k = k+1;
% 
%     % iPS_Ctrl_abs_4
%     % FastCore InFeasible
%     % k = k+1;  
%     
%     % iPS_Ctrl_norm_t1_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_norm_t1_1';
%     k = k+1;
% 
%     % iPS_Ctrl_norm_t1_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_norm_t1_2';
%     k = k+1;
% 
%     % iPS_Ctrl_norm_t1_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_norm_t1_3';
%     k = k+1;
% 
%     % iPS_Ctrl_norm_t1_4
%     % FastCore InFeasible
%     % k = k+1;  
%     
%     % iPS_Ctrl_norm_t2_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_norm_t2_1';
%     k = k+1;
% 
%     % iPS_Ctrl_norm_t2_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_norm_t2_2';
%     k = k+1;
% 
%     % iPS_Ctrl_norm_t2_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/Vadodaria_Control_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_Ctrl_norm_t2_3';
%     k = k+1;
% 
%     % iPS_Ctrl_norm_t2_4
%     % FastCore InFeasible
%     % k = k+1;  

%% iPS_BD

%     % iPS_BD_abs_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_abs_1';
%     k = k+1;
% 
%     % iPS_BD_abs_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_abs_2';
%     k = k+1;
% 
%     % iPS_BD_abs_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_abs_3';
%     k = k+1;
% 
%     % iPS_BD_abs_4
%     % FastCore InFeasible
%     % k = k+1;  
% 
%     % iPS_BD_norm_t1_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_norm_t1_1';
%     k = k+1;
% 
%     % iPS_BD_norm_t1_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_norm_t1_2';
%     k = k+1;
% 
%     % iPS_BD_norm_t1_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_norm_t1_3';
%     k = k+1;
% 
%     % iPS_BD_norm_t1_4
%     % FastCore InFeasible
%     % k = k+1;      
% 
%     % iPS_BD_norm_t2_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_norm_t2_1';
%     k = k+1;
% 
%     % iPS_BD_norm_t2_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_norm_t2_2';
%     k = k+1;
% 
%     % iPS_BD_norm_t2_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/Vadodaria_BD_Untreated.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_norm_t2_3';
%     k = k+1;
% 
%     % iPS_BD_norm_t2_4
%     % FastCore InFeasible
%     % k = k+1;              
    
%% iPS_BD_R

%     % iPS_BD_R_abs_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_abs_1';
%     k = k+1;
% 
%     % iPS_BD_R_abs_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_abs_2';
%     k = k+1;
% 
%     % iPS_BD_R_abs_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_abs_3';
%     k = k+1;
% 
%     % iPS_BD_R_abs_4
%     % FastCore InFeasible
%     % k = k+1;  
% 
%     % iPS_BD_R_norm_t1_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_norm_t1_1';
%     k = k+1;
% 
%     % iPS_BD_R_norm_t1_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_norm_t1_2';
%     k = k+1;
% 
%     % iPS_BD_R_norm_t1_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_norm_t1_3';
%     k = k+1;
% 
%     % iPS_BD_R_norm_t1_4
%     % FastCore InFeasible
%     % k = k+1;  
%     
%     % iPS_BD_R_norm_t2_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_norm_t2_1';
%     k = k+1;
% 
%     % iPS_BD_R_norm_t2_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_norm_t2_2';
%     k = k+1;
% 
%     % iPS_BD_R_norm_t2_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_Responder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_R_norm_t2_3';
%     k = k+1;
% 
%     % iPS_BD_R_norm_t2_4
%     % FastCore InFeasible
%     % k = k+1;          

%% iPS_BD_NR

%     % iPS_BD_NR_abs_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_abs_1';
%     k = k+1;
% 
%     % iPS_BD_NR_abs_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_abs_2';
%     k = k+1;
% 
%     % iPS_BD_NR_abs_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/abs/matrix2models_abs_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_abs_3';
%     k = k+1;
% 
%     % iPS_BD_NR_abs_4
%     % FastCore InFeasible
%     % k = k+1;  
% 
%     % iPS_BD_NR_norm_t1_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_norm_t1_1';
%     k = k+1;
% 
%     % iPS_BD_NR_norm_t1_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_norm_t1_2';
%     k = k+1;
% 
%     % iPS_BD_NR_norm_t1_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/matrix2models_norm_t1_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_norm_t1_3';
%     k = k+1;
% 
%     % iPS_BD_NR_norm_t1_4
%     % FastCore InFeasible
%     % k = k+1;  
% 
%     % iPS_BD_NR_norm_t2_1
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.iMAT_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_norm_t2_1';
%     k = k+1;
% 
%     % iPS_BD_NR_norm_t2_2
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.GIMME_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_norm_t2_2';
%     k = k+1;
% 
%     % iPS_BD_NR_norm_t2_3
%     expMat = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/Vadodaria_BD_Untreated_NonResponder.mat';
%     data = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t2/matrix2models_norm_t2_v3.mat');
%     model = data.MBA_model_TP;
%     modelStatsMat = printModelStats(model, expMat);
%     modelStatsMatSol(k,:) = modelStatsMat;
%     modelID{k,1} = 'iPS_BD_NR_norm_t2_3';
%     k = k+1;
% 
%     % iPS_BD_NR_norm_t2_4
%     % FastCore InFeasible
%     % k = k+1;  
    
%%


%% Final Table    
modelStatsMatSol = [modelID, modelStatsMatSol];
modelStatsMatSol.Properties.VariableNames{1} = 'model_ID';
writetable(modelStatsMatSol, 'modelStatsMatFiles/modelStatsMatSol.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except modelStatsMatSol
