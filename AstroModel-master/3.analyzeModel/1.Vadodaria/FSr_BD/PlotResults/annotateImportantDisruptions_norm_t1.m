%% Load Data
tic;
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName
% Load FluxSpanData
FileName   = 'iAstro_FluxDiff_BD.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% MODEL_NORM_T1

    % bd_lumped
    rxnList = importdata('bd_tbl_norm_t1/bd_lumped.csv'); 
    model = iAstro_iPS_Ctrl_TP_norm_t1;
    FSrTbl = FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_TP_norm_t1;
    [bd_lumped] = annotateRxnList(model, rxnList, FSrTbl);
 
    % bd_responder
    rxnList = importdata('bd_tbl_norm_t1/bd_responder.csv'); 
    model = iAstro_iPS_Ctrl_TP_norm_t1;
    FSrTbl = FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_R_TP_norm_t1;
    [bd_responder] = annotateRxnList(model, rxnList, FSrTbl);
    
    % bd_nonresponder
    rxnList = importdata('bd_tbl_norm_t1/bd_nonresponder.csv'); 
    model = iAstro_iPS_Ctrl_TP_norm_t1;
    FSrTbl = FluxDiff_iPS_Ctrl_TP_norm_t1_vs_iPS_BD_NR_TP_norm_t1;
    [bd_nonresponder] = annotateRxnList(model, rxnList, FSrTbl);
    
    writetable(bd_lumped, 'bd_tbl_norm_t1/bd_lumped.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    writetable(bd_responder, 'bd_tbl_norm_t1/bd_responder.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    writetable(bd_nonresponder, 'bd_tbl_norm_t1/bd_nonresponder.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except bd_lumped bd_responder bd_nonresponder

%%
save('bd_tbl_norm_t1/bd_tbl_norm_t1.mat');
