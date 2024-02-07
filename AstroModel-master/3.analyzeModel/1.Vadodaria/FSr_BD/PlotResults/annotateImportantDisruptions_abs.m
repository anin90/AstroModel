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

%% MODEL_ABS

    % bd_lumped
    rxnList = importdata('bd_tbl_abs/bd_lumped.csv'); 
    model = iAstro_iPS_Ctrl_TP_abs;
    FSrTbl = FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_TP_abs;
    [bd_lumped] = annotateRxnList(model, rxnList, FSrTbl);
 
    % bd_responder
    rxnList = importdata('bd_tbl_abs/bd_responder.csv'); 
    model = iAstro_iPS_Ctrl_TP_abs;
    FSrTbl = FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_R_TP_abs;
    [bd_responder] = annotateRxnList(model, rxnList, FSrTbl);
    
    % bd_nonresponder
    rxnList = importdata('bd_tbl_abs/bd_nonresponder.csv'); 
    model = iAstro_iPS_Ctrl_TP_abs;
    FSrTbl = FluxDiff_iPS_Ctrl_TP_abs_vs_iPS_BD_NR_TP_abs;
    [bd_nonresponder] = annotateRxnList(model, rxnList, FSrTbl);
    
    writetable(bd_lumped, 'bd_tbl_abs/bd_lumped.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    writetable(bd_responder, 'bd_tbl_abs/bd_responder.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    writetable(bd_nonresponder, 'bd_tbl_abs/bd_nonresponder.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except bd_lumped bd_responder bd_nonresponder

%%
save('bd_tbl_abs/bd_tbl_abs.mat');
