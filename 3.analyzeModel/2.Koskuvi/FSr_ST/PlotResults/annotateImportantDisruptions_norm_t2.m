%% Load Data
tic;
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName
% Load FluxSpanData
FileName   = 'iAstro_FluxDiff_SCZ.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% MODEL_NORM_T2

    % st_rxns
    rxnList = importdata('st_tbl_norm_t2/st_rxns.csv'); 
    model = iAstro_iPS_SCZ_Ctrl_TP_norm_t2;
    FSrTbl = FluxDiff_iPS_SCZ_Ctrl_TP_norm_t2_vs_iPS_SCZ_ST_TP_norm_t2;
    [st_rxns] = annotateRxnList(model, rxnList, FSrTbl);
    
    writetable(st_rxns, 'st_tbl_norm_t2/st_rxns.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except st_rxns

%%
save('st_tbl_norm_t2/st_tbl_norm_t2.mat');
