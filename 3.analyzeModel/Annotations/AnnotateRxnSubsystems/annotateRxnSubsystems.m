%% Load Data
tic;
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% write table
% abs
writetable(table(iAstro_Primary_TP_abs.rxns, [iAstro_Primary_TP_abs.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_Primary_TP_abs.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_Ctrl_TP_abs.rxns, [iAstro_iPS_Ctrl_TP_abs.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_Ctrl_TP_abs.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_TP_abs.rxns, [iAstro_iPS_BD_TP_abs.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_TP_abs.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_R_TP_abs.rxns, [iAstro_iPS_BD_R_TP_abs.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_R_TP_abs.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_NR_TP_abs.rxns, [iAstro_iPS_BD_NR_TP_abs.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_NR_TP_abs.csv', 'WriteVariableNames', true, 'Delimiter', '\t');

% norm_t1
writetable(table(iAstro_Primary_TP_norm_t1.rxns, [iAstro_Primary_TP_norm_t1.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_Primary_TP_norm_t1.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_Ctrl_TP_norm_t1.rxns, [iAstro_iPS_Ctrl_TP_norm_t1.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_Ctrl_TP_norm_t1.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_TP_norm_t1.rxns, [iAstro_iPS_BD_TP_norm_t1.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_TP_norm_t1.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_R_TP_norm_t1.rxns, [iAstro_iPS_BD_R_TP_norm_t1.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_R_TP_norm_t1.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_NR_TP_norm_t1.rxns, [iAstro_iPS_BD_NR_TP_norm_t1.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_NR_TP_norm_t1.csv', 'WriteVariableNames', true, 'Delimiter', '\t');

% norm_t2
writetable(table(iAstro_Primary_TP_norm_t2.rxns, [iAstro_Primary_TP_norm_t2.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_Primary_TP_norm_t2.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_Ctrl_TP_norm_t2.rxns, [iAstro_iPS_Ctrl_TP_norm_t2.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_Ctrl_TP_norm_t2.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_TP_norm_t2.rxns, [iAstro_iPS_BD_TP_norm_t2.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_TP_norm_t2.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_R_TP_norm_t2.rxns, [iAstro_iPS_BD_R_TP_norm_t2.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_R_TP_norm_t2.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_NR_TP_norm_t2.rxns, [iAstro_iPS_BD_NR_TP_norm_t2.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_NR_TP_norm_t2.csv', 'WriteVariableNames', true, 'Delimiter', '\t');

%%
