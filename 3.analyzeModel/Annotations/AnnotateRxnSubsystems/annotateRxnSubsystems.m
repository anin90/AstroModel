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

% norm
writetable(table(iAstro_Primary_TP_norm.rxns, [iAstro_Primary_TP_norm.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_Primary_TP_norm.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_Ctrl_TP_norm.rxns, [iAstro_iPS_Ctrl_TP_norm.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_Ctrl_TP_norm.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_TP_norm.rxns, [iAstro_iPS_BD_TP_norm.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_TP_norm.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_R_TP_norm.rxns, [iAstro_iPS_BD_R_TP_norm.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_R_TP_norm.csv', 'WriteVariableNames', true, 'Delimiter', '\t');
writetable(table(iAstro_iPS_BD_NR_TP_norm.rxns, [iAstro_iPS_BD_NR_TP_norm.subSystems{:}]', 'VariableNames',{'Rxn','SubSystem'}), 'Rxns_iPS_BD_NR_TP_norm.csv', 'WriteVariableNames', true, 'Delimiter', '\t');

%%
