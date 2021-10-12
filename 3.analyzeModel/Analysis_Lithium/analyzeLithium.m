%% Load Data
tic;
% Load Primary_Models (Li+)
FileName   = 'iAstro_Lithium.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% iAstro_Primary_TP
%% Untreated vs Lithium
[FluxDiff] = iAstro_Primary_TP_Akkouh_Li;
dat = FluxDiff;
m_out = vertcat(union(dat.fluxGoneHigh.FSr_significant_H, dat.fluxGoneLow.FSr_significant_L));
T_out = table(m_out);
m_in = vertcat(union(dat.rxnList_down, union (dat.rxnList_up_f, dat.rxnList_up_b)));
T_in = table(m_in);
writetable(T_out,'iAstro_Primary_TP_Untreated_vs_Lithium_Output.csv','WriteVariableNames',0);
writetable(T_in,'iAstro_Primary_TP_Untreated_vs_Lithium_Input.csv','WriteVariableNames',0);
clear FluxDiff dat m_out T_out m_in T_in

%%
clearvars -except iAstro_Primary_TP_Akkouh_Li

%%
toc;