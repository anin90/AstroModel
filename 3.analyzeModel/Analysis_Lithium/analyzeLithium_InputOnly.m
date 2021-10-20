%% Load Data
tic;
% Load Primary_Models (Li+)
FileName   = 'iAstro_Lithium.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% iAstro_Primary_TP
%% Untreated vs Lithium (Akkouh)
[FluxDiff] = iAstro_Primary_TP_Akkouh_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'Primary_TP_Akkouh_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% Untreated vs Lithium (GSE66276)
[FluxDiff] = iAstro_Primary_TP_GSE66276_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'Primary_TP_GSE66276_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% Untreated vs Lithium (GSE132397)
[FluxDiff] = iAstro_Primary_TP_GSE132397_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'Primary_TP_GSE132397_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% iAstro_iPS_Ctrl_TP
%% Untreated vs Lithium (Akkouh)
[FluxDiff] = iAstro_iPS_Ctrl_TP_Akkouh_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'iPS_Ctrl_TP_Akkouh_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% Untreated vs Lithium (GSE66276)
[FluxDiff] = iAstro_iPS_Ctrl_TP_GSE66276_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'iPS_Ctrl_TP_GSE66276_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% Untreated vs Lithium (GSE132397)
[FluxDiff] = iAstro_iPS_Ctrl_TP_GSE132397_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'iPS_Ctrl_TP_GSE132397_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% iAstro_iPS_BD_TP
%% Untreated vs Lithium (Akkouh)
[FluxDiff] = iAstro_iPS_BD_TP_Akkouh_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'iPS_BD_TP_Akkouh_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% Untreated vs Lithium (GSE66276)
[FluxDiff] = iAstro_iPS_BD_TP_GSE66276_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'iPS_BD_TP_GSE66276_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%% Untreated vs Lithium (GSE132397)
[FluxDiff] = iAstro_iPS_BD_TP_GSE132397_Li;
dat = FluxDiff;
m_in = vertcat(union(dat.rxnList_down, dat.rxnList_up));
T_in = table(m_in);
writetable(T_in,'iPS_BD_TP_GSE132397_Input.csv','WriteVariableNames', true, 'Delimiter','\t');
clear FluxDiff dat m_in T_in

%%
clear all

%%
toc;