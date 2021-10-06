%% Load Data
tic;
% Load iPS_Models
FileName   = 'iAstro_iPS_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
% Load Primary_Models
FileName   = 'iAstro_Primary_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% iAstro_Primary (T vs TP)
model_A = iAstro_Primary_T; model_B = iAstro_Primary_TP;

rxnList = setdiff(model_A.rxns, model_B.rxns);
[tf,loc] = ismember(model_A.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_A.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_a = table(rxnList, subSystem);
tbl_subsys_a = setdiff([model_A.subSystems{:}]',[model_B.subSystems{:}]');

rxnList = setdiff(model_B.rxns, model_A.rxns);
[tf,loc] = ismember(model_B.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_B.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_b = table(rxnList, subSystem);
tbl_subsys_b = setdiff([model_B.subSystems{:}]',[model_A.subSystems{:}]');

iAstro_Primary_T_vs_TP = struct('rxn_specific_T', {tbl_rxn_a}, 'subsys_specific_T', {tbl_subsys_a}, 'rxn_specific_TP', {tbl_rxn_b}, 'subsys_specific_TP', {tbl_subsys_b});

clear model_A model_B rxnList tf loc p idx subSystem tbl_rxn_a tbl_subsys_a tbl_rxn_b tbl_subsys_b

%% iAstro_iPS_Ctrl (T vs TP)
model_A = iAstro_iPS_Ctrl_T; model_B = iAstro_iPS_Ctrl_TP;

rxnList = setdiff(model_A.rxns, model_B.rxns);
[tf,loc] = ismember(model_A.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_A.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_a = table(rxnList, subSystem);
tbl_subsys_a = setdiff([model_A.subSystems{:}]',[model_B.subSystems{:}]');

rxnList = setdiff(model_B.rxns, model_A.rxns);
[tf,loc] = ismember(model_B.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_B.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_b = table(rxnList, subSystem);
tbl_subsys_b = setdiff([model_B.subSystems{:}]',[model_A.subSystems{:}]');

iAstro_iPS_Ctrl_T_vs_TP = struct('rxn_specific_T', {tbl_rxn_a}, 'subsys_specific_T', {tbl_subsys_a}, 'rxn_specific_TP', {tbl_rxn_b}, 'subsys_specific_TP', {tbl_subsys_b});

clear model_A model_B rxnList tf loc p idx subSystem tbl_rxn_a tbl_subsys_a tbl_rxn_b tbl_subsys_b

%% iAstro_iPS_BD (T vs TP)
model_A = iAstro_iPS_BD_T; model_B = iAstro_iPS_BD_TP;

rxnList = setdiff(model_A.rxns, model_B.rxns);
[tf,loc] = ismember(model_A.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_A.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_a = table(rxnList, subSystem);
tbl_subsys_a = setdiff([model_A.subSystems{:}]',[model_B.subSystems{:}]');

rxnList = setdiff(model_B.rxns, model_A.rxns);
[tf,loc] = ismember(model_B.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_B.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_b = table(rxnList, subSystem);
tbl_subsys_b = setdiff([model_B.subSystems{:}]',[model_A.subSystems{:}]');

iAstro_iPS_BD_T_vs_TP = struct('rxn_specific_T', {tbl_rxn_a}, 'subsys_specific_T', {tbl_subsys_a}, 'rxn_specific_TP', {tbl_rxn_b}, 'subsys_specific_TP', {tbl_subsys_b});

clear model_A model_B rxnList tf loc p idx subSystem tbl_rxn_a tbl_subsys_a tbl_rxn_b tbl_subsys_b

%% iAstro_iPS_BD_R (T vs TP)
model_A = iAstro_iPS_BD_R_T; model_B = iAstro_iPS_BD_R_TP;

rxnList = setdiff(model_A.rxns, model_B.rxns);
[tf,loc] = ismember(model_A.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_A.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_a = table(rxnList, subSystem);
tbl_subsys_a = setdiff([model_A.subSystems{:}]',[model_B.subSystems{:}]');

rxnList = setdiff(model_B.rxns, model_A.rxns);
[tf,loc] = ismember(model_B.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_B.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_b = table(rxnList, subSystem);
tbl_subsys_b = setdiff([model_B.subSystems{:}]',[model_A.subSystems{:}]');

iAstro_iPS_BD_R_T_vs_TP = struct('rxn_specific_T', {tbl_rxn_a}, 'subsys_specific_T', {tbl_subsys_a}, 'rxn_specific_TP', {tbl_rxn_b}, 'subsys_specific_TP', {tbl_subsys_b});

clear model_A model_B rxnList tf loc p idx subSystem tbl_rxn_a tbl_subsys_a tbl_rxn_b tbl_subsys_b

%% iAstro_iPS_BD_NR (T vs TP)
model_A = iAstro_iPS_BD_NR_T; model_B = iAstro_iPS_BD_NR_TP;

rxnList = setdiff(model_A.rxns, model_B.rxns);
[tf,loc] = ismember(model_A.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_A.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_a = table(rxnList, subSystem);
tbl_subsys_a = setdiff([model_A.subSystems{:}]',[model_B.subSystems{:}]');

rxnList = setdiff(model_B.rxns, model_A.rxns);
[tf,loc] = ismember(model_B.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model_B.subSystems(idx);
subSystem = [subSystem{:}]';
tbl_rxn_b = table(rxnList, subSystem);
tbl_subsys_b = setdiff([model_B.subSystems{:}]',[model_A.subSystems{:}]');

iAstro_iPS_BD_NR_T_vs_TP = struct('rxn_specific_T', {tbl_rxn_a}, 'subsys_specific_T', {tbl_subsys_a}, 'rxn_specific_TP', {tbl_rxn_b}, 'subsys_specific_TP', {tbl_subsys_b});

clear model_A model_B rxnList tf loc p idx subSystem tbl_rxn_a tbl_subsys_a tbl_rxn_b tbl_subsys_b

%%
writetable(table(iAstro_Primary_T.rxns, [iAstro_Primary_T.subSystems{:}]'),'Rxns_iAstro_Primary_T.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_Ctrl_T.rxns, [iAstro_iPS_Ctrl_T.subSystems{:}]'),'Rxns_iAstro_iPS_Ctrl_T.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_BD_T.rxns, [iAstro_iPS_BD_T.subSystems{:}]'),'Rxns_iAstro_iPS_BD_T.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_BD_R_T.rxns, [iAstro_iPS_BD_R_T.subSystems{:}]'),'Rxns_iAstro_iPS_BD_R_T.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_BD_NR_T.rxns, [iAstro_iPS_BD_NR_T.subSystems{:}]'),'Rxns_iAstro_iPS_BD_NR_T.csv','WriteVariableNames', true, 'Delimiter','\t');

writetable(table(iAstro_Primary_TP.rxns, [iAstro_Primary_TP.subSystems{:}]'),'Rxns_iAstro_Primary_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_Ctrl_TP.rxns, [iAstro_iPS_Ctrl_TP.subSystems{:}]'),'Rxns_iAstro_iPS_Ctrl_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_BD_TP.rxns, [iAstro_iPS_BD_TP.subSystems{:}]'),'Rxns_iAstro_iPS_BD_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_BD_R_TP.rxns, [iAstro_iPS_BD_R_TP.subSystems{:}]'),'Rxns_iAstro_iPS_BD_R_TP.csv','WriteVariableNames', true, 'Delimiter','\t');
writetable(table(iAstro_iPS_BD_NR_TP.rxns, [iAstro_iPS_BD_NR_TP.subSystems{:}]'),'Rxns_iAstro_iPS_BD_NR_TP.csv','WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except iAstro_Primary_T_vs_TP iAstro_iPS_Ctrl_T_vs_TP iAstro_iPS_BD_T_vs_TP iAstro_iPS_BD_R_T_vs_TP iAstro_iPS_BD_NR_T_vs_TP
save('iAstro_T_vs_TP_Comparison.mat');

