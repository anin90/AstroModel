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
% Load BD_Significant Rxns
bd_significant_lumped = importfile('bd_lumped_rxns_fdr.csv');
bd_significant_r = importfile('bd_r_rxns_fdr.csv');
bd_significant_nr = importfile('bd_nr_rxns_fdr.csv');

%% bd_significant_lumped
% rxns
model = iAstro_iPS_BD_TP; dat_rxns = bd_significant_lumped.rxnList;
[rxnFormula_bd_lumped]= deal(repmat({''},size(dat_rxns))');
for i = 1:length(dat_rxns);
    [rxnFormula_bd_lumped{i}] = printRxnFormula(model, dat_rxns{i})';
end
rxnFormula_bd_lumped = rxnFormula_bd_lumped';

% mets
[mets_bd_lumped] = printUniqueMets(model, dat_rxns);

%% bd_significant_r
model = iAstro_iPS_BD_R_TP; dat_rxns = bd_significant_r.rxnList;
[rxnFormula_bd_r]= deal(repmat({''},size(dat_rxns))');
for i = 1:length(dat_rxns);
    [rxnFormula_bd_r{i}] = printRxnFormula(model, dat_rxns{i})';
end
rxnFormula_bd_r = rxnFormula_bd_r';

% mets
[mets_bd_r] = printUniqueMets(model, dat_rxns);

%% bd_significant_nr
model = iAstro_iPS_BD_NR_TP; dat_rxns = bd_significant_nr.rxnList;
[rxnFormula_bd_nr]= deal(repmat({''},size(dat_rxns))');
for i = 1:length(dat_rxns);
    [rxnFormula_bd_nr{i}] = printRxnFormula(model, dat_rxns{i})';
end
rxnFormula_bd_nr = rxnFormula_bd_nr';

% mets
[mets_bd_nr] = printUniqueMets(model, dat_rxns);

%% write to table
% rxns
bd_significant_lumped = [bd_significant_lumped rxnFormula_bd_lumped];
bd_significant_lumped.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization', 'rxnFormula'};

bd_significant_r = [bd_significant_r rxnFormula_bd_r];
bd_significant_r.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization', 'rxnFormula'};

bd_significant_nr = [bd_significant_nr rxnFormula_bd_nr];
bd_significant_nr.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization', 'rxnFormula'};

% mets
writetable(mets_bd_lumped, 'mets_bd_lumped.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(mets_bd_r, 'mets_bd_r.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(mets_bd_nr, 'mets_bd_nr.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except bd_significant_lumped bd_significant_r bd_significant_nr mets_bd_lumped mets_bd_r mets_bd_nr
%%
save('bd_significant_mets.mat');
%%
toc;
