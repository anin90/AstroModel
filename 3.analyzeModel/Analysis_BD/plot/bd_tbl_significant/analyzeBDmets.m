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
model = iAstro_iPS_BD_TP; dat_rxns = bd_significant_lumped.rxnList;
[mets_bd_lumped, rxnFormula_bd_lumped]= deal(repmat({''},size(dat_rxns))');
for i = 1:length(dat_rxns);
    [mets_bd_lumped{i}] = findMetsFromRxns(model, dat_rxns{i})';
    [rxnFormula_bd_lumped{i}] = printRxnFormula(model, dat_rxns{i})';
end
mets_bd_lumped = cellfun(@(x) strjoin(x, ', '), mets_bd_lumped', 'UniformOutput', false);
rxnFormula_bd_lumped = rxnFormula_bd_lumped';

[mets, uniqueMetNames, uniqueMetIDs] = printUniqueMets(model, dat_rxns);

%% bd_significant_r
model = iAstro_iPS_BD_R_TP; dat_rxns = bd_significant_r.rxnList;
[mets_bd_r, rxnFormula_bd_r]= deal(repmat({''},size(dat_rxns))');
for i = 1:length(dat_rxns);
    [mets_bd_r{i}] = findMetsFromRxns(model, dat_rxns{i})';
    [rxnFormula_bd_r{i}] = printRxnFormula(model, dat_rxns{i})';
end
mets_bd_r = cellfun(@(x) strjoin(x, ', '), mets_bd_r', 'UniformOutput', false);
rxnFormula_bd_r = rxnFormula_bd_r';
clear model dat_rxns i

%% bd_significant_nr
model = iAstro_iPS_BD_NR_TP; dat_rxns = bd_significant_nr.rxnList;
[mets_bd_nr, rxnFormula_bd_nr]= deal(repmat({''},size(dat_rxns))');
for i = 1:length(dat_rxns);
    [mets_bd_nr{i}] = findMetsFromRxns(model, dat_rxns{i})';
    [rxnFormula_bd_nr{i}] = printRxnFormula(model, dat_rxns{i})';
end
mets_bd_nr = cellfun(@(x) strjoin(x, ', '), mets_bd_nr', 'UniformOutput', false);
rxnFormula_bd_nr = rxnFormula_bd_nr';
clear model dat_rxns i

%% write to table
bd_significant_lumped = [bd_significant_lumped rxnFormula_bd_lumped];
bd_significant_lumped.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization', 'rxnFormula'};

bd_significant_r = [bd_significant_r rxnFormula_bd_r];
bd_significant_r.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization', 'rxnFormula'};

bd_significant_nr = [bd_significant_nr rxnFormula_bd_nr];
bd_significant_nr.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization', 'rxnFormula'};

%%
clearvars -except bd_significant_lumped bd_significant_r bd_significant_nr uniqueMetIDs
%%
save('bd_significant_mets.mat');
%%
toc;
