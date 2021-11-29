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
model = iAstro_iPS_BD_TP; dat_rxns = dat.rxnList;
[mets, rxnFormula]= deal(repmat({''},size(dat_rxns))');
for i = 1:length(dat_rxns);
    [mets{i}] = findMetsFromRxns(model, dat_rxns{i})';
    [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
end
mets = cellfun(@(x) strjoin(x, ', '), mets', 'UniformOutput', false);
rxnFormula = rxnFormula';
clear model dat_rxns i

%%
clearvars -except bd_significant_lumped bd_significant_r bd_significant_nr

