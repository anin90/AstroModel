function [mta_tbl_filt, del_rxnID, mta_score] = generateMTAscores_b(model, exprs)

%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MTA_ImNotaGit/matlab/iMAT/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MTA_ImNotaGit/matlab/MTA/')

% Load MetabolicUnits from Recon3DModel_MetabolicUnits (Thiele et al. 2020)
FileName   = 'MetabolicUnits.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MetabolicUnits/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% Step-0: Invoke CobraTbx and Gurobi
initCobraToolbox(false); changeCobraSolver('gurobi','all'); 

%% Step-1: Import and Prepare Source Model
[model.rowlb, model.rowub] = deal(zeros(size(model.mets)));
[model.int_vars] = zeros(size(model.rxns));

genes = cellfun(@(x) strsplit(x, '.'), model.genes(:,1), 'UniformOutput', false);
for i = 1:length(genes)
    genes{i} = genes{i,1}{1,1};
end
[model.genes_unique] = unique(genes);
[model.genes] = sort(genes);
[model.genes_unique_map] = findgroups(genes);

%% Step-2: Import Expression data
source = exprs.source; target = exprs.target;
smean = exprs.smean; tmean = exprs.tmean;

%% Step-3: Discretize the DEGs & Generate Vref
[bin] = get_binary_expH(model, source, 1);
[~, Vref] = sampleiMAT(model, bin, 'gurobi', 'gurobi');

%% Step-4: Generate Discretized Rxn Vector
[discrete_rxns_vector, ~] = createDiscreteRxns_human(source, target, smean, tmean, model, 1);

%% Step-5: Run MTA
[mta_tbl] = MTA_ImNotaGit(model, Vref, discrete_rxns_vector, 1:length(model.rxns), 'gurobi');

%% Step-6: Filter (Top 20% predictions)
[mta_tbl_filt] = mta_tbl([mta_tbl.mta_score] > prctile([mta_tbl.mta_score]', 80));
[del_rxnID] = model.rxns([mta_tbl_filt.del_rxn]');
[mta_score] = [mta_tbl_filt.mta_score]';

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

end
