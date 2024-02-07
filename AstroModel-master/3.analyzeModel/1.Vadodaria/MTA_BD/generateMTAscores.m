function [mta_tbl] = generateMTAscores(model, dat_source, dat_target)

%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MTA_ImNotaGit/matlab/iMAT/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MTA_ImNotaGit/matlab/MTA/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/')

% Load MetabolicUnits from Recon3DModel_MetabolicUnits (Thiele et al. 2020)
FileName   = 'MetabolicUnits.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MetabolicUnits/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% Step-0: Invoke CobraTbx and Gurobi
initCobraToolbox(false); changeCobraSolver('gurobi','all'); 

%% Step-1: Import and Prepare Source Model
model = prepareModel(model);

%% Step-2: Import and Prepare Expression data
expression_source = prepareExpressionMatrix(dat_source);
expression_target = prepareExpressionMatrix(dat_target);
[exprs.source] = mapExpressionToModelGenes(model, expression_source); % SOURCE PHENOTYPE
[exprs.target] = mapExpressionToModelGenes(model, expression_target); % TARGET PHENOTYPE
source = exprs.source; target = exprs.target;
smean = mean(exprs.source,2); tmean = mean(exprs.target,2);

%% Step-3: Discretize the DEGs & Generate Vref
[bin] = get_binary_expH(model, source, 1);
[~, Vref] = sampleiMAT(model, bin, 'gurobi', 'gurobi');

%% Step-4: Generate Discretized Rxn Vector
[discrete_rxns_vector, ~] = createDiscreteRxns_human(source, target, smean, tmean, model, 1);

%% Step-5: Run MTA
[mta_tbl] = MTA_ImNotaGit(model, Vref, discrete_rxns_vector, 1:length(model.rxns), 'gurobi');

%% Step-6: Filter (Top 20% predictions)
% [mta_tbl_filt] = mta_tbl([mta_tbl.mta_score] > prctile([mta_tbl.mta_score]', 80));
[del_rxnID] = model.rxns([mta_tbl.del_rxn]');
[mta_score] = [mta_tbl.mta_score]';
[alt_score] = [mta_tbl.alt_score]';

%% Step-7: Annotation
rxnList = del_rxnID;

    % subSystem and GPR from Recon3D (Brunk et al. 2018)
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);

    % MetabolicUnits from HarveyHarvetta (Thiele et al. 2020)
    List1 = rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};

    % Organelle from Recon3D (Brunk et al. 2018)
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);
    
    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    Localization = repmat({''},size(rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(rxnList, rxns_subset{i})';
            Localization(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,Localization);
    Localization(idx) = {'intercompartmental'};
    
    % Reaction Formula
    [RxnFormula]= deal(repmat({''},size(rxnList))');
    for i = 1:length(rxnList);
        [RxnFormula{i}] = printRxnFormula(model, rxnList{i})';
    end
    RxnFormula = RxnFormula'; 

%% Final Table
mta_tbl = table(del_rxnID, mta_score, alt_score, subSystem, GPR, MetabolicUnits, Localization, RxnFormula);
mta_tbl.Properties.VariableNames = {'del_rxnID_KO', 'mta_score', 'alt_score', 'subSystem', 'GPR', 'MetabolicUnits', 'Localization', 'RxnFormula'};

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

end
