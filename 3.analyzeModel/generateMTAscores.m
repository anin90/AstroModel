function [mta_tbl_filt] = generateMTAscores(model, diff_exprs)

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

%% Step-1: Import Source model
[model] = changeObjective(model, 'biomass_maintenance', 0);
[model.rowlb, model.rowub] = deal(zeros(size(model.mets)));
[model.int_vars] = zeros(size(model.rxns));
% [bin] = get_binary_expH(model,GSE8671.healthy_tab(:,1),0);

%% Step-2: Discretize the DEGs
model = changeObjective(model, 'biomass_maintenance', 1);
Vref = optimizeCbModel(model,'max');
Vref = Vref.full;
rxnFBS = diffexprs2rxnFBS(model, diff_exprs, Vref);

%% Step-3: Generate Vref
[sample_points, Vref] = sampleiMAT(model, rxnFBS, 'gurobi', 'gurobi');

%% Step-4: Run MTA
[mta_tbl] = MTA_ImNotaGit(model, Vref, rxnFBS, 1:length(model.rxns), 'gurobi');

%% Step-5: Filter (Top 20% predictions)
[mta_tbl_filt] = mta_tbl([mta_tbl.mta_score] > prctile([mta_tbl.mta_score]', 80));
[del_rxnID] = model.rxns([mta_tbl_filt.del_rxn]');
[mta_score] = [mta_tbl_filt.mta_score]';

%% Step-6: Annotation
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
    MU = repmat({''},size(nameidx_L));
    MU(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MU(nameidx_L==0) = {'NA'};

    % Organelle from Recon3D (Brunk et al. 2018)
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);
    
    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};

%% Final Table
mta_tbl_filt = table(del_rxnID, mta_score, subSystem, GPR, MU, D);
mta_tbl_filt.Properties.VariableNames = {'del_rxnID_KO', 'mta_score', 'subSystem', 'GPR', 'MetabolicUnits', 'Localization'};

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

end
