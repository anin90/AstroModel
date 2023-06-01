function [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList)

%%
tStart = tic;

% Load Recon3DModel_301_MU
FileName   = 'Recon3DModel_301_MU.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MetabolicUnits/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

% Load MetabolicUnits from Recon3DModel_MetabolicUnits (Thiele et al. 2020)
FileName   = 'MetabolicUnits.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MetabolicUnits/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% Invoke CobraTbx and 'gurobi' (or) 'ibm_cplex'
% initCobraToolbox(false); changeCobraSolver('ibm_cplex','all'); 

%% Load model:
model = Recon3DModel_MetabolicUnits;
    
%% Annotate
    
    % MetabolicUnits from HarveyHarvetta (Thiele et al. 2020)
    List1 = rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};

    %% THIS BLOCK DOES NOT WORK FOR NON-UNIQUE rxnList
    % DONT USE THE OUTPUTS..
    
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
    
    %%
    
    % Reaction Formula
    [RxnFormula]= deal(repmat({''},size(rxnList))');
    for i = 1:length(rxnList);
        [RxnFormula{i}] = printRxnFormula(model, rxnList{i})';
    end
    RxnFormula = RxnFormula';

    % Reaction Name/ECNumbers
    for i = 1:length(rxnList);
        [tmp{i},~] = find(ismember(model.rxns,rxnList{i}));
        [RxnName{i}] = model.rxnNames(tmp{i});
        [RxnECNumbers{i}] = model.rxnECNumbers(tmp{i});
    end
    RxnName = RxnName';
    RxnECNumbers = RxnECNumbers';

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    
end
