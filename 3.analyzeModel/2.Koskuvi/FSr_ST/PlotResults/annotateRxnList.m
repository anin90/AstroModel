function [Tbl] = annotateRxnList(model, rxnList, FSrTbl)
%%
tStart = tic;

% Load MetabolicUnits from Recon3DModel_MetabolicUnits (Thiele et al. 2020)
FileName   = 'MetabolicUnits.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MetabolicUnits/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% Step-0: Invoke CobraTbx and Gurobi
initCobraToolbox(false); changeCobraSolver('gurobi','all'); 

%% Step-1: Annotate

    % subSystem and GPR from Recon3D (Brunk et al. 2018)
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    % FSr_a / FSr_b
    dat = FSrTbl.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];
    FluxSpanTbl = tbl_temp_fst(:,2:end);

    % High / Low
    dat = FSrTbl.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    Flux = repmat({''},size(nameidx_L));
    Flux(nameidx_L~=0) = {'L'};
    Flux(nameidx_L==0) = {'H'};

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
    
%% Step-2: FINAL TABLE
Tbl = [rxnList, subSystem, GPR, FluxSpanTbl, Flux, MetabolicUnits, Localization, RxnFormula];
Tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'Fluxspan_a', 'Fluxspan_b', 'FluxSpanRatio', 'Flux', 'MetabolicUnits', 'Localization', 'RxnFormula'};

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    
end