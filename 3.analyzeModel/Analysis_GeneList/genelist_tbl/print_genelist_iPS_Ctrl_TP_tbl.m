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
% Load FluxSpanData
FileName   = 'iAstro_FluxDiff_GeneList_TP.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName
% Load MetabolicUnits from Recon3DModel_MetabolicUnits
FileName   = 'MetabolicUnits.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations_MetabolicUnits/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% iAstro_iPS_Ctrl_TP 
%% CKMT2 (1160)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_CKMT2.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    
    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_CKMT2.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_CKMT2.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_CKMT2;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    CKMT2_tbl = [gene_tbl D rxnFormula];
    CKMT2_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(CKMT2_tbl, 'CKMT2_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% INPP5A (3632)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_INPP5A.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    
    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_INPP5A.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_INPP5A.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_INPP5A;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    INPP5A_tbl = [gene_tbl D rxnFormula];
    INPP5A_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(INPP5A_tbl, 'INPP5A_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% PLA2G6 (8398)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PLA2G6.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_PLA2G6.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_PLA2G6.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_PLA2G6;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    PLA2G6_tbl = [gene_tbl D rxnFormula];
    PLA2G6_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(PLA2G6_tbl, 'PLA2G6_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% GBA (2629)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_GBA.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_GBA.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_GBA.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_GBA;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    GBA_tbl = [gene_tbl D rxnFormula];
    GBA_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(GBA_tbl, 'GBA_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% PLCB4 (5332)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PLCB4.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_PLCB4.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_PLCB4.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_PLCB4;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    PLCB4_tbl = [gene_tbl D rxnFormula];
    PLCB4_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(PLCB4_tbl, 'PLCB4_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% NMNAT1 (64802)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_NMNAT1.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_NMNAT1.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_NMNAT1.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_NMNAT1;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    NMNAT1_tbl = [gene_tbl D rxnFormula];
    NMNAT1_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(NMNAT1_tbl, 'NMNAT1_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% ASPA (443)
	% No mapped rxns in 'iAstro_iPS_Ctrl_TP'
    
%% SLC5A8 (160728)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SLC5A8.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_SLC5A8.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_SLC5A8.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_SLC5A8;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    SLC5A8_tbl = [gene_tbl D rxnFormula];
    SLC5A8_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(SLC5A8_tbl, 'SLC5A8_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% PLA2G12A (81579)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PLA2G12A.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_PLA2G12A.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_PLA2G12A.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_PLA2G12A;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    PLA2G12A_tbl = [gene_tbl D rxnFormula];
    PLA2G12A_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(PLA2G12A_tbl, 'PLA2G12A_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% SLC29A2 (3177)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SLC29A2.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_SLC29A2.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_SLC29A2.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_SLC29A2;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    SLC29A2_tbl = [gene_tbl D rxnFormula];
    SLC29A2_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(SLC29A2_tbl, 'SLC29A2_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% COASY (80347)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_COASY.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_COASY.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_COASY.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_COASY;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    COASY_tbl = [gene_tbl D rxnFormula];
    COASY_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(COASY_tbl, 'COASY_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% SYNJ2 (8871)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SYNJ2.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_SYNJ2.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_SYNJ2.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_SYNJ2;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    SYNJ2_tbl = [gene_tbl D rxnFormula];
    SYNJ2_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(SYNJ2_tbl, 'SYNJ2_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% PYGL (5836)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PYGL.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_PYGL.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_PYGL.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_PYGL;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    PYGL_tbl = [gene_tbl D rxnFormula];
    PYGL_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(PYGL_tbl, 'PYGL_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% NEU1 (4758)
	% No mapped rxns in 'iAstro_iPS_Ctrl_TP'
	
%% UGT2A1 (10941)
	% No mapped rxns in 'iAstro_iPS_Ctrl_TP'	
    
%% SLC22A9 (114571)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SLC22A9.csv'); 
    model = iAstro_iPS_Ctrl_TP;

    % Annotate subSystem and GPR
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

	% Print FluxSpan
    dat = FluxDiff_iPS_Ctrl_TP_SLC22A9.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

	% Print Flux (High/Low, when compared to WT)
    dat = FluxDiff_iPS_Ctrl_TP_SLC22A9.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};
    
    % Annotate Metabolic Units based on HarveyHarvetta
    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MetabolicUnits = repmat({''},size(nameidx_L));
    MetabolicUnits(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MetabolicUnits(nameidx_L==0) = {'NA'};
    
    % PerturbationClass (DirectMapping or Cascade):
    dat = FluxDiff_iPS_Ctrl_TP_SLC22A9;
    List1 = A.rxnList;
    List2 = dat.GeneToRxnList;
    nameidx_L = getnameidx(List2, List1)';
    PerturbationClass = repmat({''},size(nameidx_L));
    PerturbationClass(nameidx_L~=0) = {'DirectMapping'};
    PerturbationClass(nameidx_L==0) = {'Cascade'};
    
    % Annotate Reaction formula:
    dat_rxns = A.rxnList;
    [rxnFormula]= deal(repmat({''},size(dat_rxns))');
    for i = 1:length(dat_rxns);
        [rxnFormula{i}] = printRxnFormula(model, dat_rxns{i})';
    end
    rxnFormula = rxnFormula';

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) PerturbationClass C MetabolicUnits];

    % Annotate Organelle Info
    model = iAstro_iPS_Ctrl_TP;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(gene_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), gene_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(gene_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    SLC22A9_tbl = [gene_tbl D rxnFormula];
    SLC22A9_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'PerturbationClass', 'Flux', 'MetabolicUnits', 'localization', 'rxnFormula'};
    writetable(SLC22A9_tbl, 'SLC22A9_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Retain Primary and iPS_Ctrl_TP models
clearvars -except iAstro_Primary_TP iAstro_iPS_Ctrl_TP


