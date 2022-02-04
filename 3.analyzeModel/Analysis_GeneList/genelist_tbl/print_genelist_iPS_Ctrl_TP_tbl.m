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

%% iAstro_iPS_Ctrl_TP 
%% CKMT2 (1160)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_CKMT2.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_CKMT2.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_CKMT2.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    CKMT2_tbl = [gene_tbl D];
    CKMT2_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% INPP5A (3632)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_INPP5A.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_INPP5A.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_INPP5A.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    INPP5A_tbl = [gene_tbl D];
    INPP5A_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% PLA2G6 (8398)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PLA2G6.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_PLA2G6.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_PLA2G6.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    PLA2G6_tbl = [gene_tbl D];
    PLA2G6_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% GBA (2629)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_GBA.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_GBA.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_GBA.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    GBA_tbl = [gene_tbl D];
    GBA_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):
    
%% PLCB4 (5332)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PLCB4.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_PLCB4.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_PLCB4.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    PLCB4_tbl = [gene_tbl D];
    PLCB4_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):    
    
%% NMNAT1 (64802)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_NMNAT1.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_NMNAT1.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_NMNAT1.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    NMNAT1_tbl = [gene_tbl D];
    NMNAT1_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% ASPA (443)
	% No mapped rxns in 'iAstro_iPS_Ctrl_TP'
        
%% SLC5A8 (160728)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SLC5A8.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_SLC5A8.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_SLC5A8.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    SLC5A8_tbl = [gene_tbl D];
	SLC5A8_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):    

%% PLA2G12A (81579)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PLA2G12A.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_PLA2G12A.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_PLA2G12A.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    PLA2G12A_tbl = [gene_tbl D];
    PLA2G12A_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% SLC29A2 (3177)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SLC29A2.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_SLC29A2.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_SLC29A2.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    SLC29A2_tbl = [gene_tbl D];
    SLC29A2_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% COASY (80347)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_COASY.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_COASY.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_COASY.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    COASY_tbl = [gene_tbl D];
    COASY_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% SYNJ2 (8871)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SYNJ2.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_SYNJ2.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_SYNJ2.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    SYNJ2_tbl = [gene_tbl D];
    SYNJ2_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% PYGL (5836)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_PYGL.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_PYGL.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_PYGL.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    PYGL_tbl = [gene_tbl D];
    PYGL_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):

%% NEU1 (4758)
	% No mapped rxns in 'iAstro_iPS_Ctrl_TP'
	
%% UGT2A1 (10941)
	% No mapped rxns in 'iAstro_iPS_Ctrl_TP' 			
    
%% SLC22A9 (114571)
    rxnList = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/FSR_iPS_Ctrl_TP_SLC22A9.csv'); 
    model = iAstro_iPS_Ctrl_TP;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model, rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_SLC22A9.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_SLC22A9.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    gene_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];

    % annotate rxns by organelle information
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
    SLC22A9_tbl = [gene_tbl D];
    SLC22A9_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'fluxspan_a', 'fluxspan_b', 'FluxSpanRatio', 'direction', 'localization'};
    % annotate rxns by mapping (direct or consequence):    

%% write to csv
writetable(CKMT2_tbl, 'CKMT2_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(INPP5A_tbl, 'INPP5A_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(PLA2G6_tbl, 'PLA2G6_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(GBA_tbl, 'GBA_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(PLCB4_tbl, 'PLCB4_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(NMNAT1_tbl, 'NMNAT1_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
% writetable(ASPA_tbl, 'ASPA_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(SLC5A8_tbl, 'SLC5A8_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(PLA2G12A_tbl, 'PLA2G12A_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(SLC29A2_tbl, 'SLC29A2_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(COASY_tbl, 'COASY_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(SYNJ2_tbl, 'SYNJ2_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(PYGL_tbl, 'PYGL_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
% writetable(NEU1_tbl, 'NEU1_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
% writetable(UGT2A1_tbl, 'UGT2A1_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(SLC22A9_tbl, 'SLC22A9_iPS_Ctrl_TP_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Retain Primary and iPS_Ctrl_TP models
clearvars -except iAstro_Primary_TP iAstro_iPS_Ctrl_TP


