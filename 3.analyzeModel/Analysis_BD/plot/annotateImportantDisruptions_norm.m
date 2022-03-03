%% Load Data
tic;
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName
% Load FluxSpanData
FileName   = 'iAstro_FluxDiff_BD.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName
% Load MetabolicUnits from Recon3DModel_MetabolicUnits
FileName   = 'MetabolicUnits.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/MetabolicUnits/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% bd_30_norm

    rxnList = importdata('bd_tbl_norm/bd_30_norm.csv'); 
    model = iAstro_iPS_Ctrl_TP_norm;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MU = repmat({''},size(nameidx_L));
    MU(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MU(nameidx_L==0) = {'NA'};

    bd_30_tbl = [tbl_temp tbl_temp_fst(:,2:end) C MU];
    clear rxnList model tf loc p idx subSystem GPR tbl_temp
    clear dat tbl_temp_fst A B ii List1 List2 List3 nameidx_L C MU

%% bd_38_norm

    rxnList = importdata('bd_tbl_norm/bd_38_norm.csv');
    model = iAstro_iPS_Ctrl_TP_norm;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MU = repmat({''},size(nameidx_L));
    MU(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MU(nameidx_L==0) = {'NA'};

    bd_38_tbl = [tbl_temp tbl_temp_fst(:,2:end) C MU];
    clear rxnList model tf loc p idx subSystem GPR tbl_temp
    clear dat tbl_temp_fst A B ii List1 List2 List3 nameidx_L C MU

%% bd_36_norm

    rxnList = importdata('bd_tbl_norm/bd_36_norm.csv');
    model = iAstro_iPS_Ctrl_TP_norm;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MU = repmat({''},size(nameidx_L));
    MU(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MU(nameidx_L==0) = {'NA'};

    bd_36_tbl = [tbl_temp tbl_temp_fst(:,2:end) C MU];
    clear rxnList model tf loc p idx subSystem GPR tbl_temp
    clear dat tbl_temp_fst A B ii List1 List2 List3 nameidx_L C MU

%% bd_95_norm

    rxnList = importdata('bd_tbl_norm/bd_95_norm.csv');
    model = iAstro_iPS_Ctrl_TP_norm;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_TP_norm.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MU = repmat({''},size(nameidx_L));
    MU(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MU(nameidx_L==0) = {'NA'};

    bd_95_tbl = [tbl_temp tbl_temp_fst(:,2:end) C MU];
    clear rxnList model tf loc p idx subSystem GPR tbl_temp
    clear dat tbl_temp_fst A B ii List1 List2 List3 nameidx_L C MU

%% bd_r_155_norm

    rxnList = importdata('bd_tbl_norm/bd_r_155_norm.csv');
    model = iAstro_iPS_Ctrl_TP_norm;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_R_TP_norm.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_R_TP_norm.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MU = repmat({''},size(nameidx_L));
    MU(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MU(nameidx_L==0) = {'NA'};

    bd_r_155_tbl = [tbl_temp tbl_temp_fst(:,2:end) C MU];
    clear rxnList model tf loc p idx subSystem GPR tbl_temp
    clear dat tbl_temp_fst A B ii List1 List2 List3 nameidx_L C MU

%% bd_nr_148_norm

    rxnList = importdata('bd_tbl_norm/bd_nr_148_norm.csv');
    model = iAstro_iPS_Ctrl_TP_norm;
    [tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
    subSystem = model.subSystems(idx);
    subSystem = [subSystem{:}]';
    GPR = findGPRFromRxns(model,rxnList);
    tbl_temp = table(rxnList, subSystem, GPR);

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_NR_TP_norm.fluxSpanTable;
    tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
    A = tbl_temp; B = tbl_temp_fst;
    [~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
    tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
    tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

    dat = FluxDiff_iPS_Ctrl_TP_norm_vs_iPS_BD_NR_TP_norm.fluxGoneLow;
    List1 = A.rxnList;
    List2 = dat.FSr_significant_L;
    nameidx_L = getnameidx(List2, List1)';
    C = repmat({''},size(nameidx_L));
    C(nameidx_L~=0) = {'L'};
    C(nameidx_L==0) = {'H'};

    List1 = A.rxnList;
    List2 = MetabolicUnits_ACS.HarvettaRxns;
    List3 = MetabolicUnits_ACS.MetabolicUnits;
    nameidx_L = getnameidx(List2, List1)';
    MU = repmat({''},size(nameidx_L));
    MU(nameidx_L~=0) = List3(nonzeros(nameidx_L));
    MU(nameidx_L==0) = {'NA'};

    bd_nr_148_tbl = [tbl_temp tbl_temp_fst(:,2:end) C MU];
    clear rxnList model tf loc p idx subSystem GPR tbl_temp
    clear dat tbl_temp_fst A B ii List1 List2 List3 nameidx_L C MU

%% merge bd_lumped; bd_r; bd_nr 
    bd_199_tbl = vertcat(bd_30_tbl, bd_38_tbl, bd_36_tbl, bd_95_tbl);
    bd_r_185_tbl = vertcat(bd_r_155_tbl, bd_30_tbl);
    bd_nr_184_tbl = vertcat(bd_nr_148_tbl, bd_36_tbl);

%% annotate rxns by organelle information
    model = iAstro_iPS_Ctrl_TP_norm;
    [transportRxnBool] = transportReactionBool(model);
    trspRxns = model.rxns(transportRxnBool==1);
    trspSubSys = model.subSystems(transportRxnBool==1);

% bd_199_tbl
    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(bd_199_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), bd_199_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(bd_199_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    bd_199_tbl = [bd_199_tbl D];
    clear transportRxnBool patterns i rxns rxns_subset nameidx idx D

% bd_r_185_tbl
    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(bd_r_185_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), bd_r_185_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(bd_r_185_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    bd_r_185_tbl = [bd_r_185_tbl D];
    clear transportRxnBool patterns i rxns rxns_subset nameidx idx D

% bd_nr_184_tbl
    patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
    D = repmat({''},size(bd_nr_184_tbl.rxnList));
    for i = 1:length(patterns)
        [rxns{i}] = findRxnFromCompartment(model,patterns{i});
        [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), bd_nr_184_tbl.rxnList);
        if ~isempty(rxns_subset{i})
            nameidx{i} = getnameidx(bd_nr_184_tbl.rxnList, rxns_subset{i})';
            D(nameidx{1,i}) = {patterns{i}};
        else % do nothing
        end
    end
    idx = cellfun(@isempty,D);
    D(idx) = {'intercompartmental'};
    bd_nr_184_tbl = [bd_nr_184_tbl D];
    clear transportRxnBool patterns i rxns rxns_subset nameidx idx D

%% write to csv
bd_199_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'Fluxspan_a', 'Fluxspan_b', 'FluxSpanRatio', 'Flux', 'MetabolicUnits', 'Localization'};
bd_r_185_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'Fluxspan_a', 'Fluxspan_b', 'FluxSpanRatio', 'Flux', 'MetabolicUnits', 'Localization'};
bd_nr_184_tbl.Properties.VariableNames = {'rxnList', 'subSystem', 'GPR', 'Fluxspan_a', 'Fluxspan_b', 'FluxSpanRatio', 'Flux', 'MetabolicUnits', 'Localization'};

writetable(bd_199_tbl, 'bd_tbl_norm/bd_199_tbl_norm.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(bd_r_185_tbl, 'bd_tbl_norm/bd_r_185_tbl_norm.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(bd_nr_184_tbl, 'bd_tbl_norm/bd_nr_184_tbl_norm.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except bd_199_tbl bd_r_185_tbl bd_nr_184_tbl

%%
save('bd_tbl_norm/bd_tbl_norm.mat');