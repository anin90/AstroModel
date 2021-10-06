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
FileName   = 'iAstro_FluxDiff_BD.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% bd_17
rxnList = importdata('bd_17.csv'); 
model = iAstro_iPS_Ctrl_TP;
[tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model.subSystems(idx);
subSystem = [subSystem{:}]';
GPR = findGPRFromRxns(model,rxnList);
tbl_temp = table(rxnList, subSystem, GPR);

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxSpanTable;
tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
A = tbl_temp; B = tbl_temp_fst;
[~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxGoneLow;
List1 = A.rxnList;
List2 = dat.FSr_significant_L;
nameidx_L = getnameidx(List2, List1)';
C = repmat({''},size(nameidx_L));
C(nameidx_L~=0) = {'L'};
C(nameidx_L==0) = {'H'};

bd_17_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];
clear rxnList model tf loc p idx subSystem GPR tbl_temp
clear dat tbl_temp_fst A B ii List1 List2 nameidx_L C

%% bd_33
rxnList = importdata('bd_33.csv');
model = iAstro_iPS_Ctrl_TP;
[tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model.subSystems(idx);
subSystem = [subSystem{:}]';
GPR = findGPRFromRxns(model,rxnList);
tbl_temp = table(rxnList, subSystem, GPR);

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxSpanTable;
tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
A = tbl_temp; B = tbl_temp_fst;
[~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxGoneLow;
List1 = A.rxnList;
List2 = dat.FSr_significant_L;
nameidx_L = getnameidx(List2, List1)';
C = repmat({''},size(nameidx_L));
C(nameidx_L~=0) = {'L'};
C(nameidx_L==0) = {'H'};

bd_33_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];
clear rxnList model tf loc p idx subSystem GPR tbl_temp
clear dat tbl_temp_fst A B ii List1 List2 nameidx_L C

%% bd_59
rxnList = importdata('bd_59.csv');
model = iAstro_iPS_Ctrl_TP;
[tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model.subSystems(idx);
subSystem = [subSystem{:}]';
GPR = findGPRFromRxns(model,rxnList);
tbl_temp = table(rxnList, subSystem, GPR);

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxSpanTable;
tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
A = tbl_temp; B = tbl_temp_fst;
[~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxGoneLow;
List1 = A.rxnList;
List2 = dat.FSr_significant_L;
nameidx_L = getnameidx(List2, List1)';
C = repmat({''},size(nameidx_L));
C(nameidx_L~=0) = {'L'};
C(nameidx_L==0) = {'H'};

bd_59_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];
clear rxnList model tf loc p idx subSystem GPR tbl_temp
clear dat tbl_temp_fst A B ii List1 List2 nameidx_L C

%% bd_103
rxnList = importdata('bd_103.csv');
model = iAstro_iPS_Ctrl_TP;
[tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model.subSystems(idx);
subSystem = [subSystem{:}]';
GPR = findGPRFromRxns(model,rxnList);
tbl_temp = table(rxnList, subSystem, GPR);

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxSpanTable;
tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
A = tbl_temp; B = tbl_temp_fst;
[~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_TP.fluxGoneLow;
List1 = A.rxnList;
List2 = dat.FSr_significant_L;
nameidx_L = getnameidx(List2, List1)';
C = repmat({''},size(nameidx_L));
C(nameidx_L~=0) = {'L'};
C(nameidx_L==0) = {'H'};

bd_103_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];
clear rxnList model tf loc p idx subSystem GPR tbl_temp
clear dat tbl_temp_fst A B ii List1 List2 nameidx_L C

%% bd_r_92
rxnList = importdata('bd_r_92.csv');
model = iAstro_iPS_Ctrl_TP;
[tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model.subSystems(idx);
subSystem = [subSystem{:}]';
GPR = findGPRFromRxns(model,rxnList);
tbl_temp = table(rxnList, subSystem, GPR);

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_R_TP.fluxSpanTable;
tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
A = tbl_temp; B = tbl_temp_fst;
[~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_R_TP.fluxGoneLow;
List1 = A.rxnList;
List2 = dat.FSr_significant_L;
nameidx_L = getnameidx(List2, List1)';
C = repmat({''},size(nameidx_L));
C(nameidx_L~=0) = {'L'};
C(nameidx_L==0) = {'H'};

bd_r_92_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];
clear rxnList model tf loc p idx subSystem GPR tbl_temp
clear dat tbl_temp_fst A B ii List1 List2 nameidx_L C

%% bd_nr_670
rxnList = importdata('bd_nr_670.csv');
model = iAstro_iPS_Ctrl_TP;
[tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem = model.subSystems(idx);
subSystem = [subSystem{:}]';
GPR = findGPRFromRxns(model,rxnList);
tbl_temp = table(rxnList, subSystem, GPR);

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_NR_TP.fluxSpanTable;
tbl_temp_fst = dat(ismember(dat.Rxn, tbl_temp.rxnList)==1,:);
A = tbl_temp; B = tbl_temp_fst;
[~,ii] = ismember(A.rxnList,B.Rxn); tbl_temp_fst = B(ii,:);
tbl_temp_fst.minFlux_a=[]; tbl_temp_fst.maxFlux_a=[];
tbl_temp_fst.minFlux_b=[]; tbl_temp_fst.maxFlux_b=[];

dat = FluxDiff_iPS_Ctrl_TP_vs_iPS_BD_NR_TP.fluxGoneLow;
List1 = A.rxnList;
List2 = dat.FSr_significant_L;
nameidx_L = getnameidx(List2, List1)';
C = repmat({''},size(nameidx_L));
C(nameidx_L~=0) = {'L'};
C(nameidx_L==0) = {'H'};

bd_nr_670_tbl = [tbl_temp tbl_temp_fst(:,2:end) C];
clear rxnList model tf loc p idx subSystem GPR tbl_temp
clear dat tbl_temp_fst A B ii List1 List2 nameidx_L C

%% merge bd_17_33_59_103
bd_212_tbl = vertcat(bd_17_tbl, bd_33_tbl, bd_59_tbl, bd_103_tbl);

%% annotate rxns by organelle information
model = iAstro_iPS_Ctrl_TP;
[transportRxnBool] = transportReactionBool(model);
trspRxns = model.rxns(transportRxnBool==1);
trspSubSys = model.subSystems(transportRxnBool==1);

% bd_212_tbl
patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
D = repmat({''},size(bd_212_tbl.rxnList));
for i = 1:length(patterns)
    [rxns{i}] = findRxnFromCompartment(model,patterns{i});
    [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), bd_212_tbl.rxnList);
    if ~isempty(rxns_subset{i})
        nameidx{i} = getnameidx(bd_212_tbl.rxnList, rxns_subset{i})';
        D(nameidx{1,i}) = {patterns{i}};
    else % do nothing
    end
end
idx = cellfun(@isempty,D);
D(idx) = {'intercompartmental'};
bd_212_tbl = [bd_212_tbl D];
clear transportRxnBool patterns i rxns rxns_subset nameidx idx D

% bd_r_92_tbl
patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
D = repmat({''},size(bd_r_92_tbl.rxnList));
for i = 1:length(patterns)
    [rxns{i}] = findRxnFromCompartment(model,patterns{i});
    [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), bd_r_92_tbl.rxnList);
    if ~isempty(rxns_subset{i})
        nameidx{i} = getnameidx(bd_r_92_tbl.rxnList, rxns_subset{i})';
        D(nameidx{1,i}) = {patterns{i}};
    else % do nothing
    end
end
idx = cellfun(@isempty,D);
D(idx) = {'intercompartmental'};
bd_r_92_tbl = [bd_r_92_tbl D];
clear transportRxnBool patterns i rxns rxns_subset nameidx idx D

% bd_nr_670_tbl
patterns = {'[c]' '[l]' '[n]' '[m]' '[i]' '[x]' '[e]' '[g]' '[r]'};
D = repmat({''},size(bd_nr_670_tbl.rxnList));
for i = 1:length(patterns)
    [rxns{i}] = findRxnFromCompartment(model,patterns{i});
    [rxns_subset{i}] = intersect(setdiff(rxns{i}(:,1), trspRxns), bd_nr_670_tbl.rxnList);
    if ~isempty(rxns_subset{i})
        nameidx{i} = getnameidx(bd_nr_670_tbl.rxnList, rxns_subset{i})';
        D(nameidx{1,i}) = {patterns{i}};
    else % do nothing
    end
end
idx = cellfun(@isempty,D);
D(idx) = {'intercompartmental'};
bd_nr_670_tbl = [bd_nr_670_tbl D];
clear transportRxnBool patterns i rxns rxns_subset nameidx idx D

%% write to csv
writetable(bd_212_tbl, 'bd_212_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(bd_r_92_tbl, 'bd_r_92_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');
writetable(bd_nr_670_tbl, 'bd_nr_670_tbl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except bd_212_tbl bd_r_92_tbl bd_nr_670_tbl

%%
save('bd_tbl.mat');