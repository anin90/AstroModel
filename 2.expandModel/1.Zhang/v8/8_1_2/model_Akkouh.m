function[rxnList_down, rxnList_up, model_Akkouh, FSr_significant_L, FSr_significant_H, subSystem_L, subSystem_H] = expose(model, filename);
tic;

%% LOAD WT MODEL
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% DATA QUALITY CONTROL

% IMPORT DATA
Akkouh = importdata(filename);

% FILTER_1: EXCLUDE DUPLICATE GENES in DATA
A = Akkouh.Gene_Symbol;
% A = str2double(A);
[~,X,Z] = unique(A,'stable');
Y = histc(Z,1:numel(X))<2;
set = A(X(Y))';

[intis_1,ind] = ismember(A',(set));
filter_1.Entrez_ID = Akkouh.Entrez_ID(intis_1==1);
filter_1.Gene_Symbol = Akkouh.Gene_Symbol(intis_1==1);
filter_1.Status = Akkouh.Status(intis_1==1);
filter_1.HPA_Level = Akkouh.HPA_Level(intis_1==1);

% FILTER_2: EXCLUDE GENES WITHOUT ENTREZ_ID
[intis_inter] = ~(isnan(filter_1.Entrez_ID));
[filter_2.Entrez_ID] = filter_1.Entrez_ID(find(intis_inter));
[filter_2.Gene_Symbol] = filter_1.Gene_Symbol(find(intis_inter));
[filter_2.Status] = filter_1.Status(find(intis_inter));
[filter_2.HPA_Level] = filter_1.HPA_Level(find(intis_inter));

% FILTER_3: EXCLUDE GENES WITHOUT HPA EVIDENCE
HPA_Level_NA = strfind(filter_2.HPA_Level,'NA');
HPA_Level_YES = find(cellfun(@isempty,HPA_Level_NA));
gene = (filter_2.Entrez_ID(HPA_Level_YES)); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
[filter_3.Entrez_ID] = filter_2.Entrez_ID(HPA_Level_YES);
[filter_3.Gene_Symbol] = filter_2.Gene_Symbol(HPA_Level_YES);
[filter_3.Status] = filter_2.Status(HPA_Level_YES);
[filter_3.HPA_Level] = filter_2.HPA_Level(HPA_Level_YES);

% SPLIT GENES (UP/DOWN)
up = strfind(filter_3.Status,'up');
up_yes = find(~cellfun(@isempty,up));
genes_up = filter_3.Entrez_ID(up_yes); genes_up = num2cell(genes_up); genes_up = cellfun(@num2str,genes_up,'uni',0);
down = strfind(filter_3.Status,'down');
down_yes = find(~cellfun(@isempty,down));
genes_down = filter_3.Entrez_ID(down_yes); genes_down = num2cell(genes_down); genes_down = cellfun(@num2str,genes_down,'uni',0);
genes_up = strcat(genes_up,'.1');
genes_down = strcat(genes_down,'.1');

% FIND REACTIONS FROM GENES
list_up = findRxnsFromGenes(model, genes_up);
list_up = struct2cell(list_up);
Aina_up = vertcat(list_up{:,1});
rxnList_up = Aina_up(:,1);
rxnList_up = unique(rxnList_up);

list_down = findRxnsFromGenes(model, genes_down);
list_down = struct2cell(list_down);
Aina_down = vertcat(list_down{:,1});
rxnList_down = Aina_down(:,1);
rxnList_down = unique(rxnList_down);

% REMOVE INTERSECTING REACTIONS
c = intersect(rxnList_up,rxnList_down);
rxnList_up(ismember(rxnList_up,c)) = [];
rxnList_down(ismember(rxnList_down,c)) = [];

% REMOVE BLOCKED-UP RXNS (DOESN'T HELP UPREGULATION PROBLEM)
% b = intersect(rxnList_up,fluxInconsistentRxns_ASM_BBB);
% rxnList_up(ismember(rxnList_up,b)) = [];

%% CODE TESTING
% SWAP TEST (T13_SWAP)
% lol1 = rxnList_up; lol2 = rxnList_down;
% clear rxnList_up rxnList_down
% rxnList_up = lol2; rxnList_down = lol1;
% clear lol1 lol2

% RANDOM SAMPLE TEST (T14)
% clear rxnList_up rxnList_down;
% rxnList_up = unique(datasample(model.rxns,1000));
% rxnList_down = unique(datasample(model.rxns,1000));
% d = intersect(rxnList_up,rxnList_down);
% rxnList_up(ismember(rxnList_up,d)) = [];
% rxnList_down(ismember(rxnList_down,d)) = [];

%% SPLIT REV RXNS (ONLY THOSE TO BE UPREGULATED) TO IRREV
% BiDirRxns = model.rxns(model.lb<0 & model.ub>0);
% BiDir_up = unique(intersect(rxnList_up,BiDirRxns));
% [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model, 'sRxns', BiDir_up);
% rxns_f_b = modelIrrev.rxns(matchRev~=0);

% UP_FWD_RXNS
% [id_f] = strfind(rxns_f_b,'_f');
% [ind_f] = find(not(cellfun('isempty',id_f)));
% rxnList_up_f = rxns_f_b(ind_f);

% UP_BACK_RXNS
% [id_b] = strfind(rxns_f_b,'_b');
% [ind_b] = find(not(cellfun('isempty',id_b)));
% rxnList_up_b = rxns_f_b(ind_b);

% UPDATE INPUT RXN LIST
% rxnList_down = rxnList_down; %DOWN
% UniDirRxns = model.rxns(model.lb==0 & model.ub>0);
% UniDirSubset = unique(intersect(rxnList_up,UniDirRxns));
% rxnList_up_f = vertcat(UniDirSubset,rxnList_up_f); %UP_FWD
% rxnList_up_b = rxnList_up_b; %UP_BACK

%% UPDATE INPUT MODEL
% model = modelIrrev;
% FVA of WT MODEL
[minFlux_WT, maxFlux_WT] = fluxVariability(model);

%% MODEL DOWNREGULATION [UB = MAX/2]
model_temp_all_MAX_DOWN = cell(rxnList_down);
model_temp_opt_all_MAX_DOWN = cell(rxnList_down);
model_ub_all_MAX_DOWN = cell(rxnList_down);

for k = 1:length(rxnList_down');
    model_temp_all_MAX_DOWN{k} = changeObjective(model, rxnList_down(k), 1);
    model_temp_opt_all_MAX_DOWN{k} = optimizeCbModel(model_temp_all_MAX_DOWN{k},'max');
    model_ub_all_MAX_DOWN{k} = model_temp_opt_all_MAX_DOWN{k}.f/2;
end

model_all_MAX_DOWN = model; [intis,k] = ismember(model_all_MAX_DOWN.rxns,(rxnList_down));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
model_all_MAX_DOWN.ub(idx) = cell2mat(model_ub_all_MAX_DOWN);

%% MODEL DOWNREGULATION [UB = MIN/2]
model_temp_all_MIN_DOWN = cell(rxnList_down);
model_temp_opt_all_MIN_DOWN = cell(rxnList_down);
model_lb_all_MIN_DOWN = cell(rxnList_down);

for k = 1:length(rxnList_down');
    model_temp_all_MIN_DOWN{k} = changeObjective(model, rxnList_down(k), 1);
    model_temp_opt_all_MIN_DOWN{k} = optimizeCbModel(model_temp_all_MIN_DOWN{k},'min');
    model_lb_all_MIN_DOWN{k} = model_temp_opt_all_MIN_DOWN{k}.f/2;
end

model_all_MIN_DOWN = model_all_MAX_DOWN; [intis,k] = ismember(model_all_MIN_DOWN.rxns,(rxnList_down));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
model_all_MIN_DOWN.lb(idx) = cell2mat(model_lb_all_MIN_DOWN);

%% MODEL UPREGULATION [LB = MAX*2]
% model_temp_all_MAX_UP = cell(rxnList_up_f);
% model_temp_opt_all_MAX_UP = cell(rxnList_up_f);
% model_lb_all_MAX_UP = cell(rxnList_up_f);
% 
% for k = 1:length(rxnList_up_f');
%     model_temp_all_MAX_UP{k} = changeObjective(model, rxnList_up_f(k), 1);
%     model_temp_opt_all_MAX_UP{k} = optimizeCbModel(model_temp_all_MAX_UP{k},'max');
%     model_lb_all_MAX_UP{k} = model_temp_opt_all_MAX_UP{k}.f*2;
% end
% 
% model_all_MAX_UP = model_all_MIN_DOWN; [intis,k] = ismember(model_all_MAX_UP.rxns,(rxnList_up_f));
% [~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
% model_all_MAX_UP.lb(idx) = cell2mat(model_lb_all_MAX_UP);

%% MODEL UPREGULATION [LB = MAX*2]
% model_temp_all_MIN_UP = cell(rxnList_up_b);
% model_temp_opt_all_MIN_UP = cell(rxnList_up_b);
% model_lb_all_MIN_UP = cell(rxnList_up_b);
% 
% for k = 1:length(rxnList_up_b');
%     model_temp_all_MIN_UP{k} = changeObjective(model, rxnList_up_b(k), 1);
%     model_temp_opt_all_MIN_UP{k} = optimizeCbModel(model_temp_all_MIN_UP{k},'max');
%     model_lb_all_MIN_UP{k} = model_temp_opt_all_MIN_UP{k}.f*2;
% end
% 
% model_all_MIN_UP = model_all_MAX_UP; [intis,k] = ismember(model_all_MIN_UP.rxns,(rxnList_up_b));
% [~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
% model_all_MIN_UP.lb(idx) = cell2mat(model_lb_all_MIN_UP);

%% MODEL UPREGULATION [LB = MIN*2]
model_temp_all_MIN_UP = cell(rxnList_up);
model_temp_opt_all_MIN_UP = cell(rxnList_up);
model_lb_all_MIN_UP = cell(rxnList_up);

for k = 1:length(rxnList_up');
    model_temp_all_MIN_UP{k} = changeObjective(model, rxnList_up(k), 1);
    model_temp_opt_all_MIN_UP{k} = optimizeCbModel(model_temp_all_MIN_UP{k},'min');
    model_lb_all_MIN_UP{k} = model_temp_opt_all_MIN_UP{k}.f*2;
end

model_all_MIN_UP = model_all_MIN_DOWN; [intis,k] = ismember(model_all_MIN_UP.rxns,(rxnList_up));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
model_all_MIN_UP.lb(idx) = cell2mat(model_lb_all_MIN_UP);

%% UPDATE CHANGES & FVA
model_Akkouh = model_all_MIN_UP;
[minFlux_Akkouh, maxFlux_Akkouh] = fluxVariability(model_Akkouh);

%% CALCULATE FLUX SPAN RATIO (FSr)
FS_WT = maxFlux_WT-minFlux_WT;
FS_Akkouh = bsxfun(@minus, maxFlux_Akkouh, minFlux_Akkouh);
FSr_Akkouh = bsxfun(@rdivide, abs(FS_WT), abs(FS_Akkouh));

%% REACTIONS WITH DIFFERENTIAL FSr
FSr_significant_L = model.rxns(FSr_Akkouh>1.5);
FSr_significant_H = model.rxns(FSr_Akkouh<0.8);

% REMOVE INPUT RXNS
rxnList_down_up = vertcat(rxnList_down,rxnList_up);
FSr_significant_L_H = vertcat(FSr_significant_L,FSr_significant_H);
c = intersect(rxnList_down_up,FSr_significant_L_H);
FSr_significant_L(ismember(FSr_significant_L,c)) = [];
FSr_significant_H(ismember(FSr_significant_L,c)) = [];

%% SUBSYSTEMS DISRUPTED
% FSr_significant_L
[intis, ~] = ismember(model.rxns, FSr_significant_L);
subSystem_L = model.subSystems(intis==1);

% FSr_significant_H
[intis, ~] = ismember(model.rxns, FSr_significant_H);
subSystem_H = model.subSystems(intis==1);

%%
toc;
end