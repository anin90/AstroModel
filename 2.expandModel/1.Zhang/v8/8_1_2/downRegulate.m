function[rxnList, models_single_down, model_double_down, model_all_down, FSr_significant] = downRegulate(model,filename);
tic;

%% Load WT model
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% perform FVA for WT model:
[minFlux_WT, maxFlux_WT] = fluxVariability(model);

%% Load rxns for downregulation
delimiterIn = '\t';
rxnList = importdata(filename,delimiterIn);

%% model_single downregulation:
model_temp = cell(rxnList);
model_temp_opt_max = cell(rxnList);
model_temp_opt_min = cell(rxnList);

for i = 1:length(rxnList');
    model_temp{i} = changeObjective(model, rxnList(i), 1);
    % constrain ub
    model_temp_opt_max{i} = optimizeCbModel(model_temp{i},'max');
    model_single_down_max{i} = changeRxnBounds(model,rxnList{i},model_temp_opt_max{i}.f/2,'u');
    % constrain lb
    model_temp_opt_min{i} = optimizeCbModel(model_temp{i},'min');
    model_single_down_min{i} = changeRxnBounds(model_single_down_max{i},rxnList{i},model_temp_opt_min{i}.f/2,'l');
    % save model and FVA
    models_single_down = model_single_down_min';
    [minFlux_single_down(:,i), maxFlux_single_down(:,i)] = fluxVariability(models_single_down{i});
end

%% model_double downregulation:
rxnList2 = {'MI1PP','INSTt4'};
model_temp_double = cell(rxnList2);
model_temp_opt_double_max = cell(rxnList2);
model_ub_double_max = cell(rxnList2);
model_temp_opt_double_min = cell(rxnList2);
model_lb_double_min = cell(rxnList2);

for j = 1:length(rxnList2');
    model_temp_double{j} = changeObjective(model, rxnList2(j), 1);
    % find max
    model_temp_opt_double_max{j} = optimizeCbModel(model_temp_double{j},'max');
    model_ub_double_max{j} = model_temp_opt_double_max{j}.f/2;
    % find min
    model_temp_opt_double_min{j} = optimizeCbModel(model_temp_double{j},'min');
    model_lb_double_min{j} = model_temp_opt_double_min{j}.f/2;
end

model_double_down = model; [intis,j] = ismember(model_double_down.rxns,(rxnList2));
[~,p] = sort(j(intis)); idx = find(intis); idx = idx(p);
model_double_down.ub(idx) = cell2mat(model_ub_double_max); % constrain ub
model_double_down.lb(idx) = cell2mat(model_lb_double_min); % constrain lb
[minFlux_double_down, maxFlux_double_down] = fluxVariability(model_double_down);

%% model_all downregulation:
model_temp_all = cell(rxnList);
model_temp_opt_all_max = cell(rxnList);
model_ub_all_max = cell(rxnList);
model_temp_opt_all_min = cell(rxnList);
model_ub_all_min = cell(rxnList);

for k = 1:length(rxnList');
    model_temp_all{k} = changeObjective(model, rxnList(k), 1);
    % find max
    model_temp_opt_all_max{k} = optimizeCbModel(model_temp_all{k},'max');
    model_ub_all_max{k} = model_temp_opt_all_max{k}.f/2;
    % find min
    model_temp_opt_all_min{k} = optimizeCbModel(model_temp_all{k},'min');
    model_ub_all_min{k} = model_temp_opt_all_min{k}.f/2;
end

model_all_down = model; [intis,k] = ismember(model_all_down.rxns,(rxnList));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
model_all_down.ub(idx) = cell2mat(model_ub_all_max);
model_all_down.lb(idx) = cell2mat(model_ub_all_min);
[minFlux_all_down, maxFlux_all_down] = fluxVariability(model_all_down);

%% FSr_single_down/double_down/all_down
FS_WT = maxFlux_WT-minFlux_WT;
FS_single_down = bsxfun(@minus, maxFlux_single_down, minFlux_single_down);
FSr_single_down = bsxfun(@rdivide, abs(FS_WT), abs(FS_single_down));

FS_double_down = bsxfun(@minus, maxFlux_double_down, minFlux_double_down);
FSr_double_down = bsxfun(@rdivide, abs(FS_WT), abs(FS_double_down));

FS_all_down = bsxfun(@minus, maxFlux_all_down, minFlux_all_down);
FSr_all_down = bsxfun(@rdivide, abs(FS_WT), abs(FS_all_down));

%% FSr_significant
% SINGLE_DOWN
for a = 1:length(models_single_down);
    FSr_significant(a).a = models_single_down{a}.rxns(FSr_single_down(:,a)>1.5);
    % FSr_significant(a).b = model_custom_KO.rxns(FSr_custom_KO>2 & FSr_custom_KO~=inf);
end;
FSr_significant = struct(FSr_significant); clear a;

% DOUBLE_DOWN
for a = 1:length(model_double_down);
%     FSr_significant(a).b = model_double_down{a}.rxns(FSr_double_down(:,a)>1.5);
    FSr_significant(a).b = model_double_down.rxns(FSr_double_down(:,a)>1.5);
end;
FSr_significant = struct(FSr_significant); clear a;

% ALL_DOWN
for a = 1:length(model_all_down);
%     FSr_significant(a).c = model_all_down{a}.rxns(FSr_all_down(:,a)>1.5);
    FSr_significant(a).c = model_all_down.rxns(FSr_all_down(:,a)>1.5);
end;
FSr_significant = struct(FSr_significant); clear a;
%%
toc;
end