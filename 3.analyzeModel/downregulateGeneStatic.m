function[gene, rxnList, model_down] = downregulateGeneStatic(model, EntrezID)
tic;

%% Load WT model
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% Load genes for downregulation
gene = strcat(EntrezID,'.1');

%% map genes to reactions in model
[rxnList] = findRxnsActiveWithGenes(model, gene)';
[rxnList] = unique(rxnList);

%% model_all downregulation:
model_temp = cell(rxnList);
model_temp_opt_max = cell(rxnList);
model_ub_max = cell(rxnList);
model_temp_opt_min = cell(rxnList);
model_ub_min = cell(rxnList);

for k = 1:length(rxnList');
    model_temp{k} = changeObjective(model, rxnList(k), 1);
    
    % find max
    model_temp_opt_max{k} = optimizeCbModel(model_temp{k},'max');
    model_ub_max{k} = model_temp_opt_max{k}.f/2;
    
    % find min
    model_temp_opt_min{k} = optimizeCbModel(model_temp{k},'min');
    model_ub_min{k} = model_temp_opt_min{k}.f/2;
end

model_down = model; [intis,k] = ismember(model_down.rxns,(rxnList));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
model_down.ub(idx) = cell2mat(model_ub_max);
model_down.lb(idx) = cell2mat(model_ub_min);

toc;
end
