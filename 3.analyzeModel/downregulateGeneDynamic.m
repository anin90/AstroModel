function[rxnList, activity, model_down] = downregulateGeneDynamic(model)
tic;

%% Load WT Model
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'Biomass_Ecoli_core_w_GAM', 0); %model

%% Defile ActivityLevel (null, 25%, 50%, 75%, wt)
[activity] = [0.00001 25/100 50/100 75/100 100/100];

%% Map Genes to Reactions in Model
[rxnList] = {'EX_glc(e)', 'FBP'}; 

%% Model Downregulation
model_down = model; 

% Preallocations
[model_temp, model_opt_max, model_newMax, model_opt_min, model_newMin] = deal(cell(rxnList));
model_down_newBounds = cell(size(activity)); 

for act = 1:length(activity);
    
    for r = 1:length(rxnList');
        
        % change objective to rxnList
        model_temp{r,act} = changeObjective(model, rxnList{r}, 1);
        
        % find max
        model_opt_max{r,act} = optimizeCbModel(model_temp{r,act}, 'max');
        model_newMax{r,act} = model_opt_max{r,act}.f*act;
    
        % find min
        model_opt_min{r,act} = optimizeCbModel(model_temp{r,act}, 'min');
        model_newMin{r,act} = model_opt_min{r,act}.f*act;
        
    end
    
    % update bounds on rxnList and create a model for each activity level (n=5)
    [intis,k] = ismember(model_down.rxns,(rxnList));
    [~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
    model_down_newBounds{act} = model;
    model_down_newBounds{act}.ub(idx) = cell2mat(model_newMax(:,act));
    model_down_newBounds{act}.lb(idx) = cell2mat(model_newMin(:,act));

end

%% Final changes
model_down = model_down_newBounds;
activity = activity';

end