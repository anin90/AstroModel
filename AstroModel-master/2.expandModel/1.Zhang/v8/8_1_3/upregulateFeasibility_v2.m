function[deltaModel, deltaModelRelaxed, idx, relaxed_UpRxns_LB, relaxed_UpRxns_UB, relaxed_AllRxns_LB, relaxed_AllRxns_UB] = upregulate_v2(model, rxnList_up);
%%
tic;
changeCobraSolver('gurobi','all'); clear ans;
%%
% model = iMAT_model_TP_EXP_ASM_BBB;
newModel = cell(rxnList_up);
newModelOpt = cell(rxnList_up);
newModelOpt_LB = cell(rxnList_up);

for k = 1:length(rxnList_up');
    newModel{k} = changeObjective(model, rxnList_up(k), 1);
    newModelOpt{k} = optimizeCbModel(newModel{k}, 'max');
    newModelOpt_LB{k} = newModelOpt{k}.f*2;
end

deltaModel = model;
[intis,k] = ismember(deltaModel.rxns,(rxnList_up));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
deltaModel.lb(idx) = cell2mat(newModelOpt_LB);

for k = 1:length(idx)
    if deltaModel.lb(idx(k)) > deltaModel.ub(idx(k));
        deltaModel.lb(idx(k)) = deltaModel.ub(idx(k));
    else
        % do nothing
    end
end

printFluxBounds(deltaModel,rxnList_up);

%% RUN RELAXED FBA

% add objective
deltaModel.biomassBool=strcmp(deltaModel.rxns,'biomass_maintenance');
deltaModel.c(deltaModel.biomassBool)=1;

% check if objective is feasible
FBAsolution = optimizeCbModel(deltaModel,'max');
if FBAsolution.stat == 1
    disp('Relaxed model is feasible');
    bioMassProductionRate=FBAsolution.x(deltaModel.biomassBool);
    fprintf('%g%s\n', bioMassProductionRate, ' is the biomass production rate')
else
    disp('Relaxed model is infeasible');
end

% add parameters
relaxOption.internalRelax = 2; % allow to relax bounds on all internal reactions
relaxOption.exchangeRelax = 0; % do not allow to relax bounds on exchange reactions
relaxOption.steadyStateRelax = 0; % do not allow to relax the steady state constraint S*v = b
feasTol = getCobraSolverParams('LP', 'feasTol');
relaxOption.epsilon = feasTol/100;%*100;

% run relaxedFBA
tic;
[solution] = relaxedFBA(deltaModel, relaxOption)
timeTaken=toc;
[v,r,p,q] = deal(solution.v,solution.r,solution.p,solution.q);

% incorporate relaxed constraints onto the model
if solution.stat == 1
    modelRelaxed=deltaModel;
    delta=0;%can be used for debugging, in case more relaxation is necessary
    modelRelaxed.lb = deltaModel.lb - p - delta;
    modelRelaxed.ub = deltaModel.ub + q + delta;
    modelRelaxed.b = deltaModel.b - r;
    FBAsolution = optimizeCbModel(modelRelaxed,'max', 0, true);
    if FBAsolution.stat == 1
        disp('Relaxed model is feasible');
    else
        disp('Relaxed model is infeasible');
        solutionRelaxed = relaxedFBA(modelRelaxed,relaxOption);
    end
end

%% IDENTIFY REACTIONS WHO'S BOUNDS ARE RELAXED
if modelRelaxed.lb(idx) == deltaModel.lb(idx)
    disp('lower bounds are same')
else
    
    disp('lower bounds are different')
    
    % UP_RXNS
    relaxed_UpIdx_LB = deltaModel.lb(idx) == modelRelaxed.lb(idx);
    relaxed_UpRxns_LB = modelRelaxed.rxns(idx(relaxed_UpIdx_LB~=1));
    relaxed_UpIdx_UB = deltaModel.ub(idx) == modelRelaxed.ub(idx);
    relaxed_UpRxns_UB = modelRelaxed.rxns(idx(relaxed_UpIdx_UB~=1));
    
    % ALL RXNS
    relaxed_AllIdx_LB = deltaModel.lb == modelRelaxed.lb;
    relaxed_AllRxns_LB = modelRelaxed.rxns(relaxed_AllIdx_LB~=1);
    relaxed_AllIdx_UB = deltaModel.ub(idx) == modelRelaxed.ub(idx);
    relaxed_AllRxns_UB = modelRelaxed.rxns(relaxed_AllIdx_UB~=1);
    
end

%% MAKE SURE NO OBJECTIVES ARE PLACED IN THE FINAL MODEL
model = changeObjective(model, 'biomass_maintenance', 0);
deltaModel = changeObjective(deltaModel, 'biomass_maintenance', 0);
deltaModelRelaxed = changeObjective(modelRelaxed, 'biomass_maintenance', 0);

%%

[min, max] = fluxVariability(deltaModelRelaxed);

% clear model newModel newModelOpt newModelOpt_LB deltaModel k

%%
toc; 
end