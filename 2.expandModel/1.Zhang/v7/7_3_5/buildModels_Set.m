function[rxnList, models_single_KO, model_custom_KO, FSr_single_KO, FSr_custom_KO, FSr_significant, FSr_mat_significant, FSr_significant_subSystems_Pooled, FSr_significant_KO_Specific_subSystems, FSr_significant_KO_Specific_Exchanges, FSr_significant_KO_Specific_Exchanges_Index] = ko_model(model,filename);
tic;

%% Load WT model
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% perform FVA for WT model:
[minFlux_WT, maxFlux_WT] = fluxVariability(model);

%% Load rxns for single_KO
delimiterIn = '\t';
rxnList = importdata(filename,delimiterIn);

%% build single_KO models and perform FVA
models_single_KO=cell(rxnList);
parfor ind = 1:length(rxnList');
    model_single_KO = changeRxnBounds(model,rxnList(ind),0,'b');
    models_single_KO{ind} = model_single_KO;
    [minFlux_single_KO(:,ind), maxFlux_single_KO(:,ind)] = fluxVariability(model_single_KO);    
    % save(['model_' num2str(ind) '.mat'],'model');
    % printFluxBounds(allModels{1},{'INSTt4'});
end

%% build custom_KO models and perform FVA
rxnList_all = rxnList;
model_custom_KO = changeRxnBounds(model,rxnList_all,0,'b');
[minFlux_custom_KO, maxFlux_custom_KO] = fluxVariability(model_custom_KO);

%% FSr_single/custom_KO
FS_WT = maxFlux_WT-minFlux_WT;
FS_single_KO = bsxfun(@minus, maxFlux_single_KO, minFlux_single_KO);
FSr_single_KO = bsxfun(@rdivide, abs(FS_WT), abs(FS_single_KO));
FS_custom_KO = bsxfun(@minus, maxFlux_custom_KO, minFlux_custom_KO);
FSr_custom_KO = bsxfun(@rdivide, abs(FS_WT), abs(FS_custom_KO));

%% FSr_significant
% SINGLE_KO
for a = 1:length(models_single_KO);
    FSr_significant(a).a = models_single_KO{a}.rxns(FSr_single_KO(:,a)>2);
    % FSr_significant(a).b = model_custom_KO.rxns(FSr_custom_KO>2 & FSr_custom_KO~=inf);
end;
FSr_significant = struct(FSr_significant); clear a;

% CUSTOM_KO
FSr_significant(1).b = model_custom_KO.rxns(FSr_custom_KO>2);

%% cat all FSr_significant rxns
%Pooled
FSr_cat_unique = unique(cat(1,FSr_significant(:).a,FSr_significant(:).b));
[~, FSr_cat_unique_indices] = intersect(model.rxns,FSr_cat_unique);
FSr_mat_all = cat(2,FSr_single_KO,FSr_custom_KO);
FSr_mat_significant= FSr_mat_all(FSr_cat_unique_indices,:);
FSr_significant_subSystems_Pooled = model.subSystems(FSr_cat_unique_indices);

%KO_specific
FSr_significant_KO_Specific_subSystems = struct2cell(FSr_significant); %% was cell(rxnList);
ans_temp = struct2cell(FSr_significant(:,1:end));
for x = 1:length(ans_temp);
    [~,ind_temp{x}] = intersect(model.rxns,ans_temp{x});
    FSr_significant_KO_Specific_subSystems{x} = model.subSystems(ind_temp{x});
end

%% identify the EX/DM/sink rxns in each KO model
% define exchanges
% modelexchanges1 = strmatch('Ex_',model.rxns);
% modelexchanges4 = strmatch('EX_',model.rxns);
% modelexchanges = unique([modelexchanges1;modelexchanges4]);
% ExRxns = unique(model.rxns(modelexchanges));
% 
% FSr_significant_KO_Specific_Exchanges = cell(rxnList);
% FSr_significant_KO_Specific_Exchanges_Index = cell(rxnList);
% temp2 = struct2cell(FSr_significant(:,1:end)); 
% for x = 1:length(temp2);
%     [~,FSr_significant_KO_Specific_Exchanges_Index{x}] = intersect(ExRxns,temp2{x});
%     FSr_significant_KO_Specific_Exchanges{x} = ExRxns(FSr_significant_KO_Specific_Exchanges_Index{x});
% end
% clear ans_temp ans x
% clear model modelexchanges1 modelexchanges4 modelexchanges2 modelexchanges3 selExc modelexchanges


%% plot Test_stats
dat = sum(FSr_mat_significant(1:end,:)>2);
bar(dat);
set(gca,'xticklabel',{'MI1PP','INSTt4','PGS','PDE1','BPNT','All'});
ax = gca;
ax.XTickLabelRotation = 45; ylim([0,1650]);
xlabel('KO-Model'); ylabel('#Rxns with FSr>2');
text(1:length(dat),dat,num2str(dat'),'vert','bottom','horiz','center');
hold on; breakyaxis([16 1600]); clear ans ax dat;

%%
toc;
end