function[models_single_KO, model_custom_KO, FSr_single_KO, FSr_custom_KO, FSr_significant, FSr_mat_significant, FSr_significant_subSystems] = multiModels(model,filename1);
tic;

%% Load WT model
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% perform FVA for WT model:
[minFlux_WT, maxFlux_WT] = fluxVariability(model);

%% Load rxns for single_KO
filename1 = filename1;
delimiterIn = '\t';
rxnList = importdata(filename1,delimiterIn);

%% build single_KO models and perform FVA
models_single_KO=cell(rxnList);
parfor ind = 1:length(rxnList');
    model_single_KO = changeRxnBounds(model,rxnList(ind),0,'b');
    models_single_KO{ind} = model_single_KO;
    [minFlux_single_KO(:,ind), maxFlux_single_KO(:,ind)] = fluxVariability(model_single_KO);    
    % save(['model_' num2str(ind) '.mat'],'model');
    % printFluxBounds(allModels{1},{'INSTt4'});
end

%% build double_KO models and perform FVA
% rxnList2 = [{'MI1PP','INSTt4'};{'CDIPTr','MI1PS'};{'PIK4','RE2223M'}];
% rxnList2 = rxnList2';
% models_double_KO=cell(rxnList2(1,:)');
% parfor indj = 1:length(rxnList2);
%     model_double_KO = changeRxnBounds(model,rxnList2(:,indj),0,'b');
%     models_double_KO{indj} = model_double_KO;
%     [minFlux_double_KO(:,indj), maxFlux_double_KO(:,indj)] = fluxVariability(model_double_KO);   
% end

%% build custom_KO models and perform FVA
rxnList3 = rxnList;
model_custom_KO = changeRxnBounds(model,rxnList3,0,'b');
[minFlux_custom_KO, maxFlux_custom_KO] = fluxVariability(model_custom_KO);

%% FSr_single/custom_KO
FS_WT = maxFlux_WT-minFlux_WT;
FS_single_KO = bsxfun(@minus, maxFlux_single_KO, minFlux_single_KO);
FSr_single_KO = bsxfun(@rdivide, abs(FS_WT), abs(FS_single_KO));
FS_custom_KO = bsxfun(@minus, maxFlux_custom_KO, minFlux_custom_KO);
FSr_custom_KO = bsxfun(@rdivide, abs(FS_WT), abs(FS_custom_KO));

%% FSr_significant
% SINGLE_KO
parfor a = 1:length(models_single_KO);
    FSr_significant(a).a = models_single_KO{a}.rxns(FSr_single_KO(:,a)>2);
    % FSr_significant(a).b = model_custom_KO.rxns(FSr_custom_KO>2 & FSr_custom_KO~=inf);
end;
FSr_significant = struct(FSr_significant); clear a;

% CUSTOM_KO
FSr_significant(1).b = model_custom_KO.rxns(FSr_custom_KO>2);

% DOUBLE_KO
%parfor a = 1:length(models_double_KO);
    % FSr_significant(a).b = models_double_KO{a}.rxns(FSr_double_KO(:,a)>2 & FSr_double_KO(:,a)~=inf);
%end;
%FSr_significant = struct(FSr_significant); clear a;

%% cat all FSr_significant rxns
% FSr_cat_unique = unique(cat(1,FSr_significant(:).a,FSr_significant(:).b,FSr_significant(:).c));
FSr_cat_unique = unique(cat(1,FSr_significant(:).a,FSr_significant(:).b));
[~, FSr_cat_unique_indices] = intersect(model.rxns,FSr_cat_unique);
FSr_mat_all = cat(2,FSr_single_KO,FSr_custom_KO);
FSr_mat_significant= FSr_mat_all(FSr_cat_unique_indices,:);
FSr_significant_subSystems = model.subSystems(FSr_cat_unique_indices);

%% plot Test_stats
dat = sum(FSr_mat_significant(1:end,:)> 2);
bar(dat);
set(gca,'xticklabel',{'MI1PP','INSTt4','MI1PS','CDIPTr','PIK4','PI4P5K','PI45PLC','MI145PP','MI14P4P','DAGK_hs','CDS','MI145PK','MI1345PKc','MI13456PK','All'});
ax = gca;
ax.XTickLabelRotation = 45; clear ax; ylim([0,205]);
xlabel('KO-Model'); ylabel('#Rxns with FSr>2');
text(1:length(dat),dat,num2str(dat'),'vert','bottom','horiz','center'); clear dat;
hold on; breakyaxis([20 180]); clear ans;

%%
toc;
end