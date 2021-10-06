function[PowRxnList, model_KO, FSr_KO, FSr_significant, FSr_mat_significant, FSr_significant_subSystems_Pooled, FSr_significant_KO_Specific_subSystems, FSr_significant_KO_Specific_Exchanges, FSr_significant_KO_Specific_Exchanges_Index] = KO_Model_Pow(model,filename);
tic;

%% Load WT model
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% perform FVA for WT model:
[minFlux_WT, maxFlux_WT] = fluxVariability(model);

%% Load rxns for single_KO
delimiterIn = '\t';
rxnList = importdata(filename,delimiterIn);
PowRxnList = PowerSet(rxnList)';

%% Generate KO_models
model_KO = cell(PowRxnList);
parfor i = 1:length(PowRxnList);
    model_temp = changeRxnBounds(model,PowRxnList{(i)}, 0, 'b'); % also works with {i}
    model_KO{i} = model_temp;
    [minFlux_KO(:,i), maxFlux_KO(:,i)] = fluxVariability(model_KO{i});
end

%% FSr_KO
FS_WT = maxFlux_WT-minFlux_WT;
FS_KO = bsxfun(@minus, maxFlux_KO, minFlux_KO);
FSr_KO = bsxfun(@rdivide, abs(FS_WT), abs(FS_KO));

%% FSr_significant
for a = 1:length(model_KO);
    FSr_significant(a).a = model_KO{a}.rxns(FSr_KO(:,a)>1.5);
    % FSr_significant(a).b = model_custom_KO.rxns(FSr_custom_KO>2 & FSr_custom_KO~=inf);
end;
FSr_significant = struct(FSr_significant); clear a;

%% cat all FSr_significant rxns, identify disrupted subSystems
%Pooled
FSr_cat_unique = unique(cat(1,FSr_significant(:).a));
[~, FSr_cat_unique_indices] = intersect(model.rxns,FSr_cat_unique);
FSr_mat_significant= FSr_KO(FSr_cat_unique_indices,:);
FSr_significant_subSystems_Pooled = model.subSystems(FSr_cat_unique_indices);

%KO_specific
FSr_significant_KO_Specific_subSystems = cell(PowRxnList);
temp1 = struct2cell(FSr_significant(:,1:end));
for x = 1:length(temp1);
    [~,ind_temp{x}] = intersect(model.rxns,temp1{x});
    FSr_significant_KO_Specific_subSystems{x} = model.subSystems(ind_temp{x});
end

%% identify the EX/DM/sink rxns in each KO model
% define exchanges
modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelexchanges = unique([modelexchanges1;modelexchanges4]);
ExRxns = unique(model.rxns(modelexchanges));

FSr_significant_KO_Specific_Exchanges = cell(PowRxnList);
FSr_significant_KO_Specific_Exchanges_Index = cell(PowRxnList);
temp2 = struct2cell(FSr_significant(:,1:end));
for x = 1:length(temp2);
    [~,FSr_significant_KO_Specific_Exchanges_Index{x}] = intersect(ExRxns,temp2{x});
    FSr_significant_KO_Specific_Exchanges{x} = ExRxns(FSr_significant_KO_Specific_Exchanges_Index{x});
end
% clear ans_temp ans x
% clear model modelexchanges1 modelexchanges4 modelexchanges2 modelexchanges3 selExc modelexchanges

%% plot Test_stats
dat = sum(FSr_mat_significant(1:end,:)>1.5);
bar(dat, 0.5);
set(gca,'xticklabel',{'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10'});
ax = gca;
% ax.XTickLabelRotation = 45;
ylim([0,1650]);
xlabel('KO-Model'); ylabel('#Rxns with FSr>1.5');
% text(1:length(dat),dat,num2str(dat'),'vert','bottom','horiz','center');
hold on;
breakyaxis([18 1600]); clear ans ax dat;

%%
toc;
end