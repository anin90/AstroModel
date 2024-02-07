function[query_rxns,rxns, ind] = checkRxnsInModel(model,filename);
%% usage
% [~,rxns,~] = checkRxnsInModel(MBA_model_7_3_4,'rxnList_expandModel.csv');
% rxns will return a boolean vector (input rxns present/ absent in model)..

%% tic;

delimiterIn = '\t';
query_rxns = importdata(filename,delimiterIn);
[rxns, ind] = ismember(query_rxns,model.rxns);

toc;
end