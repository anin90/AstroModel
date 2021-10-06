tic;
%%
model = MBA_model_7_3_2;
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, Elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
%%
Rxns_imBalancedCharge = model.rxns(imBalancedCharge~=0); %Reactions with ImBalanced "charge"
Rxns_imBalancedMass = model.rxns(find(~cellfun(@isempty,imBalancedMass))); %Reactions with ImBalanced "mass"
Rxns_imBalanced = model.rxns(imBalancedRxnBool~=0); % all Imbalanced reactions
Mets_ImBalancedRxns = model.mets(balancedMetBool==0); %metabolites involved in ImBalanced reactions
Mets_MissingFormulae = model.mets(missingFormulaeBool~=0); %metabolites with missing formulae
checkMassChargeBalance_All = struct('Rxns_imBalancedCharge', {Rxns_imBalancedCharge},'Rxns_imBalancedMass', {Rxns_imBalancedMass}, 'Rxns_imBalanced',{Rxns_imBalanced} ,'Mets_ImBalancedRxns', {Mets_ImBalancedRxns}, 'Mets_MissingFormulae', {Mets_MissingFormulae});

%%
NonExcRxns_imBalancedCharge = intersect(Rxns_imBalancedCharge,model.rxns(findExcRxns(model)==0));
NonExcRxns_imBalancedMass = intersect(Rxns_imBalancedMass,model.rxns(findExcRxns(model)==0));
NonExcRxns_imBalanced = intersect(Rxns_imBalanced,model.rxns(findExcRxns(model)==0));
checkMassChargeBalance_NonExcRxns = struct('NonExcRxns_imBalancedCharge', {NonExcRxns_imBalancedCharge},'NonExcRxns_imBalancedMass', {NonExcRxns_imBalancedMass}, 'NonExcRxns_imBalanced',{NonExcRxns_imBalanced});

%% Below results same as "NonExcRxns_imBalanced", just to verify!
Imbalanced_NonExcRxns = unique(cat(1,NonExcRxns_imBalancedCharge,NonExcRxns_imBalancedMass,NonExcRxns_imBalanced));
%%
clearvars -except AstroModel_Lewis_2010 MBA_model_7_3_1 MBA_model_7_3_2 TestSolutionName_Brain TestSolution_Brain Imbalanced_NonExcRxns
%%
toc;

