tic;
%%
model = MBA_model_7_3_2;
changeCobraSolver('gurobi','all');

model = removeRxns(model, model.rxns(findExcRxns(model)==1));
%% add DMEM missing exchanges into the model
%% 
%model=addExchangeRxn(model,{'thm[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'cortsn[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'tststerone[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'fol[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'bilirub[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'prgstrn[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%%
model=addExchangeRxn(model,{'pydxn[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'PYDXNtr','reactionFormula','pydxn[e] <=> pydxn[c] ', 'subSystem','Transport, extracellular', 'geneRule','');

model=addExchangeRxn(model,{'prostge2[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'PROSTGE2t','reactionFormula','hco3[c] + prostge2[e] <=> hco3[e] + prostge2[c] ', 'subSystem','Transport, extracellular', 'geneRule','28231.1 or 6578.1 or 28232.1 or 6579.1');

model=addExchangeRxn(model,{'val_L[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'VALt4','reactionFormula','na1[e] + val_L[e] -> na1[c] + val_L[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1 or (340024.1 and 57393.1) or (340024.1 and 59272.1)');

model=addExchangeRxn(model,{'lys_L[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'LYSt4','reactionFormula','lys_L[e] + na1[e] -> lys_L[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1 or 6584.1');

model=addExchangeRxn(model,{'arg_L[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'ARGt4','reactionFormula','arg_L[e] + na1[e] -> arg_L[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1');

model=addExchangeRxn(model,{'chsterol[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'r1050','reactionFormula','chsterol[e] <=> chsterol[c] ', 'subSystem','Transport, extracellular', 'geneRule','');

model=addExchangeRxn(model,{'inost[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'INSTt4','reactionFormula','inost[e] + na1[e] <=> inost[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','6526.1');

model=addExchangeRxn(model,{'Nacasp[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'NACASPt','reactionFormula','3.0 na1[c] + Nacasp[c] <=> 3.0 na1[e] + Nacasp[e] ', 'subSystem','Transport, extracellular', 'geneRule','');

%model=addExchangeRxn(model,{'retinol_cis_11[e]'},-1000,1000); %NA_Ex;NA_t
%model = addReaction(model,'','reactionFormula','', 'subSystem','Transport, extracellular', 'geneRule','');

model=addExchangeRxn(model,{'chol[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'CHOLt4','reactionFormula','chol[e] + na1[e] <=> chol[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','60482.1');
%%
model=addExchangeRxn(model,{'na1[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'cl[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'hco3[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'avite1[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'ser_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'prostge1[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'gln_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'thr_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'leu_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'asn_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'glu_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'his_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'tyr_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'ribflv[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'pnto_R[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'so4[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'ca2[e]'},-1000,1000); %NA_Ex
%%
model=addExchangeRxn(model,{'glc_D[e]'},-1000,1000);
model=addExchangeRxn(model,{'gly[e]'},-1000,1000);
model=addExchangeRxn(model,{'creat[e]'},-1000,1000);
model=addExchangeRxn(model,{'ala_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'pro_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'cys_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'ile_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'asp_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'met_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'phe_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'pi[e]'},-1000,1000);
model=addExchangeRxn(model,{'trp_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'k[e]'},-1000,1000);
model=addExchangeRxn(model,{'urea[e]'},-1000,1000);

%% add synaptic exchanges (from addSynapseRxns_7_3_2.m)
model = addExchangeRxn(model, {'h[s]'}, -1000, 1000);
model = addExchangeRxn(model, {'glu_L[s]'}, -1000, 0);  %only uptake
model = addExchangeRxn(model, {'k[s]'}, -1000, 1000);
model = addExchangeRxn(model, {'gln_L[s]'}, 0, 1000);   %only release
model = addExchangeRxn(model, {'lac_L[s]'}, 0, 1000);   %only release
model = addExchangeRxn(model, {'ala_L[s]'}, 0, 1000);   %only release
model = addExchangeRxn(model, {'cit[s]'}, 0, 1000);     %only release  
model = addExchangeRxn(model, {'h2o[s]'}, -1000, 0);
model = addExchangeRxn(model, {'nh4[s]'}, 0, 1000);     %only release
model = addExchangeRxn(model, {'co2[e]'}, 0, 1000);     %only release   
model = addExchangeRxn(model, {'4mop[e]'}, -1000, 1000);
%%

%%
%fluxConsistency = verifyModel(model,'fluxConsistency',true); 
%fluxInconsistentRxns = model.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1)
%%
DMEM_model = model;
optimizeCbModel(DMEM_model,'max')
%%
clearvars -except AstroModel_Lewis_2010 MBA_model_7_3_1 MBA_model_7_3_2 TestSolutionName_Brain TestSolution_Brain Imbalanced_NonExcRxns DMEM_model
%%
toc;
