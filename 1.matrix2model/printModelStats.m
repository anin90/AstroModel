function[modelStats] = printModelStats(model, expMat)
%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/')

%% Load Lewis et al. rxns
    astroModelLewis = xls2model('AstroModel_Lewis_et_al_2010.xlsx');
    
%% Print

    % # of METS, RXNS, GENES 
    modelMets = length(model.mets);
    modelRxns = length(model.rxns);
    modelGenes = length(model.genes);
    
    % I: fluxConsistency
    fluxConsistency = verifyModel(model,'fluxConsistency',true); 
    fluxInconsistentRxns = length(model.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1));
    fluxInconsistentRxnsPrct = (fluxInconsistentRxns/modelRxns)*100;
    
    % II: CoreRxns
    coreRxns =  printCoreRxns(expMat);
    overlapCoreRxns = length(intersect(model.rxns, coreRxns));
    overlapCoreRxnsPrct = (overlapCoreRxns/length(coreRxns))*100;
    
    % III: overlapLewis
    astroModelLewisRxns = length(astroModelLewis.rxns);
    overlapLewis = length(intersect(lower(model.rxns),lower(astroModelLewis.rxns)));
    overlapLewisPrct = (overlapLewis/astroModelLewisRxns)*100;
   
%% Final Table
modelStats = table(modelMets, modelRxns, modelGenes, fluxInconsistentRxns, fluxInconsistentRxnsPrct, length(coreRxns), overlapCoreRxns, overlapCoreRxnsPrct, astroModelLewisRxns, overlapLewis, overlapLewisPrct);
modelStats.Properties.VariableNames = {'modelMets', 'modelRxns', 'modelGenes', 'fluxInconsistentRxns', 'fluxInconsistentRxnsPrct', 'coreRxns', 'overlapCoreRxns', 'overlapCoreRxnsPrct', 'astroModelLewisRxns', 'overlapLewis', 'overlapLewisPrct'};

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

end