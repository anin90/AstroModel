%delete all synapseRxns from model_7.3.1
[SynapseReactions] = findRxnFromCompartment(MBA_model_7_3_1,'[s]');
model = removeRxns(MBA_model_7_3_1, SynapseReactions);
%fluxConsistency = verifyModel(model,'fluxConsistency',true); model.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1)
findRxnFromCompartment(model,'[s]');

%add transport Rxns
model = addReaction(model,'GLUt6[s]','reactionFormula','h[s] + 3 na1[s] + k[c] + glu_L[s]  -> h[c] + 3 na1[c] + glu_L[c] + k[s] ', 'subSystem','Transport, Synapse', 'geneRule','6511.1 or 6505.1 or 6507.1 or 6506.1 or 6512.1');
model = addReaction(model,'NAt[s]','reactionFormula','na1[s] <=> na1[c] ', 'subSystem','Transport, Synapse', 'geneRule','6526.1 or 6523.1 or 6528.1');
model = addReaction(model,'r2525[s]','reactionFormula','gln_L[c] <=> gln_L[s] ', 'subSystem','Transport, Synapse', 'geneRule','8501.1');
model = addReaction(model,'L_LACt2r[s]','reactionFormula','h[c] + lac_L[c] <=> h[s] + lac_L[s] ', 'subSystem','Transport, Synapse', 'geneRule','9194.1 or 23539.1 or 6566.1 or 9123.1 or 682.1');
model = addReaction(model,'ALAPAT4te[s]','reactionFormula','ala_L[c] <=> ala_L[s] ', 'subSystem','Transport, Synapse', 'geneRule','120103.1');
model = addReaction(model,'CITt4_2[s]','reactionFormula','cit[c] + 2.0 na1[c] <=> cit[s] + 2.0 na1[s] ', 'subSystem','Transport, Synapse', 'geneRule','9058.1');
model = addReaction(model,'HMR_9802[s]','reactionFormula','gln_L[s] + h2o[s] -> glu_L[s] + nh4[s] ', 'subSystem','Glutamate metabolism', 'geneRule','27165.1');
model = addReaction(model,'CO2t','reactionFormula','co2[e] <=> co2[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'PGI','reactionFormula','g6p[c] <=> f6p[c] ', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','2821.1');
model = addReaction(model,'TPI','reactionFormula','dhap[c] <=> g3p[c] ', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','7167.1 or 286016.1');
model = addReaction(model,'PGM','reactionFormula','2pg[c] <=> 3pg[c] ', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','669.1 or 5223.1 or 5224.2 or 5224.1 or 669.2');
model = addReaction(model,'PYK','reactionFormula','adp[c] + h[c] + pep[c] -> atp[c] + pyr[c] ', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','5315.2 or 5313.2 or 5315.3 or 5313.1 or 5315.1');
model = addReaction(model,'ME2','reactionFormula','mal_L[c] + nadp[c] -> co2[c] + nadph[c] + pyr[c] ', 'subSystem','Pyruvate metabolism', 'geneRule','4199.1');
model = addReaction(model,'PCm','reactionFormula','atp[m] + hco3[m] + pyr[m] -> adp[m] + h[m] + oaa[m] + pi[m] ', 'subSystem','Pyruvate metabolism', 'geneRule','5091.1 or 5091.2');
model = addReaction(model,'r0941','reactionFormula','hco3[c] <=> hco3[m] ', 'subSystem','Transport, mitochondrial', 'geneRule','');
model = addReaction(model,'PYRt2m','reactionFormula','h[c] + pyr[c] <=> h[m] + pyr[m] ', 'subSystem','Transport, mitochondrial', 'geneRule','6566.1');
model = addReaction(model,'ICDHxm','reactionFormula','icit[m] + nad[m] -> akg[m] + co2[m] + nadh[m] ', 'subSystem','Citric acid cycle', 'geneRule','(3421.2 and 3420.3 and 3419.1) or (3420.1 and 3421.1 and 3419.1) or (3420.1 and 3421.1 and 3419.1) or (3420.3 and 3421.1 and 3419.1) or (3421.1 and 3419.1 and 3420.2) or (3421.2 and 3419.1 and 3420.2)');
model = addReaction(model,'ICDHy','reactionFormula','icit[c] + nadp[c] -> akg[c] + co2[c] + nadph[c] ', 'subSystem','Citric acid cycle', 'geneRule','3417.1');
model = addReaction(model,'SUCD1m','reactionFormula','fad[m] + succ[m] <=> fadh2[m] + fum[m] ', 'subSystem','Citric acid cycle', 'geneRule','6389.1 and 6392.1 and 6391.1 and 6390.1');
model = addReaction(model,'CITtam','reactionFormula','cit[c] + mal_L[m] <=> cit[m] + mal_L[c] ', 'subSystem','Transport, mitochondrial', 'geneRule','6576.1');
model = addReaction(model,'r0081','reactionFormula','akg[m] + ala_L[m] <=> glu_L[m] + pyr[m] ', 'subSystem','Citric acid cycle', 'geneRule','84706.1 or 2875.1');
model = addReaction(model,'ILETA','reactionFormula','akg[c] + ile_L[c] <=> 3mop[c] + glu_L[c] ', 'subSystem','Valine, leucine, and isoleucine metabolism', 'geneRule','586.1');
model = addReaction(model,'LEUTA','reactionFormula','akg[c] + leu_L[c] <=> 4mop[c] + glu_L[c] ', 'subSystem','Valine, leucine, and isoleucine metabolism', 'geneRule','586.1');
model = addReaction(model,'4MOPte','reactionFormula','4mop[e] <=> 4mop[c] ', 'subSystem','Transport, extracellular', 'geneRule','');

%add exchange Rxns
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

%remove Rxns from model
model = removeRxns(model, {'EX_melatn[e]','EX_hxan[e]','HYXNt','EX_gua[e]','sink_melatn[c]','r1437','CITt4_2'});

MBA_model = model;

%check fluxConsistency and identify DeadEnds
fluxConsistency = verifyModel(MBA_model,'fluxConsistency',true); MBA_model.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1);
DeadEndMets = MBA_model.mets(detectDeadEnds(MBA_model));

clear SynapseReactions;
clear fluxConsistency;
clear model;
