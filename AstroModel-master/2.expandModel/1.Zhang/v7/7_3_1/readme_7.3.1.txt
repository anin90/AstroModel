Build on model v7.3.1:
Last_update: 04/02/2020
**********************************************************************
>> MBA_model = addReaction(MBA_model,'EX_glc_D[e]','reactionFormula','glc_D[e] <=> ', 'subSystem','Exchange/demand reaction', 'geneRule','');
>> MBA_model = addReaction(MBA_model,'EX_o2[e]','reactionFormula','o2[e] <=> ', 'subSystem','Exchange/demand reaction', 'geneRule','');
>> MBA_model = addReaction(MBA_model,'GLCt4','reactionFormula','glc_D[e] + na1[e] <=> glc_D[c] + na1[c] ', 'subSystem','Exchange/demand reaction', 'geneRule','115584.1 or 6526.1 or 6524.1 or 125206.1 or 159963.1 or 200010.1');
>> MBA_model = addReaction(MBA_model,'DM_atp_c_','reactionFormula','atp[c] + h2o[c] -> adp[c] + h[c] + pi[c] ', 'subSystem','Exchange/demand reaction', 'geneRule','');
>> MBA_model = addReaction(MBA_model,'EX_gln_L[s]','reactionFormula','gln_L[s] <=> ', 'subSystem','Exchange/demand reaction', 'geneRule','');
>> MBA_model = addReaction(MBA_model,'EX_glu_L[s]','reactionFormula','glu_L[s] <=> ', 'subSystem','Exchange/demand reaction', 'geneRule','');
>> MBA_model = changeRxnBounds(MBA_model, 'EX_o2[e]', 0, 'u');
>> MBA_model = addReaction(MBA_model,'CITt4_2','reactionFormula','cit[e] + 2.0 na1[e] <=> cit[c] + 2.0 na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','9058.1');
>> MBA_model = addReaction(MBA_model,'CITt4_2[s]','reactionFormula','cit[c] + 2.0 na1[c] <=> cit[s] + 2.0 na1[s] ', 'subSystem','Transport, synapse', 'geneRule','9058.1');
**********************************************************************
