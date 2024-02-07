function[model_EXP_Unconstrained, fluxInconsistentRxns_non_media, model_EXP_Constrained, fluxInconsistentRxns_media] = expandModel(model);
tic;
changeCobraSolver('gurobi','all'); clear ans;
%% change lb of all sink/EX_rxns to 10:

%number of EX_, DM_ & sink_ rxns;
modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelexchanges2 = strmatch('DM_',model.rxns);
modelexchanges3 = strmatch('sink_',model.rxns);
selExc = (find( full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))))';
modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
Total_ExRxns = unique(model.rxns(modelexchanges));

% input: EX, sink rxns
EX_Rxns = Total_ExRxns(strmatch('EX_',Total_ExRxns));
Sink_Rxns = Total_ExRxns(strmatch('sink_',Total_ExRxns));
model=changeRxnBounds(model,EX_Rxns,-10,'l');
model=changeRxnBounds(model,Sink_Rxns,-10,'l');

%% add biomass_maintenance
model = addReaction(model,'biomass_maintenance','reactionFormula','0.50563 ala_L[c] + 0.35926 arg_L[c] + 0.27942 asn_L[c] + 0.35261 asp_L[c] + 20.7045 atp[c] + 0.020401 chsterol[c] + 0.011658 clpn_hs[c] + 0.039036 ctp[c] + 0.046571 cys_L[c] + 0.27519 g6p[c] + 0.326 gln_L[c] + 0.38587 glu_L[c] + 0.53889 gly[c] + 0.036117 gtp[c] + 20.6508 h2o[c] + 0.12641 his_L[c] + 0.28608 ile_L[c] + 0.54554 leu_L[c] + 0.59211 lys_L[c] + 0.15302 met_L[c] + 0.023315 pail_hs[c] + 0.15446 pchol_hs[c] + 0.055374 pe_hs[c] + 0.002914 pglyc_hs[c] + 0.25947 phe_L[c] + 0.41248 pro_L[c] + 0.005829 ps_hs[c] + 0.39253 ser_L[c] + 0.017486 sphmyln_hs[c] + 0.31269 thr_L[c] + 0.013306 trp_L[c] + 0.15967 tyr_L[c] + 0.053446 utp[c] + 0.35261 val_L[c] -> 20.6508 adp[c] + 20.6508 h[c] + 20.6508 pi[c]', 'subSystem','Exchange/demand reaction', 'geneRule','');
model = addReaction(model,'biomass_maintenance_noTrTr','reactionFormula','20.7045 atp[c] + 0.020401 chsterol[c] + 0.011658 clpn_hs[c] + 0.27519 g6p[c] + 20.6508 h2o[c] + 0.023315 pail_hs[c] + 0.15446 pchol_hs[c] + 0.055374 pe_hs[c] + 0.002914 pglyc_hs[c] + 0.005829 ps_hs[c] + 0.017486 sphmyln_hs[c] -> 20.6508 adp[c] + 20.6508 h[c] + 20.6508 pi[c]', 'subSystem','Exchange/demand reaction', 'geneRule','');
model = changeObjective(model, 'biomass_maintenance', 0);

%% add essential rxns
model= addExchangeRxn(model,{'glc_D[e]'},-10,1000);
model = addReaction(model,'DM_atp_c_','reactionFormula','atp[c] + h2o[c] -> adp[c] + h[c] + pi[c]', 'subSystem','Exchange/demand reaction', 'geneRule','');
model = addExchangeRxn(model,'o2[e]', -10, 0);
model = addExchangeRxn(model,'co2[e]', 0, 1000);
model = addExchangeRxn(model,{'4mop[e]'}, -10, 1000);
model=changeRxnBounds(model,'EX_o2[e]',0,'u');
model=changeRxnBounds(model,'EX_co2[e]',0,'l');

%% add literature curated rxns
model = addReaction(model,'GTHS','reactionFormula','atp[c] + glucys[c] + gly[c] -> adp[c] + gthrd[c] + h[c] + pi[c]', 'subSystem','Glutathione metabolism', 'geneRule','2937.1'); %gluthathione be secreted in to media
model = addReaction(model,'CO2t','reactionFormula','co2[e] <=> co2[c]', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'PGI','reactionFormula','g6p[c] <=> f6p[c]', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','2821.1');
model = addReaction(model,'TPI','reactionFormula','dhap[c] <=> g3p[c]', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','7167.1 or 286016.1');
model = addReaction(model,'PGM','reactionFormula','2pg[c] <=> 3pg[c]', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','669.1 or 5223.1 or 5224.2 or 5224.1 or 669.2');
model = addReaction(model,'PYK','reactionFormula','adp[c] + h[c] + pep[c] -> atp[c] + pyr[c]', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','5315.2 or 5313.2 or 5315.3 or 5313.1 or 5315.1');
model = addReaction(model,'ME2','reactionFormula','mal_L[c] + nadp[c] -> co2[c] + nadph[c] + pyr[c]', 'subSystem','Pyruvate metabolism', 'geneRule','4199.1');
model = addReaction(model,'PCm','reactionFormula','atp[m] + hco3[m] + pyr[m] -> adp[m] + h[m] + oaa[m] + pi[m]', 'subSystem','Pyruvate metabolism', 'geneRule','5091.1 or 5091.2');
model = addReaction(model,'r0941','reactionFormula','hco3[c] <=> hco3[m]', 'subSystem','Transport, mitochondrial', 'geneRule','');
model = addReaction(model,'PYRt2m','reactionFormula','h[c] + pyr[c] <=> h[m] + pyr[m]', 'subSystem','Transport, mitochondrial', 'geneRule','6566.1');
model = addReaction(model,'ICDHxm','reactionFormula','icit[m] + nad[m] -> akg[m] + co2[m] + nadh[m]', 'subSystem','Citric acid cycle', 'geneRule','(3421.2 and 3420.3 and 3419.1) or (3420.1 and 3421.1 and 3419.1) or (3420.1 and 3421.1 and 3419.1) or (3420.3 and 3421.1 and 3419.1) or (3421.1 and 3419.1 and 3420.2) or (3421.2 and 3419.1 and 3420.2)');
model = addReaction(model,'ICDHy','reactionFormula','icit[c] + nadp[c] -> akg[c] + co2[c] + nadph[c]', 'subSystem','Citric acid cycle', 'geneRule','3417.1');
model = addReaction(model,'SUCD1m','reactionFormula','fad[m] + succ[m] <=> fadh2[m] + fum[m]', 'subSystem','Citric acid cycle', 'geneRule','6389.1 and 6392.1 and 6391.1 and 6390.1');
model = addReaction(model,'CITtam','reactionFormula','cit[c] + mal_L[m] <=> cit[m] + mal_L[c]', 'subSystem','Transport, mitochondrial', 'geneRule','6576.1');
model = addReaction(model,'r0081','reactionFormula','akg[m] + ala_L[m] <=> glu_L[m] + pyr[m]', 'subSystem','Citric acid cycle', 'geneRule','84706.1 or 2875.1');
model = addReaction(model,'ILETA','reactionFormula','akg[c] + ile_L[c] <=> 3mop[c] + glu_L[c]', 'subSystem','Valine, leucine, and isoleucine metabolism', 'geneRule','586.1');
model = addReaction(model,'LEUTA','reactionFormula','akg[c] + leu_L[c] <=> 4mop[c] + glu_L[c]', 'subSystem','Valine, leucine, and isoleucine metabolism', 'geneRule','586.1');
model = addReaction(model,'4MOPte','reactionFormula','4mop[e] <=> 4mop[c]', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'3SALACBOXL','reactionFormula','3sala[c] + 2.0 h[c] -> co2[c] + hyptaur[c]', 'subSystem','Taurine and hypotaurine metabolism', 'geneRule','51380.1 or 2571.1 or 2572.1 or 2571.2');
model = addReaction(model,'HYPTROX','reactionFormula','2.0 hyptaur[c] + o2[c] -> 2.0 taur[c]', 'subSystem','Taurine and hypotaurine metabolism', 'geneRule','');
model = addReaction(model,'SRR','reactionFormula','ser_L[c] <=> ser_D[c]', 'subSystem','', 'geneRule','');  %not in Recon3D
model = addReaction(model,'CYSTS','reactionFormula','hcys_L[c] + ser_L[c] -> cyst_L[c] + h2o[c]', 'subSystem','Methionine and cysteine metabolism', 'geneRule','875.1');
model = addReaction(model,'ILETAm','reactionFormula','akg[m] + ile_L[m] <=> 3mop[m] + glu_L[m] ', 'subSystem','Valine, leucine, and isoleucine metabolism', 'geneRule','587.1');
model = addReaction(model,'LEUTAm','reactionFormula','akg[m] + leu_L[m] <=> 4mop[m] + glu_L[m]', 'subSystem','Valine, leucine, and isoleucine metabolism', 'geneRule','587.1');
model = addReaction(model,'ASPTA','reactionFormula','akg[c] + asp_L[c] <=> glu_L[c] + oaa[c]', 'subSystem','Alanine and aspartate metabolism', 'geneRule','2805.1');
model = addReaction(model,'ALATA_L','reactionFormula','akg[c] + ala_L[c] <=> glu_L[c] + pyr[c]', 'subSystem','Glutamate metabolism', 'geneRule','2875.1 or 84706.1');
model = addReaction(model,'GLUDxm','reactionFormula','glu_L[m] + h2o[m] + nad[m] <=> akg[m] + h[m] + nadh[m] + nh4[m]', 'subSystem','Glutamate metabolism', 'geneRule','2747.1 or 2746.1');
model = addReaction(model,'HMR_9802','reactionFormula','gln_L[c] + h2o[c] -> glu_L[c] + nh4[c]', 'subSystem','Glutamate metabolism', 'geneRule','27165.1');
model = addReaction(model,'GLUDC','reactionFormula','glu_L[c] + h[c] -> 4abut[c] + co2[c]', 'subSystem','Glutamate metabolism', 'geneRule','2571.1 or 2572.1 or 2571.2');
model = addReaction(model,'HISDC','reactionFormula','h[c] + his_L[c] -> co2[c] + hista[c]', 'subSystem','Histidine metabolism', 'geneRule','1644.1 or 3067.1');
model = addReaction(model,'GLNS','reactionFormula','atp[c] + glu_L[c] + nh4[c] -> adp[c] + gln_L[c] + h[c] + pi[c]', 'subSystem','Glutamate metabolism', 'geneRule','2752.1 or 51557.1');
model = addReaction(model,'PDHm','reactionFormula','coa[m] + nad[m] + pyr[m] -> accoa[m] + co2[m] + nadh[m]', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','(1738.1 and 8050.1) and (5161.1 and 5162.1) and (1737.1) or (1738.1 and 8050.1) and (5160.1 and 5162.1) and (1737.1)');
model = addReaction(model,'LDH_L','reactionFormula','lac_L[c] + nad[c] <=> h[c] + nadh[c] + pyr[c]', 'subSystem','Glycolysis/gluconeogenesis', 'geneRule','(3945.1 and 3939.1) or 160287.1 or 3948.2 or 3939.1 or 3948.1 or 55293.1 or 3945.1 or 92483.1');
model = addReaction(model,'LDH_Lm','reactionFormula','lac_L[m] + nad[m] <=> h[m] + nadh[m] + pyr[m]', 'subSystem','Pyruvate metabolism', 'geneRule','3939.1 or 3945.1 or 197257.1');
model = addReaction(model,'ABTArm','reactionFormula','4abut[m] + akg[m] <=> glu_L[m] + sucsal[m]', 'subSystem','Glutamate metabolism', 'geneRule','18.1 or 18.2');
model = addReaction(model,'METAT','reactionFormula','atp[c] + h2o[c] + met_L[c] -> amet[c] + pi[c] + ppi[c]', 'subSystem','Methionine and cysteine metabolism', 'geneRule','4143.1 or (27430.2 and 4144.1) or (27430.1 and 4144.1)');
model = addReaction(model,'ADNK1','reactionFormula','adn[c] + atp[c] -> adp[c] + amp[c] + h[c]', 'subSystem','Nucleotide interconversion', 'geneRule','132.1 or 132.2');
model = addReaction(model,'ADNK1m','reactionFormula','adn[m] + atp[m] -> adp[m] + amp[m] + h[m]', 'subSystem','Nucleotide interconversion', 'geneRule','');
model = addReaction(model,'HXPRT','reactionFormula','hxan[c] + prpp[c] -> imp[c] + ppi[c]', 'subSystem','Nucleotide salvage pathway', 'geneRule','3251.1');
model = addReaction(model,'GUAPRT','reactionFormula','gua[c] + prpp[c] -> gmp[c] + ppi[c]', 'subSystem','Nucleotide salvage pathway', 'geneRule','3251.1');
model = addReaction(model,'BDHm','reactionFormula','bhb[m] + nad[m] <=> acac[m] + h[m] + nadh[m]', 'subSystem','Nucleotide salvage pathway', 'geneRule','622.2 or 622.3 or 622.1');
model = addReaction(model,'ME1m','reactionFormula','mal_L[m] + nad[m] -> co2[m] + nadh[m] + pyr[m]', 'subSystem','Pyruvate metabolism', 'geneRule','4200.1');
model = addReaction(model,'CK','reactionFormula','atp[m] + creat[m] <=> adp[m] + h[m] + pcreat[m]', 'subSystem','Urea cycle', 'geneRule','1159.1 or 1160.1 or 548596.1');
model = addReaction(model,'PCREATtmdiffir','reactionFormula','pcreat[m] -> pcreat[c]', 'subSystem','Urea cycle', 'geneRule','');
model = addReaction(model,'CREATtmdiffir','reactionFormula','creat[c] -> creat[m]', 'subSystem','Urea cycle', 'geneRule','');
model = addReaction(model,'PGS','reactionFormula','arachd[c] + h[c] + nadph[c] + 2.0 o2[c] -> h2o[c] + nadp[c] + prostgh2[c]', 'subSystem','Eicosanoid metabolism', 'geneRule','5743.1 or 5742.1');
model = addReaction(model,'TAURt4_2_r','reactionFormula','2.0 na1[e] + taur[e] <=> 2.0 na1[c] + taur[c] ', 'subSystem','Transport, extracellular', 'geneRule','6533.1');
model = addReaction(model,'ACITL','reactionFormula','atp[c] + cit[c] + coa[c] -> accoa[c] + adp[c] + oaa[c] + pi[c]', 'subSystem','Citric acid cycle', 'geneRule','47.1 or 47.2');
model = addReaction(model,'CHOLt4','reactionFormula','atp[c] + cit[c] + coa[c] -> accoa[c] + adp[c] + oaa[c] + pi[c]', 'subSystem','Citric acid cycle', 'geneRule','47.1 or 47.2');
model = addReaction(model,'DASCBR','reactionFormula','dhdascb[c] + nadph[c] -> ascb_L[c] + nadp[c]', 'subSystem','Vitamin C metabolism', 'geneRule','2745.1 or 51022.2 or 51022.1');
model = addReaction(model,'LYSt4','reactionFormula','lys_L[e] + na1[e] -> lys_L[c] + na1[c]', 'subSystem','Transport, extracellular', 'geneRule','11254.1 or 6584.1'); % wasn't captured by iMAT in BD condition
model = addReaction(model,'MDHm','reactionFormula','mal_L[m] + nad[m] <=> h[m] + nadh[m] + oaa[m]', 'subSystem','Citric acid cycle', 'geneRule','4191.1 or 4190.1'); % wasn't captured by iMAT in BD condition
model = addReaction(model,'THMtrbc','reactionFormula','thm[e] <=> thm[c]', 'subSystem','Transport, extracellular', 'geneRule',''); % wasn't captured by iMAT in BD condition

% literature curated astrocytic media secretions:
model = addExchangeRxn(model,{'taur[e]'},0,1000); % secretion of taurine into culture media
model = addExchangeRxn(model,{'ser_D[e]'},0,1000); % to gapFill d-serine into culture media
model = addExchangeRxn(model,{'cit[e]'},0,1000); % secretion of cit into culture media
model = addExchangeRxn(model,{'atp[e]'},0,1000); % secretion of atp into culture media
model = addExchangeRxn(model,{'gln_L[e]'},0,1000); % secretion of glutamine into culture media
model = addExchangeRxn(model,{'lac_L[e]'},0,1000); % secretion of lactate into culture media
model = addExchangeRxn(model,{'ala_L[e]'},0,1000); % secretion of alanine into culture media
model = addExchangeRxn(model,{'gthrd[e]'},0,1000); % secretion of glutathione into culture media
model = addExchangeRxn(model,{'chsterol[e]'},0,1000); % secretion of cholesterol into culture media
model = addExchangeRxn(model,{'HC00250[e]'},0,1000); % secretion of hydrogen sulphide into culture media

% gapFill so that some of the above rxns carry a non-zero flux:
model = addSinkReactions(model,{'pi[m]'},-10,1000); % to gapFill PCm
model = addSinkReactions(model,{'bhb[m]'},-10,1000); % to gapFill BDHm
model = addSinkReactions(model,{'acac[m]'},-10,1000); % to gapFill BDHm
model = addSinkReactions(model,{'arachd[c]'},-10,1000); % to gapFill PGS
model = addSinkReactions(model,{'prostgh2[c]'},-10,1000); % to gapFill PGS
model = addReaction(model,'HMR_9797','reactionFormula','ade[c] + h[c] + h2o[c] -> hxan[c] + nh4[c]', 'subSystem','Purine catabolism', 'geneRule',''); % to gapFill HXPRT
model = addReaction(model,'GACMTRc','reactionFormula','amet[c] + gudac[c] <=> ahcys[c] + creat[c] + h[c] ', 'subSystem','Urea cycle', 'geneRule','2593.1 or 2593.2'); % to gapFill CK
model = addReaction(model,'HMR_9191_e','reactionFormula','ser_D[c] -> ser_D[e]', 'subSystem','Transport, extracellular', 'geneRule','56301.1 or 6520.1'); % not in recon; to gapFill SRR
model = addReaction(model,'L_LACtm','reactionFormula','h[c] + lac_L[c] -> h[m] + lac_L[m]', 'subSystem','Transport, mitochondrial', 'geneRule','6566.1'); % to gapFill LDH_Lm
model = addReaction(model,'FE3R2e','reactionFormula','ascb_L[e] + 2.0 fe3[e] -> dhdascb[e] + 2.0 fe2[e] + h[e] ', 'subSystem','Heme synthesis', 'geneRule','79901.1'); % to gapFill EX_fe3[e]
model = addReaction(model,'r1106','reactionFormula','ribflv[e] <=> ribflv[c]', 'subSystem','Transport, extracellular', 'geneRule',''); % to gapFill EX_ribflv[e]
model = addReaction(model,'TMDPK','reactionFormula','atp[c] + thm[c] -> amp[c] + h[c] + thmpp[c]', 'subSystem','Thiamine metabolism', 'geneRule','27010.1'); % to gapFill EX_thm[e]
model = addReaction(model,'PRGSTRNt','reactionFormula','prgstrn[e] <=> prgstrn[c]', 'subSystem','Transport, extracellular', 'geneRule',''); % to gapFill EX_prgstrn[e]
model = addReaction(model,'HMR_1981','reactionFormula','h[c] + nadph[c] + o2[c] + prgstrn[c] -> 17ahprgstrn[c] + h2o[c] + nadp[c]', 'subSystem','Steroid metabolism', 'geneRule','1586.1'); % to gapFill EX_prgstrn[e]
model = addReaction(model,'HMR_2002','reactionFormula','C03681[c] + nadp[c] <=> h[c] + nadph[c] + prgstrn[c]', 'subSystem','Steroid metabolism', 'geneRule',''); % to gapFill EX_prgstrn[e]
model = addReaction(model,'HMR_2007','reactionFormula','h[c] + nadph[c] + o2[c] + prgstrn[c] -> 11docrtstrn[c] + h2o[c] + nadp[c]', 'subSystem','Steroid metabolism', 'geneRule','1589.1 or 1589.1 or 1589.1'); % to gapFill EX_prgstrn[e]
model = addReaction(model,'17AHPRGSTRNte','reactionFormula','17ahprgstrn[c] + atp[c] + h2o[c] -> 17ahprgstrn[e] + adp[c] + h[c] + pi[c]', 'subSystem','Transport, extracellular', 'geneRule','5243.1'); % to gapFill EX_prgstrn[e]
model = addReaction(model,'HMR_1988','reactionFormula','17ahprgstrn[c] + h[c] + nadph[c] + o2[c] -> 11docrtsl[c] + h2o[c] + nadp[c]', 'subSystem','Steroid metabolism', 'geneRule','1589.1 or 1589.1 or 1589.1'); % to gapFill EX_prgstrn[e]
model = addReaction(model,'11DOCRTSLte','reactionFormula','11docrtsl[c] + atp[c] + h2o[c] -> 11docrtsl[e] + adp[c] + h[c] + pi[c]', 'subSystem','Transport, extracellular', 'geneRule','5243.1'); % to gapFill EX_prgstrn[e]
model = addSinkReactions(model,{'11docrtsl[c]'},-10,1000); % to gapFill EX_prgstrn[e]
model = addSinkReactions(model,{'17ahprgstrn[c]'},-10,1000); % to gapFill EX_prgstrn[e]
model = addReaction(model,'NCAMUP','reactionFormula','ncam[e] -> ncam[c]', 'subSystem','Transport, extracellular', 'geneRule',''); % to gapFill EX_ncam[e]
model = addReaction(model,'NMNS','reactionFormula','h[c] + ncam[c] + prpp[c] -> nmn[c] + ppi[c]', 'subSystem','NAD metabolism', 'geneRule','10135.1'); % to gapFill EX_ncam[e]
model = addReaction(model,'CHOLtu','reactionFormula','chol[e] <=> chol[c]', 'subSystem','Transport, extracellular', 'geneRule','6582.1 or 6584.1 or 6582.2'); % to gapFill EX_chol[e] in BD model

% add PI cycle (Saiardi'18) missing rxns:
model = addExchangeRxn(model,{'inost[e]'},-10,1000);
model = addReaction(model,'MI1PP','reactionFormula','h2o[c] + mi1p_D[c] -> inost[c] + pi[c] ', 'subSystem','Inositol phosphate metabolism', 'geneRule', '3612.1 or 3613.1 or 8776.1 or 8898.1');
model = addReaction(model,'INSTt4','reactionFormula','inost[e] + na1[e] <=> inost[c] + na1[c]', 'subSystem','Transport, extracellular', 'geneRule','6526.1');
model = addReaction(model,'DAGK_hs','reactionFormula','atp[c] + dag_hs[c] <=> adp[c] + h[c] + pa_hs[c]', 'subSystem','Glycerophospholipid metabolism', 'geneRule','9162.1 or 1607.1 or 8525.2 or 160851.1 or 1608.1 or 1607.2 or 8525.1 or 8527.1 or 1606.1 or 160851.2 or 8526.1 or 8525.3 or 1607.1 or 8527.2 or 1609.1');
model = addReaction(model,'MI1PS','reactionFormula','g6p[c] -> mi1p_D[c]', 'subSystem','Inositol phosphate metabolism', 'geneRule','51477.1');
model = addReaction(model,'MI14P4P','reactionFormula','h2o[c] + mi14p[c] -> mi1p_D[c] + pi[c]', 'subSystem','Inositol phosphate metabolism', 'geneRule','51477.1');
model = addReaction(model,'MI1345PKc','reactionFormula','atp[c] + mi1345p[c] -> adp[c] + h[c] + mi13456p[c]', 'subSystem','Inositol phosphate metabolism', 'geneRule','253430.1');
model = addReaction(model,'MI13456PK','reactionFormula','atp[c] + mi13456p[c] -> adp[c] + h[c] + minohp[c]', 'subSystem','Inositol phosphate metabolism', 'geneRule','');
model = addSinkReactions(model,{'minohp[c]'},-10,1000);

%% add ExchangeRxns as per media composition..
% add ASM_media_EX_unconstrained:
model=addExchangeRxn(model,{'ala_L[e]'},-10,1000);
model=addExchangeRxn(model,{'gln_L[e]'},-10,1000);
model=addExchangeRxn(model,{'arg_L[e]'},-10,1000);
model=addExchangeRxn(model,{'asn_L[e]'},-10,1000);
model=addExchangeRxn(model,{'asp_L[e]'},-10,1000);
model=addExchangeRxn(model,{'ca2[e]'},-10,1000);
model=addExchangeRxn(model,{'chol[e]'},-10,1000);
model=addExchangeRxn(model,{'cl[e]'},-10,1000);
model=addExchangeRxn(model,{'cys_L[e]'},-10,1000);
model=addExchangeRxn(model,{'fe3[e]'},-10,1000);
model=addExchangeRxn(model,{'fol[e]'},-10,1000);
model=addExchangeRxn(model,{'glc_D[e]'},-10,1000);
model=addExchangeRxn(model,{'glu_L[e]'},-10,1000);
model=addExchangeRxn(model,{'gly[e]'},-10,1000);
model=addExchangeRxn(model,{'hco3[e]'},-10,1000);
model=addExchangeRxn(model,{'his_L[e]'},-10,1000);
model=addExchangeRxn(model,{'ile_L[e]'},-10,1000);
model=addExchangeRxn(model,{'inost[e]'},-10,1000);
model=addExchangeRxn(model,{'k[e]'},-10,1000);
model=addExchangeRxn(model,{'leu_L[e]'},-10,1000);
model=addExchangeRxn(model,{'lys_L[e]'},-10,1000);
model=addExchangeRxn(model,{'met_L[e]'},-10,1000);
model=addExchangeRxn(model,{'mg2[e]'},-10,1000);
model=addExchangeRxn(model,{'na1[e]'},-10,1000);
model=addExchangeRxn(model,{'ncam[e]'},-10,1000);
model=addExchangeRxn(model,{'no3[e]'},-10,1000);
model=addExchangeRxn(model,{'phe_L[e]'},-10,1000);
model=addExchangeRxn(model,{'pi[e]'},-10,1000);
model=addExchangeRxn(model,{'pnto_R[e]'},-10,1000);
model=addExchangeRxn(model,{'prgstrn[e]'},-10,1000);
model=addExchangeRxn(model,{'pro_L[e]'},-10,1000);
model=addExchangeRxn(model,{'ptrc[e]'},-10,1000);
model=addExchangeRxn(model,{'pydx[e]'},-10,1000);
model=addExchangeRxn(model,{'pyr[e]'},-10,1000);
model=addExchangeRxn(model,{'ribflv[e]'},-10,1000);
model=addExchangeRxn(model,{'selni[e]'},-10,1000);
model=addExchangeRxn(model,{'ser_L[e]'},-10,1000);
model=addExchangeRxn(model,{'so4[e]'},-10,1000);
model=addExchangeRxn(model,{'thm[e]'},-10,1000);
model=addExchangeRxn(model,{'thr_L[e]'},-10,1000);
model=addExchangeRxn(model,{'trp_L[e]'},-10,1000);
model=addExchangeRxn(model,{'tyr_L[e]'},-10,1000);
model=addExchangeRxn(model,{'val_L[e]'},-10,1000);
model=addExchangeRxn(model,{'zn2[e]'},-10,1000);

% addition of misc. reactions for biomass_maintenance.
model = changeRxnBounds(model,'EX_o2[e]',-10,'l');
model = addSinkReactions(model,{'chsterol[c]'},-10,1000);
model = addSinkReactions(model,{'clpn_hs[c]'},-10,1000);
model = addSinkReactions(model,{'ctp[c]'},-10,1000);
model = addSinkReactions(model,{'gtp[c]'},-10,1000);
model = addSinkReactions(model,{'pail_hs[c]'},-10,1000);
model = addSinkReactions(model,{'pe_hs[c]'},-10,1000);
model = addSinkReactions(model,{'pglyc_hs[c]'},-10,1000);
model = addSinkReactions(model,{'ps_hs[c]'},-10,1000);
model = addSinkReactions(model,{'sphmyln_hs[c]'},-10,1000);
model = addSinkReactions(model,{'utp[c]'},-10,1000);

% add BBB_media_EX_unconstrained:
model=addExchangeRxn(model,{'glu_L[e]'},-10,1000);
model=addExchangeRxn(model,{'crn[e]'},-10,1000);
model=addExchangeRxn(model,{'lac_L[e]'},-10,1000);
model=addExchangeRxn(model,{'lys_L[e]'},-10,1000);
model=addExchangeRxn(model,{'arg_L[e]'},-10,1000);
model=addExchangeRxn(model,{'orn[e]'},-10,1000);
model=addExchangeRxn(model,{'his_L[e]'},-10,1000);
model=addExchangeRxn(model,{'gln_L[e]'},-10,1000);
model=addExchangeRxn(model,{'met_L[e]'},-10,1000);
model=addExchangeRxn(model,{'leu_L[e]'},-10,1000);
model=addExchangeRxn(model,{'ile_L[e]'},-10,1000);
model=addExchangeRxn(model,{'val_L[e]'},-10,1000);
model=addExchangeRxn(model,{'phe_L[e]'},-10,1000);
model=addExchangeRxn(model,{'trp_L[e]'},-10,1000);
model=addExchangeRxn(model,{'cys_L[e]'},-10,1000);
model=addExchangeRxn(model,{'asn_L[e]'},-10,1000);
model=addExchangeRxn(model,{'ala_L[e]'},-10,1000);
model=addExchangeRxn(model,{'glc_D[e]'},-10,1000);
model=addExchangeRxn(model,{'gal[e]'},-10,1000);
model=addExchangeRxn(model,{'dhdascb[e]'},-10,1000);
model=addExchangeRxn(model,{'k[e]'},-10,1000);
model=addExchangeRxn(model,{'chol[e]'},-10,1000);
model=addExchangeRxn(model,{'ade[e]'},-10,1000);
model=addExchangeRxn(model,{'adn[e]'},-10,1000);
model=addExchangeRxn(model,{'triodthy[e]'},-10,1000);
model=addExchangeRxn(model,{'bhb[e]'},-10,1000);
model=addExchangeRxn(model,{'thyox_L[e]'},-10,1000);
model=addExchangeRxn(model,{'tststerone[e]'},-10,1000);
model=addExchangeRxn(model,{'estradiol[e]'},-10,1000);
model=addExchangeRxn(model,{'crtsl[e]'},-10,1000);
model=addExchangeRxn(model,{'aldstrn[e]'},-10,1000);
model=addExchangeRxn(model,{'crtstrn[e]'},-10,1000);
model=addExchangeRxn(model,{'cortsn[e]'},-10,1000);
model=addExchangeRxn(model,{'prgstrn[e]'},-10,1000);
model=addExchangeRxn(model,{'melatn[e]'},-10,1000);
model=addExchangeRxn(model,{'hdca[e]'},-10,1000);
model=addExchangeRxn(model,{'lnlc[e]'},-10,1000);
model=addExchangeRxn(model,{'octa[e]'},-10,1000);
model=addExchangeRxn(model,{'but[e]'},-10,1000);
model=addExchangeRxn(model,{'pnto_R[e]'},-10,1000);
model=addExchangeRxn(model,{'ttdca[e]'},-10,1000);
model=addExchangeRxn(model,{'tyr_L[e]'},-10,1000);
model=addExchangeRxn(model,{'asp_L[e]'},-10,1000);
model=addExchangeRxn(model,{'taur[e]'},-10,1000);
model=addExchangeRxn(model,{'ser_L[e]'},-10,1000);

%% save non-media model
model_EXP_Unconstrained = model;
%%

%% ASTROCYTE SUSTENANCE MEDIA_START
%% define closed model
modelClosed = model;
modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
modelexchanges = unique([modelexchanges1; modelexchanges4]);
modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
model = modelClosed;

%% now update ASM bounds, with lb to 10:
model=changeRxnBounds(model,'EX_ala_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_gln_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_arg_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_asn_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_asp_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_ca2[e]',-10,'l');
model=changeRxnBounds(model,'EX_chol[e]',-10,'l');
model=changeRxnBounds(model,'EX_cl[e]',-10,'l');
model=changeRxnBounds(model,'EX_cys_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_fe3[e]',-10,'l');
model=changeRxnBounds(model,'EX_fol[e]',-10,'l');
model=changeRxnBounds(model,'EX_glc_D[e]',-10,'l');
model=changeRxnBounds(model,'EX_glu_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_gly[e]',-10,'l');
model=changeRxnBounds(model,'EX_hco3[e]',-10,'l');
model=changeRxnBounds(model,'EX_his_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_ile_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_inost[e]',-10,'l');
model=changeRxnBounds(model,'EX_k[e]',-10,'l');
model=changeRxnBounds(model,'EX_leu_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_lys_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_met_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_mg2[e]',-10,'l');
model=changeRxnBounds(model,'EX_na1[e]',-10,'l');
model=changeRxnBounds(model,'EX_ncam[e]',-10,'l');
model=changeRxnBounds(model,'EX_no3[e]',-10,'l');
model=changeRxnBounds(model,'EX_phe_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_pi[e]',-10,'l');
model=changeRxnBounds(model,'EX_pnto_R[e]',-10,'l');
model=changeRxnBounds(model,'EX_prgstrn[e]',-10,'l');
model=changeRxnBounds(model,'EX_pro_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_ptrc[e]',-10,'l');
model=changeRxnBounds(model,'EX_pydx[e]',-10,'l');
model=changeRxnBounds(model,'EX_pyr[e]',-10,'l');
model=changeRxnBounds(model,'EX_ribflv[e]',-10,'l');
model=changeRxnBounds(model,'EX_selni[e]',-10,'l');
model=changeRxnBounds(model,'EX_ser_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_so4[e]',-10,'l');
model=changeRxnBounds(model,'EX_thm[e]',-10,'l');
model=changeRxnBounds(model,'EX_thr_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_trp_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_tyr_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_val_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_zn2[e]',-10,'l');
%% ASTROCYTE SUSTENANCE MEDIA_END

%% allow o2.. 
model = changeRxnBounds(model,'EX_o2[e]',-10,'l');
%%
%% check for blocked reactions
fluxConsistency = verifyModel(model,'fluxConsistency',true); fluxInconsistentRxns_media = model.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1); clear fluxConsistency;
fluxConsistency = verifyModel(model_EXP_Unconstrained,'fluxConsistency',true); fluxInconsistentRxns_non_media = model_EXP_Unconstrained.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1); clear fluxConsistency;

%% save media model
model = changeObjective(model, 'biomass_maintenance', 0);
model_EXP_Constrained = model;

%%
toc;
end
