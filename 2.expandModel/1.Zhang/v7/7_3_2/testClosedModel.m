function [TestSolution,TestSolutionName,TestedRxns,PercTestedRxns] = test4HumanFctExt(model,test,optionSinks)
changeCobraSolver('gurobi','all');

% fixes the model met names in cases the compartments are not given with ()
model.mets = regexprep(model.mets,'[','(');
model.mets = regexprep(model.mets,']',')');
model.mets = regexprep(model.mets,'_','-');
model.mets = regexprep(model.mets,'-FSLASH-','/');
model.rxns = regexprep(model.rxns,'\[','\(');
model.rxns = regexprep(model.rxns,'\]','\)');

%define closed model
modelClosed = model;
modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
modelexchanges2 = strmatch('DM_',modelClosed.rxns);
modelexchanges3 = strmatch('sink_',modelClosed.rxns);
selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
%set bounds of all exchanges to: lb=0 & ub=1000
modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=1000;
model = modelClosed;

model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=0;

TestSolution = [];
TestedRxns =[];
k = 1;
tol = 1e-6;
% RPMI_composition={'EX_ala_L(e)','EX_arg-L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_cys-L(e)','EX_gln-L(e)','EX_glu-L(e)','EX_gly(e)','EX_his-L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys-L(e)','EX_met_L(e)','EX_phe_L(e)','EX_4HPRO','EX_pro-L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_ascb_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)','EX_fol(e)','EX_ncam(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_cbl1(e)','EX_inost(e)','EX_ca2(e)','EX_fe3(e)','EX_k(e)','EX_hco3(e)','EX_na1(e)','EX_pi(e)','EX_glc(e)','EX_hxan(e)','EX_lnlc(e)','EX_lipoate(e)','EX_ptrc(e)','EX_pyr(e)','EX_thymd(e)','EX_etha(e)','EX_gthrd(e)'};
modelOri = model;
%% "Human Recon 1 human  biomass"
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'biomass_maintenance'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Human Recon 1 test human biomass (maintenance)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
%% gln-L(c) -> nh4(c)/ glutaminase
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gln-L(c)','nh4(c)'},[-1000 -1; 0 1000]);
    model.lb(ismember(model.rxns,'EX_nh4(e)'))=0; model.ub(ismember(model.rxns,'EX_nh4(e)'))=1000;
    model.c(ismember(model.rxns,'sink_nh4(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gln-L(c) -> nh4(c)/ glutaminase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

TestSolutionName(:,2) = num2cell(TestSolution);
TestedRxns = unique(TestedRxns);
TestedRxns = intersect(model.rxns,TestedRxns); % only those reactions that are also in modelOri not those that have been added to the network
PercTestedRxns = length(TestedRxns)*100/length(model.rxns);