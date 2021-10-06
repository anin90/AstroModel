function [TestResults] = test4HumanAstro(model)

changeCobraSolver('gurobi','all'); 

%% keep original model as 'modelOri'
modelOri = model;

%% define closed model
modelClosed = model;
modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges2 = strmatch('DM_',modelClosed.rxns);
modelexchanges3 = strmatch('sink_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
selExc = (find( full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))))';
modelexchanges = unique([modelexchanges1; modelexchanges2; modelexchanges3; modelexchanges4; selExc]);
modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;

% allow oxygen and clear objectives
modelClosed = changeRxnBounds(modelClosed,'EX_o2[e]',-10,'l');
modelClosed.c(find(modelClosed.c)) = 0;
k = 1;

%% ASTROCYTE_FUNCTION_TESTS

    %% "Recon3D  biomass_maintenance"
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'biomass_maintenance'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'biomass_maintenance';
    TestedModel{k,1} = 'modelOri';
    k = k +1;
    %% "Recon3D  biomass_maintenance_noTrTr"
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'biomass_maintenance_noTrTr'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'biomass_maintenance_noTrTr';
    TestedModel{k,1} = 'modelOri';
    k = k +1;
    %% ATP max aerobic, glc, v0.05
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.rxns,'EX_glc_D[e]')))
		model.lb(ismember(model.rxns,'EX_glc_D[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc_D[e]'))=-1;
		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, glc -> atp';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;
    %% ATP max, anaerobic glc, v0.05
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.rxns,'EX_glc_D[e]')))
		model.lb(ismember(model.rxns,'EX_glc_D[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc_D[e]'))=-1;
		model.lb(ismember(model.rxns,'EX_o2[e]'))=0;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, anaerobic, glc -> atp';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;   
  %% ATP max, aerobic, pyr -> atp
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.rxns,'EX_pyr[e]')))
		model.lb(ismember(model.rxns,'EX_pyr[e]'))=-1;model.ub(ismember(model.rxns,'EX_pyr[e]'))=-1;
		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, pyr -> atp';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;
 %% ATP max, aerobic, lac_L -> atp
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.rxns,'EX_lac_L[e]')))
		model.lb(ismember(model.rxns,'EX_lac_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_lac_L[e]'))=-1;
		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, lac_L -> atp';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;       
 %% ATP max, aerobic, glu_L, synaptic -> atp
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.rxns,'EX_glu_L[s]')))
		model.lb(ismember(model.rxns,'EX_glu_L[s]'))=-1;model.ub(ismember(model.rxns,'EX_glu_L[s]'))=-1;
		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, glu_L (synaptic) -> atp';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;
 %% ATP max, aerobic, glu_L, extracellular -> atp
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.rxns,'EX_glu_L[e]')))
		model.lb(ismember(model.rxns,'EX_glu_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_glu_L[e]'))=-1;
		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, glu_L (extracellular) -> atp';
    TestedModel{k,1} = 'modelClosed';
    k = k +1; 
%% pyr[c] -> accoa[m]
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'accoa[m]')))
        [model] = addSinkReactions(model,{'pyr[c]','accoa[m]'},[-1000 -1; 0 1000]);
        model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
        FBA = optimizeCbModel(model);
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> accoa[m]';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;    
    %% pyr[c] -> oaa[m]
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'oaa[m]')))
		[model] = addSinkReactions(model,{'pyr[c]','oaa[m]'},[-1000 -1; 0 1000]);
		model.c(ismember(model.rxns,'sink_oaa[m]'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
	else
		TestSolution(k,1) = NaN;
	end	
    TestSolutionName{k,1} = 'pyr[c] -> oaa[m]';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;    
    %% pyr[c] -> etoh[c]
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'etoh[c]')))
		[model] = addSinkReactions(model,{'pyr[c]','etoh[c]'},[-1000 -1; 0 1000]);
		model.c(ismember(model.rxns,'sink_etoh[c]'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
	end	
    TestSolutionName{k,1} = 'pyr[c] -> etoh[c]';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;
    %% pyr[c] -> lac_L[s]
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'lac_L[s]')))
        [model] = addSinkReactions(model,{'pyr[c]','lac_L[s]'},[-1000 -1; 0 1000]);
        model.c(ismember(model.rxns,'sink_lac_L[s]'))=1;
        FBA = optimizeCbModel(model);
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> lac_L[s]';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;    
    %% pyr[c] -> lac_L[m]
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'lac_L[m]')))
		[model] = addSinkReactions(model,{'pyr[c]','lac_L[m]'},[-1000 -1; 0 1000]);
		model.c(ismember(model.rxns,'sink_lac_L[m]'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
	end	
    TestSolutionName{k,1} = 'pyr[c] -> lac_L[m]';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;
    %% pyr[c] -> lac_L[c]
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'lac_L[c]')))
		[model] = addSinkReactions(model,{'pyr[c]','lac_L[c]'},[-1000 -1; 0 1000]);
		model.c(ismember(model.rxns,'sink_lac_L[c]'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> lac_L[c]';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;
    %% trp_L -> melatn
    model = modelClosed;
    model.c(find(model.c)) = 0;
    if ~isempty(find(ismember(model.mets,'trp_L[c]'))) && ~isempty(find(ismember(model.mets,'melatn[c]')))
		[model] = addSinkReactions(model,{'trp_L[c]','melatn[c]'},[-1000 -1; 0 1000]);
		model.c(ismember(model.rxns,'sink_melatn[c]'))=1;
		FBA = optimizeCbModel(model);
		TestSolution(k,1) = FBA.f;
    else
		TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp_L -> melatn';
    TestedModel{k,1} = 'modelClosed';
    k = k +1;    
%%
VarNames = {'TestSolutionName', 'TestSolution'};
TestResults = table(TestSolutionName, TestSolution, 'VariableNames',VarNames);

end
