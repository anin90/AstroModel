function [TestResults] = Test4HumanAstro(model)

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
    TestSolutionGroup{k,1} = 'biomass';
    TestedMetabolite{k,1} = 'maintenance';
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
    TestSolutionGroup{k,1} = 'biomass';
    TestedMetabolite{k,1} = 'maintenance_noTrTr';
    TestedModel{k,1} = 'modelOri';
    k = k +1;

    %% PI5P4K - Inositol phosphate metabolism
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'PI5P4K'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'PI5P4K';
    TestSolutionGroup{k,1} = 'Inositol phosphate metabolism';
    TestedMetabolite{k,1} = 'NA';
    TestedModel{k,1} = 'modelOri';
    k = k +1;   

    %% MI1PP - Inositol phosphate metabolism
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'MI1PP'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'MI1PP';
    TestSolutionGroup{k,1} = 'Inositol phosphate metabolism';
    TestedMetabolite{k,1} = 'NA';
    TestedModel{k,1} = 'modelOri';
    k = k +1;   

    %% GLUDC - Glutamate metabolism
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'GLUDC'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'GLUDC';
    TestSolutionGroup{k,1} = 'Glutamate metabolism';
    TestedMetabolite{k,1} = 'NA';
    TestedModel{k,1} = 'modelOri';
    k = k +1;   
    
    %% ATP max, modelOri
%     model = modelOri;
%     model.c(find(model.c)) = 0;
%     model.c(ismember(model.rxns,'DM_atp_c_'))=1;
%     if find(model.c)>0
%         FBA = optimizeCbModel(model,'max','zero');
%         TestSolution(k,1) = FBA.f;
%     else
%         TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ATP max, modelOri';
%     TestSolutionGroup{k,1} = 'ATP demand';
%     TestedMetabolite{k,1} = 'All';
%     TestedModel{k,1} = 'modelOri';
%     k = k +1;

    %% ATP max aerobic, glc, v0.05
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.rxns,'EX_glc_D[e]')))
% 		model.lb(ismember(model.rxns,'EX_glc_D[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc_D[e]'))=-1;
% 		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
% 		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ATP max, aerobic, glc -> atp';
%     TestSolutionGroup{k,1} = 'ATP demand';
%     TestedMetabolite{k,1} = 'glc_D[e], aerobic';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;
    
    %% ATP max, anaerobic glc, v0.05
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.rxns,'EX_glc_D[e]')))
% 		model.lb(ismember(model.rxns,'EX_glc_D[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc_D[e]'))=-1;
% 		model.lb(ismember(model.rxns,'EX_o2[e]'))=0;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
% 		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ATP max, anaerobic, glc -> atp';
%     TestSolutionGroup{k,1} = 'ATP demand';
%     TestedMetabolite{k,1} = 'glc_D[e], anaerobic';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;   

  %% ATP max, aerobic, pyr -> atp
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.rxns,'EX_pyr[e]')))
% 		model.lb(ismember(model.rxns,'EX_pyr[e]'))=-1;model.ub(ismember(model.rxns,'EX_pyr[e]'))=-1;
% 		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
% 		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ATP max, aerobic, pyr -> atp';
%     TestSolutionGroup{k,1} = 'ATP demand';
%     TestedMetabolite{k,1} = 'pyr[e]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;
    
 %% ATP max, aerobic, lac_L -> atp
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.rxns,'EX_lac_L[e]')))
% 		model.lb(ismember(model.rxns,'EX_lac_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_lac_L[e]'))=-1;
% 		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
% 		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ATP max, aerobic, lac_L -> atp';
%     TestSolutionGroup{k,1} = 'ATP demand';
%     TestedMetabolite{k,1} = 'lac_L[e]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;       
    
  %% ATP max, aerobic, glu_L, synaptic -> atp
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.rxns,'EX_glu_L[s]')))
% 		model.lb(ismember(model.rxns,'EX_glu_L[s]'))=-1;model.ub(ismember(model.rxns,'EX_glu_L[s]'))=-1;
% 		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
% 		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ATP max, aerobic, glu_L (synaptic) -> atp';
%     TestSolutionGroup{k,1} = 'ATP demand';
%     TestedMetabolite{k,1} = 'glu_L[s]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;

 %% ATP max, aerobic, glu_L, extracellular -> atp
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.rxns,'EX_glu_L[e]')))
% 		model.lb(ismember(model.rxns,'EX_glu_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_glu_L[e]'))=-1;
% 		model.lb(ismember(model.rxns,'EX_o2[e]'))=-10;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
% 		model.c(ismember(model.rxns,'DM_atp_c_'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ATP max, aerobic, glu_L (extracellular) -> atp';
%     TestSolutionGroup{k,1} = 'ATP demand';
%     TestedMetabolite{k,1} = 'glu_L[e]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1; 
    
%% pyr[c] -> accoa[m]
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'accoa[m]')))
%         [model] = addSinkReactions(model,{'pyr[c]','accoa[m]'},[-1000 -1; 0 1000]);
%         model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
%         FBA = optimizeCbModel(model);
%         TestSolution(k,1) = FBA.f;
%     else
%         TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'pyr[c] -> accoa[m]';
%     TestSolutionGroup{k,1} = 'Pyruvate utilization';
%     TestedMetabolite{k,1} = 'accoa[m]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;    
    
    %% pyr[c] -> oaa[m]
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'oaa[m]')))
% 		[model] = addSinkReactions(model,{'pyr[c]','oaa[m]'},[-1000 -1; 0 1000]);
% 		model.c(ismember(model.rxns,'sink_oaa[m]'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
% 	else
% 		TestSolution(k,1) = NaN;
% 	end	
%     TestSolutionName{k,1} = 'pyr[c] -> oaa[m]';
%     TestSolutionGroup{k,1} = 'Pyruvate utilization';
%     TestedMetabolite{k,1} = 'oaa[m]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;    

    %% pyr[c] -> etoh[c]
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'etoh[c]')))
% 		[model] = addSinkReactions(model,{'pyr[c]','etoh[c]'},[-1000 -1; 0 1000]);
% 		model.c(ismember(model.rxns,'sink_etoh[c]'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
% 	end	
%     TestSolutionName{k,1} = 'pyr[c] -> etoh[c]';
%     TestSolutionGroup{k,1} = 'Pyruvate utilization';
%     TestedMetabolite{k,1} = 'etoh[c]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;
    
    %% pyr[c] -> lac_L[s]
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'lac_L[s]')))
%         [model] = addSinkReactions(model,{'pyr[c]','lac_L[s]'},[-1000 -1; 0 1000]);
%         model.c(ismember(model.rxns,'sink_lac_L[s]'))=1;
%         FBA = optimizeCbModel(model);
%         TestSolution(k,1) = FBA.f;
%     else
%         TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'pyr[c] -> lac_L[s]';
%     TestSolutionGroup{k,1} = 'Pyruvate utilization';
%     TestedMetabolite{k,1} = 'lac_L[s]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;    

    %% pyr[c] -> lac_L[m]
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'lac_L[m]')))
% 		[model] = addSinkReactions(model,{'pyr[c]','lac_L[m]'},[-1000 -1; 0 1000]);
% 		model.c(ismember(model.rxns,'sink_lac_L[m]'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
% 	end	
%     TestSolutionName{k,1} = 'pyr[c] -> lac_L[m]';
%     TestSolutionGroup{k,1} = 'Pyruvate utilization';
%     TestedMetabolite{k,1} = 'lac_L[m]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;

    %% pyr[c] -> lac_L[c]
%     model = modelClosed;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'pyr[c]'))) && ~isempty(find(ismember(model.mets,'lac_L[c]')))
% 		[model] = addSinkReactions(model,{'pyr[c]','lac_L[c]'},[-1000 -1; 0 1000]);
% 		model.c(ismember(model.rxns,'sink_lac_L[c]'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'pyr[c] -> lac_L[c]';
%     TestSolutionGroup{k,1} = 'Pyruvate utilization';
%     TestedMetabolite{k,1} = 'lac_L[c]';
%     TestedModel{k,1} = 'modelClosed';
%     k = k +1;
    
    %% ddcacoa[c] -> cit[m]
%     model = modelOri;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'ddcacoa[c]'))) && ~isempty(find(ismember(model.mets,'cit[m]')))
%         [model] = addSinkReactions(model,{'ddcacoa[c]','cit[m]'},[-1000 -1; 0 1000]);
% 		model.c(ismember(model.rxns,'sink_cit[m]'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ddcacoa[c] -> cit[m]';
%     TestSolutionGroup{k,1} = 'TCA citrate';
%     TestedMetabolite{k,1} = 'ddcacoa[c]';
%     TestedModel{k,1} = 'modelOri';
%     k = k +1;  

    %% ddcacoa[r] -> cit[m]
%     model = modelOri;
%     model.c(find(model.c)) = 0;
%     if ~isempty(find(ismember(model.mets,'ddcacoa[r]'))) && ~isempty(find(ismember(model.mets,'cit[m]')))
%         [model] = addSinkReactions(model,{'ddcacoa[r]','cit[m]'},[-1000 -1; 0 1000]);
% 		model.c(ismember(model.rxns,'sink_cit[m]'))=1;
% 		FBA = optimizeCbModel(model);
% 		TestSolution(k,1) = FBA.f;
%     else
% 		TestSolution(k,1) = NaN;
%     end
%     TestSolutionName{k,1} = 'ddcacoa[r] -> cit[m]';
%     TestSolutionGroup{k,1} = 'TCA citrate';
%     TestedMetabolite{k,1} = 'ddcacoa[r]';
%     TestedModel{k,1} = 'modelOri';
%     k = k +1;  

%%
VarNames = {'TestSolutionName', 'TestSolution', 'TestSolutionGroup', 'TestedMetabolite', 'TestedModel'};
TestResults = table(TestSolutionName, TestSolution, TestSolutionGroup, TestedMetabolite, TestedModel ,'VariableNames',VarNames);

end
