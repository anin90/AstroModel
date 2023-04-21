%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/PlotResults/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/PlotResults/plotDataOverlaysHyper_Tbl/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/PlotResults/plotDataOverlaysHyper_Tbl_Final/')

%% ST_Rxns

    % ST_S1
    Tbl = readtable('ST_S1.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/ST_S1.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % ST_S2
    Tbl = readtable('ST_S2.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/ST_S2.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    
%% DF_X_ST_Minerva

    % DF_1_ST_S1_Minerva
    Tbl = readtable('DF_1_ST_S1_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_1_ST_S1_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_2_ST_S2_Minerva
    Tbl = readtable('DF_2_ST_S2_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_2_ST_S2_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    
%% Final_ST_Minerva

    % ST_Rxns
    Tbl = readtable('ST_Rxns.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/ST_Rxns.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));    

%%
clearvars
