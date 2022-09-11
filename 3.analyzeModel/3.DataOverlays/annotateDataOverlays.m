%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl_Final/')

%% BD_Lumped

    % BD_S1
    Tbl = readtable('BD_S1.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S1.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_S2
    Tbl = readtable('BD_S2.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S2.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_S3
    Tbl = readtable('BD_S3.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S3.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% BD_R

    % BD_R_S1
    Tbl = readtable('BD_R_S1.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S1.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_R_S2
    Tbl = readtable('BD_R_S2.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S2.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_R_S3
    Tbl = readtable('BD_R_S3.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S3.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% BD_NR

    % BD_NR_S1
    Tbl = readtable('BD_NR_S1.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S1.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S2
    Tbl = readtable('BD_NR_S2.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S2.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S3
    Tbl = readtable('BD_NR_S3.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S3.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S4
    Tbl = readtable('BD_NR_S4.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula'};
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S4.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));    

%%
clearvars
