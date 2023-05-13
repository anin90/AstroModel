%%
tStart = tic;

%% Load Data
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to run 'Test4HumanAstro.m'
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/')

%% model_abs

    [TestResults_Primary] = Test4HumanAstro(iAstro_Primary_TP_abs);
    [TestResults_iPS_Ctrl_a] = Test4HumanAstro(iAstro_iPS_Ctrl_TP_abs);
    [TestResults_iPS_BD] = Test4HumanAstro(iAstro_iPS_BD_TP_abs);
    [TestResults_iPS_BD_R] = Test4HumanAstro(iAstro_iPS_BD_R_TP_abs);
    [TestResults_iPS_BD_NR] = Test4HumanAstro(iAstro_iPS_BD_NR_TP_abs);
    [TestResults_iPS_Ctrl_b] = Test4HumanAstro(iAstro_iPS_SCZ_Ctrl_TP_abs);
    [TestResults_iPS_HT] = Test4HumanAstro(iAstro_iPS_SCZ_HT_TP_abs);
    [TestResults_iPS_ST] = Test4HumanAstro(iAstro_iPS_SCZ_ST_TP_abs);
    VarNames_TP = {'TestSolutionName','Primary_Ctrl','iPS_Ctrl_a','iPS_BD','iPS_BD_R','iPS_BD_NR', 'iPS_Ctrl_b', 'iPS_HT', 'iPS_ST', 'TestSolutionGroup', 'TestedMetabolite','TestedModel'};
    Test4HumanAstro_abs = table(TestResults_Primary.TestSolutionName, TestResults_Primary.TestSolution, TestResults_iPS_Ctrl_a.TestSolution, TestResults_iPS_BD.TestSolution, TestResults_iPS_BD_R.TestSolution, TestResults_iPS_BD_NR.TestSolution, TestResults_iPS_Ctrl_b.TestSolution, TestResults_iPS_HT.TestSolution, TestResults_iPS_ST.TestSolution, TestResults_Primary.TestSolutionGroup, TestResults_Primary.TestedMetabolite, TestResults_Primary.TestedModel, 'VariableNames', VarNames_TP);
    writetable(Test4HumanAstro_abs,'TestResults_files/Test4HumanAstro_abs.csv','WriteVariableNames', true, 'Delimiter','\t');

%% model_norm_t1

    [TestResults_Primary] = Test4HumanAstro(iAstro_Primary_TP_norm_t1);
    [TestResults_iPS_Ctrl_a] = Test4HumanAstro(iAstro_iPS_Ctrl_TP_norm_t1);
    [TestResults_iPS_BD] = Test4HumanAstro(iAstro_iPS_BD_TP_norm_t1);
    [TestResults_iPS_BD_R] = Test4HumanAstro(iAstro_iPS_BD_R_TP_norm_t1);
    [TestResults_iPS_BD_NR] = Test4HumanAstro(iAstro_iPS_BD_NR_TP_norm_t1);
    [TestResults_iPS_Ctrl_b] = Test4HumanAstro(iAstro_iPS_SCZ_Ctrl_TP_norm_t1);
    [TestResults_iPS_HT] = Test4HumanAstro(iAstro_iPS_SCZ_HT_TP_norm_t1);
    [TestResults_iPS_ST] = Test4HumanAstro(iAstro_iPS_SCZ_ST_TP_norm_t1);
    VarNames_TP = {'TestSolutionName','Primary_Ctrl','iPS_Ctrl_a','iPS_BD','iPS_BD_R','iPS_BD_NR', 'iPS_Ctrl_b', 'iPS_HT', 'iPS_ST', 'TestSolutionGroup', 'TestedMetabolite','TestedModel'};
    Test4HumanAstro_norm_t1 = table(TestResults_Primary.TestSolutionName, TestResults_Primary.TestSolution, TestResults_iPS_Ctrl_a.TestSolution, TestResults_iPS_BD.TestSolution, TestResults_iPS_BD_R.TestSolution, TestResults_iPS_BD_NR.TestSolution, TestResults_iPS_Ctrl_b.TestSolution, TestResults_iPS_HT.TestSolution, TestResults_iPS_ST.TestSolution, TestResults_Primary.TestSolutionGroup, TestResults_Primary.TestedMetabolite, TestResults_Primary.TestedModel, 'VariableNames', VarNames_TP);
    writetable(Test4HumanAstro_norm_t1,'TestResults_files/Test4HumanAstro_norm_t1.csv','WriteVariableNames', true, 'Delimiter','\t');

%% model_norm_t2

    [TestResults_Primary] = Test4HumanAstro(iAstro_Primary_TP_norm_t2);
    [TestResults_iPS_Ctrl_a] = Test4HumanAstro(iAstro_iPS_Ctrl_TP_norm_t2);
    [TestResults_iPS_BD] = Test4HumanAstro(iAstro_iPS_BD_TP_norm_t2);
    [TestResults_iPS_BD_R] = Test4HumanAstro(iAstro_iPS_BD_R_TP_norm_t2);
    [TestResults_iPS_BD_NR] = Test4HumanAstro(iAstro_iPS_BD_NR_TP_norm_t2);
    [TestResults_iPS_Ctrl_b] = Test4HumanAstro(iAstro_iPS_SCZ_Ctrl_TP_norm_t2);
    [TestResults_iPS_HT] = Test4HumanAstro(iAstro_iPS_SCZ_HT_TP_norm_t2);
    [TestResults_iPS_ST] = Test4HumanAstro(iAstro_iPS_SCZ_ST_TP_norm_t2);    
    VarNames_TP = {'TestSolutionName','Primary_Ctrl','iPS_Ctrl_a','iPS_BD','iPS_BD_R','iPS_BD_NR', 'iPS_Ctrl_b', 'iPS_HT', 'iPS_ST', 'TestSolutionGroup', 'TestedMetabolite','TestedModel'};
    Test4HumanAstro_norm_t2 = table(TestResults_Primary.TestSolutionName, TestResults_Primary.TestSolution, TestResults_iPS_Ctrl_a.TestSolution, TestResults_iPS_BD.TestSolution, TestResults_iPS_BD_R.TestSolution, TestResults_iPS_BD_NR.TestSolution, TestResults_iPS_Ctrl_b.TestSolution, TestResults_iPS_HT.TestSolution, TestResults_iPS_ST.TestSolution, TestResults_Primary.TestSolutionGroup, TestResults_Primary.TestedMetabolite, TestResults_Primary.TestedModel, 'VariableNames', VarNames_TP);
    writetable(Test4HumanAstro_norm_t2,'TestResults_files/Test4HumanAstro_norm_t2.csv','WriteVariableNames', true, 'Delimiter','\t');
    
%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));    

%%
clearvars
