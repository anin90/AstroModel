%% Load Data
FileName   = 'iAstro_Lithium.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% iAstro_Primary_TP
%% Akkouh_Li
dat = iAstro_Primary_TP_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_Primary_TP_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_Primary_TP_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_Primary_TP_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_Primary_TP_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_Primary_TP_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_Ctrl_TP
%% Akkouh_Li
dat = iAstro_iPS_Ctrl_TP_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_Ctrl_TP_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_Ctrl_TP_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_Ctrl_TP_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_Ctrl_TP_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_Ctrl_TP_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_BD_TP
%% Akkouh_Li
dat = iAstro_iPS_BD_TP_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_TP_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_BD_TP_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_TP_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_BD_TP_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_TP_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_BD_R_TP
%% Akkouh_Li
dat = iAstro_iPS_BD_R_TP_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_R_TP_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_BD_R_TP_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_R_TP_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_BD_R_TP_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_R_TP_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_BD_NR_TP
%% Akkouh_Li
dat = iAstro_iPS_BD_NR_TP_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_NR_TP_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_BD_NR_TP_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_NR_TP_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_BD_NR_TP_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_NR_TP_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T


%% iAstro_Primary_T
%% Akkouh_Li
dat = iAstro_Primary_T_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_Primary_T_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_Primary_T_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_Primary_T_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_Primary_T_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_Primary_T_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_Ctrl_T
%% Akkouh_Li
dat = iAstro_iPS_Ctrl_T_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_Ctrl_T_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_Ctrl_T_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_Ctrl_T_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_Ctrl_T_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_Ctrl_T_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_BD_T
%% Akkouh_Li
dat = iAstro_iPS_BD_T_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_T_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_BD_T_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_T_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_BD_T_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_T_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_BD_R_T
%% Akkouh_Li
dat = iAstro_iPS_BD_R_T_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_R_T_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_BD_R_T_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_R_T_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_BD_R_T_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_R_T_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%% iAstro_iPS_BD_NR_T
%% Akkouh_Li
dat = iAstro_iPS_BD_NR_T_Akkouh_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_NR_T_Akkouh_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE66276_Li
dat = iAstro_iPS_BD_NR_T_GSE66276_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_NR_T_GSE66276_Li.csv','WriteVariableNames',0); 
clear dat m T
%% GSE132397_Li
dat = iAstro_iPS_BD_NR_T_GSE132397_Li; 
m = vertcat(dat.rxnList_down, dat.rxnList_up_f, dat.rxnList_up_b);
T = table(m);
writetable(T,'RxnInput_iAstro_iPS_BD_NR_T_GSE132397_Li.csv','WriteVariableNames',0); 
clear dat m T

%%
clear all
