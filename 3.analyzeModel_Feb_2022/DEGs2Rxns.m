function[DEGs2Rxns_Table] = DEGs2Rxns(model, filename)
tic;

%% LOAD WT MODEL
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% DATA QUALITY CONTROL

% IMPORT DATA
DEG = importdata(filename);

% FILTER_1: EXCLUDE DUPLICATE GENES in DATA
A = DEG.Gene_Symbol;
[~,X,Z] = unique(A,'stable');
Y = histc(Z,1:numel(X))<2;
set = A(X(Y))';

[intis_1,ind] = ismember(A',(set));
filter_1.Entrez_ID = DEG.Entrez_ID(intis_1==1);
filter_1.Gene_Symbol = DEG.Gene_Symbol(intis_1==1);
filter_1.Status = DEG.Status(intis_1==1);
filter_1.HPA_Level = DEG.HPA_Level(intis_1==1);

% FILTER_2: EXCLUDE GENES WITHOUT ENTREZ_ID
[intis_inter] = ~(isnan(filter_1.Entrez_ID));
[filter_2.Entrez_ID] = filter_1.Entrez_ID(find(intis_inter));
[filter_2.Gene_Symbol] = filter_1.Gene_Symbol(find(intis_inter));
[filter_2.Status] = filter_1.Status(find(intis_inter));
[filter_2.HPA_Level] = filter_1.HPA_Level(find(intis_inter));

% FILTER_3: EXCLUDE GENES WITHOUT HPA EVIDENCE
HPA_Level_NA = strfind(filter_2.HPA_Level,'NA');
HPA_Level_YES = find(cellfun(@isempty,HPA_Level_NA));
gene = (filter_2.Entrez_ID(HPA_Level_YES)); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
[filter_3.Entrez_ID] = filter_2.Entrez_ID(HPA_Level_YES);
[filter_3.Gene_Symbol] = filter_2.Gene_Symbol(HPA_Level_YES);
[filter_3.Status] = filter_2.Status(HPA_Level_YES);
[filter_3.HPA_Level] = filter_2.HPA_Level(HPA_Level_YES);

% SPLIT GENES (UP/DOWN) (filter_1 / filter_2 / filter_3)
up = strfind(filter_2.Status,'up');
up_yes = find(~cellfun(@isempty,up));
genes_up = filter_2.Entrez_ID(up_yes); genes_up = num2cell(genes_up); genes_up = cellfun(@num2str,genes_up,'uni',0);
down = strfind(filter_2.Status,'down');
down_yes = find(~cellfun(@isempty,down));
genes_down = filter_2.Entrez_ID(down_yes); genes_down = num2cell(genes_down); genes_down = cellfun(@num2str,genes_down,'uni',0);
genes_up = strcat(genes_up,'.1');
genes_down = strcat(genes_down,'.1');

% FIND REACTIONS 'ACTIVE' WITH GENES
[list_up] = findRxnsActiveWithGenes(model, genes_up)'; %PostCobra
rxnList_up = unique(list_up); %PostCobra

[list_down] = findRxnsActiveWithGenes(model, genes_down)'; %PostCobra
rxnList_down = unique(list_down); %PostCobra

% REMOVE INTERSECTING REACTIONS
c = intersect(rxnList_up,rxnList_down);
rxnList_up(ismember(rxnList_up,c)) = [];
rxnList_down(ismember(rxnList_down,c)) = [];

%% PRINT RESULTS TO STRUCT
DEGs2Rxns_Table = struct('filter_1', {filter_1}, 'filter_2', {filter_2}, 'filter_3', {filter_3}, 'rxnList_down', {rxnList_down}, 'rxnList_up', {rxnList_up});

%%
toc;
end