function[expression] = prepareExpressionMatrix(data)

    data.matrix = quantilenorm(log10(data.matrix+1));
    HPA_Level_NA = strfind(data.HPA_Level,'NA');
    HPA_Level_YES = find(cellfun(@isempty,HPA_Level_NA)); %YES implies the protein expression can be low/med/high
    gene = (data.Entrez_ID(HPA_Level_YES)); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
    % expression.gene = strcat(gene,'.1');
    expression.gene = gene;
    expression.matrix = data.matrix(HPA_Level_YES,:);
    expression.maxvalue = max(data.matrix(HPA_Level_YES,:),[],2);
    
end