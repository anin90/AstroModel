function[exprs] = mapExpressionToModelGenes(model, expression)

    expression_tbl = table(expression.gene, expression.matrix);
    genes = cell2table(model.genes);
    modelGeneExp = outerjoin(genes,expression_tbl,'Type','left');
    modelGeneExp.Var1_expression_tbl=[];
    % modelGeneExp_Ctrl.Var2(isnan(modelGeneExp_Ctrl.Var2))=0;
    % modelGeneExp_Ctrl.Var3(isnan(modelGeneExp_Ctrl.Var3))=0;
    exprs = modelGeneExp.Var2;

end