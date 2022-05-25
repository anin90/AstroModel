function[model] = prepareModel(model)

    [model.rowlb, model.rowub] = deal(zeros(size(model.mets)));
    [model.int_vars] = zeros(size(model.rxns));

    genes = cellfun(@(x) strsplit(x, '.'), model.genes(:,1), 'UniformOutput', false);
    for i = 1:length(genes)
        genes{i} = genes{i,1}{1,1};
    end
    [model.genes_unique] = unique(genes);
    [model.genes] = sort(genes);
    [model.genes_unique_map] = findgroups(model.genes);
    
end