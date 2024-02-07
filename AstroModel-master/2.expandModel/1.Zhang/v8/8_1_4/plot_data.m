tic;
%% literature rxns
model = iMAT_model_TP_EXP_ASM_BBB;
data = T12.FSr_merged(T12.FSr_L_index,:);
imagesc(data); colorbar;
xtick = {'MI1PP','INSTt4','PGS','PDE1','BPNT','PGMT','ADNCYC','MI1PP+INSTt4','All'};
ytick = model.rxns(T12.FSr_L_index);
xlabel('Lithium models'); ylabel('Reactions with FSr>1.5');
ax = gca; ax.YTickLabel = ytick; ax.YTickMode = 'manual'; ax.YTick = 1:1:length(ytick);
ax.XTickLabel = xtick; ax.XTickLabelRotation = 45;
hold on;
n=size(data,1);
for i = 1:n
    plot([.5,n+.5],[i-.5,i-.5],'k-');
    plot([i-.5,i-.5],[.5,n+.5],'k-');
end
clear model data ax xtick ytick n i

%% subsystems (literature + omics)
data = importdata('plot/subsystem2mat.csv');
imagesc(data.data); colorbar;
xtick = {'Literature','Akkouh-Li','Akkouh-Li+HPA','GSE66276-Li','GSE66276-Li+HPA','GSE132397-Li', 'GSE132397-Li+HPA'};
ytick = data.textdata;
%xlabel('Lithium models'); ylabel('subsystems with FSr>1.5');
ax = gca; ax.YTickLabel = ytick; ax.YTickMode = 'manual'; ax.YTick = 1:1:length(ytick);
ax.XTickLabel = xtick; ax.XTickLabelRotation = 45;
ax.XTickMode = 'manual'; ax.XTick = 1:1:length(xtick);
hold on;
n=size(data.data,1);
for i = 1:n
plot([.5,n+.5],[i-.5,i-.5],'k-');
plot([i-.5,i-.5],[.5,n+.5],'k-');
end
clear data xtick ytick ax ans n i

%%
tic