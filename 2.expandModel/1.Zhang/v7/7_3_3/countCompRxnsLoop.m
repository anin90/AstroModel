tic;
model = MBA_model_7_3_2;
model = changeObjective(model, 'biomass_maintenance', 0);

comp = {'[e]','[s]','[c]','[g]','[l]','[m]','[n]','[r]','[x]','[i]'};
for k = 1:length(comp);
    disp(length(findRxnFromCompartment(model, comp{k})));
    disp(comp{k});
end
toc
%[s] = findRxnFromCompartment(model,'[s]');
%[c] = findRxnFromCompartment(model,'[c]');
%[g] = findRxnFromCompartment(model,'[g]');
%[l] = findRxnFromCompartment(model,'[l]');
%[m] = findRxnFromCompartment(model,'[m]');
%[n] = findRxnFromCompartment(model,'[n]');
%[r] = findRxnFromCompartment(model,'[r]');
%[x] = findRxnFromCompartment(model,'[x]');
%[i] = findRxnFromCompartment(model,'[i]');

clear model
clear comp
clear k
