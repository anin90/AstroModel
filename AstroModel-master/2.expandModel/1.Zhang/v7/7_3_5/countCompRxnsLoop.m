tic;
model = MBA_model_7_3_4_ASM;
model = changeObjective(model, 'biomass_maintenance', 0);

% comp = {'[e]','[s]','[c]','[g]','[l]','[m]','[n]','[r]','[x]','[i]'};
% for k = 1:length(comp);
%     disp(length(findRxnFromCompartment(model, comp{k})));
%     disp(comp{k});
% end

[comp.s] = findRxnFromCompartment(model,'[s]');
[comp.c] = findRxnFromCompartment(model,'[c]');
[comp.g] = findRxnFromCompartment(model,'[g]');
[comp.l] = findRxnFromCompartment(model,'[l]');
[comp.m] = findRxnFromCompartment(model,'[m]');
[comp.n] = findRxnFromCompartment(model,'[n]');
[comp.r] = findRxnFromCompartment(model,'[r]');
[comp.x] = findRxnFromCompartment(model,'[x]');
[comp.i] = findRxnFromCompartment(model,'[i]');

clearvars -except MBA_model_7_3 MBA_model_7_3_4_ASM modelStats comp
toc;
