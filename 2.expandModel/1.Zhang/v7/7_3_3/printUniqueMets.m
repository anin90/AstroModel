tic;
modelMets = MBA_model_7_3_2_ASM.mets;
for str = 1:length(modelMets);
    str = modelMets{str};
    opens = str == '[';
    closes = str == ']';
    nestingcount = cumsum(opens - [0 closes(1:end-1)]);
    outstr = str(nestingcount == 0);
    disp(outstr);
end
clear modelMets
clear opens
clear closes
clear nestingcount
clear outstr
clear str
toc;