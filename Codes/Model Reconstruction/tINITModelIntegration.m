function [modelUI,modelI]=tINITModelIntegration(val_UI,val_I,ihuman)

refModelUI = ihuman;  % the reference model from which the GEM will be extracted
refModelI=ReactionAddO(ihuman); %The metabolic model with Virus VBOF added

%Additional purification model

refModelI=changeObjective(refModelI,{'Cov2VBOF'},1);

arrayData_UI.genes=refModelUI.genes;
arrayData_UI.tissues={'UI'}; %The name of the tissues involved
arrayData_UI.levels=val_UI; %Gene expression data
arrayData_UI.threshold=1; %The threshold ~ 75th percentile

arrayData_I.genes=refModelI.genes;
arrayData_I.tissues={'I'}; %The name of the tissues involved
arrayData_I.levels=val_I;%Gene expression data
arrayData_I.threshold=1; %The threshold ~ 75th percentile
%Define the parameters in getINITModel2 

essentialTasksUI=parseTaskList('metabolicTasks_Essential.xlsx');
essentialTasksI=parseTaskList('metabolicTasks_EssentialV.xlsx');

tissueUI = 'UI';% must match the tissue name in data_struct.tissues
tissueI='I';
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)

metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we alreawwwdy loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen

taskStructureUI = essentialTasksUI;  % metabolic task structure (used instead "taskFile")
taskStructureI = essentialTasksI;  % metabolic task structure (used instead "taskFile")

params.TimeLimit=5000;  % additional optimization parameters for the INIT algorithm. 5000 s is optimal. 
paramsFT = [];  % additional optimization parameters for the fit-tasks algorithm

%Run the getINITModel2 command. The execution takes around 2 hours.
modelUI = getINITModel2(refModelUI, tissueUI, celltype, hpaData, arrayData_UI, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructureUI, params, paramsFT);

modelI = getINITModel2(refModelI, tissueI, celltype, hpaData, arrayData_I, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructureI, params, paramsFT);

ui_bin=ismember(ihuman.rxns,modelUI.rxns);
i_bin=ismember(ihuman.rxns,modelI.rxns);
sum(~(ui_bin==i_bin));

%Check if the model can perform all the expected essential tasks. There
%shouldn't be any exchange reactions open. 
checkTasks(modelUI, [], true, false, false, essentialTasksUI);

modelI=removeMets(modelI,'sarscov2s');
modelI=ReactionAddO(modelI);
modelI=addBoundaryMets(modelI);
checkTasks(modelI, [], true, false, false, essentialTasksI);

%Simplify the model to add the exchange reaction and make it working. 
modelUI=simplifyModel(modelUI);
modelI=simplifyModel(modelI);

%The resulting model.b array is erroneus. This is to correct it. 
modelUI.b=repelem(0,length(modelUI.b))';
modelI.b=repelem(0,length(modelI.b))';

%Check the objective of the model using biomass objective (biomass)
checkObjective(modelUI);

modelI=changeObjective(modelI,'Cov2VBOF',1)
checkObjective(modelI);

%Run further analysis FVA, Flux sampling, Flux balance analysis

media=readtable('HumanGEMHAM.csv');
metList=table2array(media(:,8));
lb=table2array(media(:,9));
% model=simplifyModel(ihuman)
exc=modelI.rxns(findExcRxns(modelI));
mediacomps=exc(ismember(exc,findRxnsFromMets(modelI,metList)));

%Constrain media by HAM media

modelUI=changeRxnBounds(modelUI,exc,0,'l');
modelUI=changeRxnBounds(modelUI,mediacomps,lb,'l');

modelI=changeRxnBounds(modelI,exc,0,'l');
modelI=changeRxnBounds(modelI,mediacomps,-1,'l');

end