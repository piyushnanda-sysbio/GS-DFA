%Add path of Human-GEM and RAVEN functions 

addpath('Human-GEM\')
addpath('RAVEN\')

%Load the Human_GEM metabolic model in .mat format

load 'HumanGEM.mat';

%Change solver to gurobi

changeCobraSolver('gurobi','all');
setRavenSolver('gurobi');

%Load the essential task list and check whether the model is able to
%satisfy the tasks

essentialTasks=parseTaskList('metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

%Load gene expression/proteomics data

[val_UI_NHBES,val_I_NHBES]=expressionDataPrepare(1);
[val_UI_NHBEI,val_I_NHBEI]=expressionDataPrepare(2);
[val_UI_Biop,val_I_Biop]=expressionDataPrepare(3);

%Run tINIT algorithm to obtain models

[modelUI_NHBES,modelI_NHBES]=tINITModelIntegration(val_UI_NHBES,val_I_NHBES,ihuman);
[modelUI_NHBEI,modelI_NHBEI]=tINITModelIntegration(val_UI_NHBEI,val_I_NHBEI,ihuman);
[modelUI_Biop,modelI_Biop]=tINITModelIntegration(val_UI_Biop,val_I_Biop,ihuman);

modelI_NHBES=rmfield(modelI_NHBES,'rules')
modelI_NHBEI=rmfield(modelI_NHBEI,'rules')
modelI_Biop=rmfield(modelI_Biop,'rules')

writeCbModel(modelI_NHBES,'fileName','HumanGEMNHBESARS1.mat','format','mat')
writeCbModel(modelUI_NHBES,'fileName','HumanGEMNHBESARSMock1.mat','format','mat')
writeCbModel(modelI_NHBEI,'fileName','HumanGEMNHBEIAV1.mat','format','mat')
writeCbModel(modelUI_NHBEI,'fileName','HumanGEMNHBEIAVMock1.mat','format','mat')
writeCbModel(modelI_Biop,'fileName','HumanGEMBiopsySARS1.mat','format','mat')
writeCbModel(modelUI_Biop,'fileName','HumanGEMBiopsyMock1.mat','format','mat')
