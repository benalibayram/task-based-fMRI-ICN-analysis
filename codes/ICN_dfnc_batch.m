%% Run dfnc, modified version for HUBAL
clear

icatb_defaults; % GIFT toolbox function
global DETRENDNUMBER;

%% Initialise
prefix = 'test_dfnc'; 
TR = 2; % TR in seconds
outputDir = pwd;

% Component Time Courses (tc)
tc_file = '../data/ICA_timeseries_loaded_C15.mat';
cond_names = {'pview', 'smotor', 'srtt', 'gonogo', 'oneback', 'twoback', 'threeback'};
labels = {'VN_Medial' 'VN_Lateral' 'SMN_Lateral' 'DMN_PCC'...
    'FPN_Right' 'BG' 'CN' 'FPN_Left' 'DMN_MPFC'...
    'SMN_Superior' 'DAN' 'LN' 'SMN_Left'...
    'DMN_IPL' 'VAN'};
nRSN = length(labels);

% Preprocess params
preprocess_params.detrend_no = DETRENDNUMBER; % default 3
preprocess_params.doDespike = 'yes';
preprocess_params.tc_filter = 0.15; % Cutoff in Hz

% Windowing params
windowing_params.method = 'correlation';
windowing_params.wsize = 22; % 2*22 = 44 s
windowing_params.window_alpha = 3;
windowing_params.window_type = 'gaussian';
windowing_params.num_L1_repetitions = 10;

covariates = 'none';
covariateInfo = [];

% model time course
spm_file = '../data/GLM_design_mdes.mat';
selectedRegressors = {'pview', 1;'smotor', 2;'srtt', 3;'gonogo', 3;'1geri', 4;'2geri', 4;'3geri', 4};
task_onset  = [1   141 281 568];
task_offset = [140 280 567 995];

%% ----------------Time course load-------------------

load(tc_file, 'tcourses')
tc = tcourses;


%% ----------------modelX load---------------------
load(spm_file, 'SPM'); % ../data/ICA_timeseries_loaded_C15.mat
XX = SPM.xX.X;
allRegressors = {SPM.Sess.Fc.name};
indRegressors = {SPM.Sess.Fc.i};
numOfReg = size(selectedRegressors,1); % 7

%% ----------------dynamic connectivity calc.---------------------
numOfSub = 28;
numOfSess = 4;
nT = [140 140 287 428]; % number of time points for each session.
numComp = size(tc, 3); % 15 components
minTR = TR; % 2 sn

disp('----------------------------------------------------');
disp('-------------- STARTING DYNAMIC FNC --------------');
disp('----------------------------------------------------');

TC = cell(numOfSub, numOfSess); % 28x4
FNCdyn = cell(numOfSub, numOfSess); % 28x4
task_connectivity_R2 = cell(numOfSub, numOfSess); % 28x4
task_connectivity = cell(numOfSub, numOfSess); % 28x4
model_tcwin = cell(numOfSub, numOfSess);  % 28x4
XWin = cell(numOfSub, numOfSess);  % 28x4

selectedRegressors2 = cell2mat(selectedRegressors(:,2)); % [1 2 3 3 4 4 4]
task_connectivity_all = cell(0);

for iTask = 1:numOfSess % 1:4
    % Session time point indexes
    taskTp = task_onset(iTask):task_offset(iTask); % 1:140, 141:280, 281:567, 568:995
    
    % Regressors correspond this session
    iselectedRs = find(selectedRegressors2 == iTask); % 1, 2, [3 4], [5 6 7]
    if (isempty(iselectedRs))
        continue % skip, if no regressor in this session
    end
    
    % Regressor models
    modelX = [];
    for iselectedR = iselectedRs' % 1, 2, [3 4]', [5 6 7]'
        iR = strcmp(selectedRegressors{iselectedR,1}, allRegressors); % ['pview'], ['smotor'], ['srtt', 'gonogo'], ['1geri', '2geri', '3geri']
        model_iR = indRegressors(iR); 
        modelX = [modelX, XX(taskTp, model_iR{:}(1))];% without model derivatives, 140x1, 140x1, 287x2, 428x3
    end
    clear iR model_iR 
    
    minTP = nT(iTask); %[140 140 287 428]
    for iSub = 1:numOfSub
        current_tc = squeeze(tc(iSub, taskTp, :)); % 140x15,  140x15, 287x15, 428x15
        corrInfo = icatb_compute_task_dfnc_HUBAL(current_tc, TR, 'mintr', minTR, 'minTP', minTP, 'preprocess_params', preprocess_params, ...
        'windowing_params', windowing_params, 'modelTC', modelX, 'covariateInfo', covariateInfo);
        
        disp(['.... current Task:' 9 num2str(iTask) 9 'Subject:' 9 num2str(iSub)]);
        TC(iSub, iTask) = {corrInfo.X}; 
        FNCdyn(iSub, iTask) = {corrInfo.FNCdyn}; % 28x4 cell
        % 118x105, 118x105, 265x105, 406x105 time points wsize(22) shortened
        % 105 = 15x14/2
        task_connectivity_R2(iSub, iTask) = {corrInfo.task_connectivity_R2};
        task_connectivity(iSub, iTask) = {corrInfo.task_connectivity}; % 28x4 cell
        % 1x105, 1x105, 2x105, 3x105
        model_tcwin(iSub, iTask) = {corrInfo.model_tcwin}; % 28x4 cell 
        % 118x1, 118x1, 265x2, 406x3
        XWin(iSub, iTask) = {corrInfo.XWin};
    end
    
    task_connectivity_current = cell2mat(task_connectivity(:, iTask)); % (nR x 28) x 105
    % Separating: 'pview', 'smotor', 'srtt', 'gonogo', '1geri', '2geri', '3geri'
    nR = length(iselectedRs); % 1, 1, 2, 3
    if nR > 1
        for j = 1:nR % 1:2, 1:3
            task_connectivity_beta_all(end+1) = {task_connectivity_current(j:nR:end, :)};
        end
    else
        task_connectivity_beta_all(iTask) = {task_connectivity_current};
    end
end

%% executive Regressors {'pview', 1;'smotor', 2;'srtt', 3;'gonogo', 3}
exeReg = 1:4;
exeNames = cond_names(exeReg)';

perm = nchoosek(1:nRSN, 2);

couplenames_exe = join([repelem(exeNames, length(perm), 1),...
                        repmat(labels(perm), length(exeReg), 1)],'_')';

exe_task_connectivity_beta = cell2mat(task_connectivity_beta_all(exeReg));

exe_task_table_beta = cell2table(num2cell(exe_task_connectivity_beta), 'VariableNames',couplenames_exe);

writetable(exe_task_table_beta,'Beta_task_connectivity_Executive.xlsx');
disp('Beta_task_connectivity_Executive.xlsx')


%% memoload Regressors {'srtt', 3;'1geri', 4;'2geri', 4;'3geri', 4}
wmReg = [3 5 6 7];
wmNames = cond_names(wmReg)';

% perm

couplenames_wm = join([repelem(wmNames, length(perm), 1),...
                        repmat(labels(perm), length(wmReg), 1)],'_')';
 
wm_task_connectivity_beta = cell2mat(task_connectivity_beta_all(wmReg));

wm_task_table_beta = cell2table(num2cell(wm_task_connectivity_beta), 'VariableNames',couplenames_wm);

writetable(wm_task_table_beta,'Beta_task_connectivity_Memoload.xlsx');
disp('Beta_task_connectivity_Memoload.xlsx')



%%
fprintf('\n');
disp('Analysis Complete');
fprintf('\n');

