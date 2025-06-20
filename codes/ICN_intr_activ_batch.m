clear all

marsbar('on') % SPM12 toolbox marsbar patched with 'mars_arm_call.m'
spm('defaults', 'fmri');

% **********************Design creation: D *************************
load('../data/GLM_design_mdes.mat') % stores SPM.mat
nder = 2; % 2:1derivatives
% 1deriv = ( 9 x 2 ) + 3 = 21
cont_default = 21;

D = des_struct(mardo_5, SPM);

SPM.xsDes = struct(...
	'Basis_functions',	'hrf (with time derivative)',...
	'Number_of_sessions',	1,...
	'Trials_per_session',	length(SPM.Sess.U),...
	'Global_calculation',	'none',...
	'Grand_mean_scaling',	'none',...
	'Global_normalisation',	'none',...
    'Interscan_interval',   sprintf('%0.2f {s}',SPM.xY.RT));
K = struct(	'HParam',	512,...
			'row',		SPM.Sess.row,...
			'RT',		SPM.xY.RT);

k = length(K.row);
n = fix(2*(k*K.RT)/K.HParam + 1);
X0 = spm_dctmtx(k,n);
K.X0 = X0(:,2:end);

SPM.xX.K = K;
SPM.xsDes.High_pass_Filter = sprintf('Cutoff: %d {s}', 512);

SPM.xVi.V  = speye(995);
f2cl       = 'Vi';
SPM.xVi.form = 'i.i.d';

if isfield(SPM.xVi, f2cl)
    SPM.xVi = rmfield(SPM.xVi, f2cl);
end

xsDes = struct('Serial_correlations', SPM.xVi.form);
SPM.xsDes = mars_struct('ffillmerge', SPM.xsDes, xsDes);

D = des_struct(D,SPM);
mars_arm('update', 'def_design', D);

%%
% ICA_timeseries_load, 15 Components representing ICNs, 28 subjects
load('../data/ICA_timeseries_loaded_C15.mat')
clearvars -except tcourses D nder cont_default

labels = {'VN_Medial' 'VN_Lateral' 'SMN_Lateral' 'DMN_PCC'...
    'FPN_Right' 'BG' 'CN' 'FPN_Left' 'DMN_MPFC'...
    'SMN_Superior' 'DAN' 'LN' 'SMN_Left'...
    'DMN_IPL' 'VAN'};

nsubj = size(tcourses,1);

cont_names = {'pview', 'smotor', 'srtt', 'gonogo', 'oneback', 'twoback', 'threeback'};
cont_types = repmat({'T'}, 1, length(cont_names));
cont_vals = {[1 zeros(1,cont_default)],...
            [zeros(1, 1*nder) 1 zeros(1, cont_default-1*nder)],...
            [zeros(1, 2*nder) 1 zeros(1, cont_default-2*nder)],...
            [zeros(1, 3*nder) 1 zeros(1, cont_default-3*nder)],...
            [zeros(1, 4*nder) 1 zeros(1, cont_default-4*nder)],...
            [zeros(1, 5*nder) 1 zeros(1, cont_default-5*nder)],...
            [zeros(1, 6*nder) 1 zeros(1, cont_default-6*nder)]};
tcourses(29,:,:) = squeeze(mean(tcourses(:,:,:),1)); %subject mean saved as 29th subj, for visualisation purposes.
for subj_ind = 1:nsubj+1
    Y = squeeze(tcourses(subj_ind, :,:));
    for rno = 1:length(labels)
        r_data{rno} = Y(:,rno);
        r_info{rno} = struct(...
            'name', labels{rno},...
            'descrip', '',...
            'Y', Y(:,rno),...
            'weights', [],...
            'info', [],...
            'vXYZ', [],...
            'mat', []);
    end
    s_info = struct(...
        'sumfunc', 'unknown', ...
        'descrip', 'Data from matlab',...
        'info', []);
    marsY = marsy(r_data, r_info, s_info);
    
    E = estimate(D, marsY);
    
    [E, Ic] = add_contrasts(E, cont_names, cont_types, cont_vals);
    
    stat_struct(subj_ind) = compute_contrasts(E, Ic);
end


%%
n=1;
for cont_ex = 1:length(cont_names)
    for ICN_ex = 1:length(labels)
        names_exe(n) = {sprintf('%s_%s',cont_names{cont_ex},labels{ICN_ex})};
        n=n+1;
    end
end

%%

for i=1:28
    dummy = stat_struct(i).stat';
    f512_none_v5(i,:) = dummy(:)';
end

T_f512_none_v5 = cell2table(num2cell(f512_none_v5),'VariableNames',names_exe);

writetable(T_f512_none_v5,'C15_ICN_intr_activ_output.xls');

