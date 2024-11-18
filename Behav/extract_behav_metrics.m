clear all
clc

% Set path
main_dir = 'C:\Users\Utente\Desktop\Work\MR-Methods\EIB\';
input_dir = fullfile(main_dir,'fMRS_fMRI_NBACK_project');
% input_dir = '\\mri-storage.cimec.unitn.it\mumi-ext-001\mumi-data\USERS\Asia_Ferrari\fMRS_fMRI_NBACK_project';

% Load demographics info 
info = readtable(fullfile(main_dir,'info_demographics_35s.csv'));

% Specify IDs of partecipants to be included in the analysis
ID_fMRI = [1:3,5:36];
ID_fMRS = [12,13,14,15,19,20,23,24,25,28,34,36];
ID_list = {ID_fMRI,ID_fMRS};

% Loop through acquisiton (fMRI -> acq_iter = 1; fMRS -> acq_iter = 2)
for acq_iter = 1:2
    sub_id_list = ID_list{acq_iter};
    all_rt_wide_0back = [];
    all_rt_wide_1back = [];
    all_rt_wide_2back = [];

    % Loop through participants
    for sub_iter = 1:length(sub_id_list)
        if sub_id_list(sub_iter) < 10
            sub_id = strcat("sub0",string(sub_id_list(sub_iter)));
        else
            sub_id = strcat("sub",string(sub_id_list(sub_iter)));
        end

        if acq_iter == 1
            input_0back_file_name = dir(fullfile(input_dir,sub_id,"sub*001_task-0nback*xlsx"));
            input_1back_file_name = dir(fullfile(input_dir,sub_id,"sub*003_task-1nback*xlsx"));
            input_2back_file_name = dir(fullfile(input_dir,sub_id,"sub*005_task-2nback*xlsx"));
        else
            input_0back_file_name = dir(fullfile(input_dir,sub_id,"sub*002_task-0nback*xlsx"));
            input_1back_file_name = dir(fullfile(input_dir,sub_id,"sub*004_task-1nback*xlsx"));
            input_2back_file_name = dir(fullfile(input_dir,sub_id,"sub*006_task-2nback*xlsx"));
        end
        input_0back_file = readtable(fullfile(input_dir,sub_id,input_0back_file_name.name));
        input_1back_file = readtable(fullfile(input_dir,sub_id,input_1back_file_name.name));
        input_2back_file = readtable(fullfile(input_dir,sub_id,input_2back_file_name.name));

        % Loop through n-back task
        input_file_list = {input_0back_file,input_1back_file,input_2back_file};
        for task_iter = 1:length(input_file_list)
            input_file = input_file_list{task_iter};

            % Loop through block
            n_blocks = 4;
            block_cols = 1:5:20;
            hit_cnt = 0;
            miss_cnt = 0;
            fa_cnt = 0;
            cr_cnt = 0;
            rt_by_block = nan(9,4);

            for block_iter = 1:n_blocks

                % Extract number of hits
                hit_idx = find(input_file{:,block_cols(block_iter)} == 1 & input_file{:,block_cols(block_iter)+1} == 1);
                hit_cnt = hit_cnt + length(hit_idx);

                % Extract number of false alarms
                fa_idx = find(input_file{:,block_cols(block_iter)} == 1 & input_file{:,block_cols(block_iter)+1} == 0);
                fa_cnt = fa_cnt + length(fa_idx);

                % Extract reaction times only for congruent trials
                rt_congruent = input_file{hit_idx,block_cols(block_iter)+2};
                rt_by_block(1:length(rt_congruent),block_iter) = rt_congruent; 

            end

            % Compute and store accuracy across blocks
            n_congr = 36; % number of congruent trials in each task by design
            acc = (hit_cnt*100)/n_congr; 
            all_acc{sub_iter,task_iter} = acc;

            % Compute and store mean reaction times across blocks
            rt_mean = mean(rt_by_block,"all", "omitnan");
            all_rt_mean{sub_iter,task_iter} = rt_mean;

            % Compute d prime:
            % 1) Adjust perfect scores by employing a log-linear approach
            % (https://stats.stackexchange.com/questions/134779/d-prime-with-100-hit-rate-probability-and-0-false-alarm-probability)
            signal_trials_ratio = n_congr/(size(input_file,1)*n_blocks); % 44*4 = number of total trials in each task by design  
            hit_cnt_adjusted = hit_cnt + signal_trials_ratio;
            n_signal_trials_adjusted = n_congr + 2*signal_trials_ratio;

            noise_trials_ratio = 1 - signal_trials_ratio;
            fa_cnt_adjusted = fa_cnt + (noise_trials_ratio);
            n_noise_trials = size(input_file,1)*n_blocks - n_congr;
            n_noise_trials_adjusted = n_noise_trials + 2*noise_trials_ratio;

            % 2) Compute hit rate and false alarm rate
            hit_rate = hit_cnt_adjusted/n_signal_trials_adjusted;
            fa_rate = fa_cnt_adjusted/n_noise_trials_adjusted;

            % 3) Then z trasform hit and false alarm rates by computing the inverse of the standard normal cumulative distribution
            % (https://doi.org/10.1080/13803391003596421)
            zhr = norminv(hit_rate,0,1);
            zfar = norminv(fa_rate,0,1);

            % 4) Finalize (https://doi.org/10.1080/13803391003596421) and store d prime computation
            dprime = zhr - zfar;
            all_d{sub_iter,task_iter} = dprime;
            
            % Prepare output structures for saving reaction times for each congruent trial  
            rt_wide_values = reshape(rt_by_block,1,[]);
            rt_wide_values = rt_wide_values(~isnan(rt_wide_values));
            n_risp = length(rt_wide_values);
            n_block = ~isnan(rt_by_block) .* [1,2,3,4];
            n_block = reshape(n_block,1,[]);
            n_block = n_block(n_block > 0);
            info_row = find(strcmp(sub_id,info{:,1}));
            rt_wide = [repmat(sub_id,n_risp,1), repmat(info{info_row,2},n_risp,1), ...
                repmat(info{info_row,3},n_risp,1), n_block', repmat(n_risp,n_risp,1), ...
                rt_wide_values'];
            if task_iter == 1
                rt_wide = [rt_wide repmat("0back",n_risp,1)];
                all_rt_wide_0back = [all_rt_wide_0back; rt_wide];
            elseif task_iter == 2
                rt_wide = [rt_wide repmat("1back",n_risp,1)];
                all_rt_wide_1back = [all_rt_wide_1back; rt_wide];
            else
                rt_wide = [rt_wide repmat("2back",n_risp,1)];
                all_rt_wide_2back = [all_rt_wide_2back; rt_wide];
            end
        end
    end

    % Create output structure
    if acq_iter == 1 
        output_table = info;
    else
        output_table = info([11,12,13,14,18,19,22,23,24,27,33,35],:);
    end
    output_table.Acc_0back = all_acc(:,1);
    output_table.Acc_1back = all_acc(:,2);
    output_table.Acc_2back = all_acc(:,3);
    output_table.dprime_0back = all_d(:,1);
    output_table.dprime_1back = all_d(:,2);
    output_table.dprime_2back = all_d(:,3);
    output_table.RT_0back = all_rt_mean(:,1);
    output_table.RT_1back = all_rt_mean(:,2);
    output_table.RT_2back = all_rt_mean(:,3);

    header = ["Subject", "Gender", "Age", "Block", "Nrisp", "RT", "Runs"];
    all_rt_wide_0back = [header; all_rt_wide_0back];
    all_rt_wide_1back = [header; all_rt_wide_1back];
    all_rt_wide_2back = [header; all_rt_wide_2back];

    % Export output structures
    if acq_iter == 1
        writetable(output_table, fullfile(main_dir,'EIB_behav_acc_dprime_meanRT_fMRI_35s.csv'));
        writematrix(all_rt_wide_0back, fullfile(main_dir,'EIB_behav_RT_wide_0back_fMRI_35s.csv'));
        writematrix(all_rt_wide_1back, fullfile(main_dir,'EIB_behav_RT_wide_1back_fMRI_35s.csv'));
        writematrix(all_rt_wide_2back, fullfile(main_dir,'EIB_behav_RT_wide_2back_fMRI_35s.csv'));
    else
        writetable(output_table, fullfile(main_dir,'EIB_behav_acc_dprime_meanRT_fMRS_12s.csv'));
        writematrix(all_rt_wide_0back, fullfile(main_dir,'EIB_behav_RT_wide_0back_fMRS_12s.csv'));
        writematrix(all_rt_wide_1back, fullfile(main_dir,'EIB_behav_RT_wide_1back_fMRS_12s.csv'));
        writematrix(all_rt_wide_2back, fullfile(main_dir,'EIB_behav_RT_wide_2back_fMRS_12s.csv'));
    end

    clear("all_acc","all_d","all_rt_mean")
end
