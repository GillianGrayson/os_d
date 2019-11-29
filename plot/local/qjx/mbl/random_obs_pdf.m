clear all;

num_bins_obs = 100;

obs_start_id = 10;
obs_num = 10;

sys_id = 3;
task_id = 4;
prop_id = 0;

seed = 1;
mns = 1000000;

Nc = 8;

W_seed = 10;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

prm_W = 20.0;
prm_U = 1.0;
prm_J = 1.0;

start_type = 0;
start_state = 49;

data_path = '../../../../source/cpp/QJX/QJX';


for tr_id = 0 : 9
    
    fig = figure;
    for obs_id = obs_start_id : obs_start_id + obs_num - 1
        
        suffix = sprintf('%d_%d_setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
            tr_id, ...
            obs_id, ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            seed, ...
            mns, ...
            Nc, ...
            W_seed, ...
            W_mns, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            prm_W, ...
            prm_U, ...
            prm_J, ...
            start_type, ...
            start_state);
        
        fn = sprintf('%s/random_obs_evo_%s.txt', data_path, suffix);
        data = importdata(fn);
        obs = zeros(size(data, 1), 1);
        for id = 1:size(data, 1)
            str = string(data(id));
            curr_data = sscanf(str, '(%e,%e)', 2);
            obs(id) = curr_data(1);
        end
        
        max_obs = max(obs) + 1e-8;
        min_obs = min(obs) - 1e-8;
        
        shift_obs = (max_obs - min_obs) / num_bins_obs;
        bins_obs = linspace(min_obs + 0.5 * shift_obs, max_obs - 0.5 * shift_obs, num_bins_obs)';
        pdf_obs = zeros(num_bins_obs, 1);
        
        for id = 1:size(data, 1)
            curr_obs = obs(id);
            bin_id = floor((curr_obs - min_obs) / (max_obs - min_obs + eps) * num_bins_obs) + 1;
            pdf_obs(bin_id) = pdf_obs(bin_id) + 1;
        end
        
        sum_pdf = sum(pdf_obs);
        pdf_obs = pdf_obs / (sum_pdf * shift_obs);
        norm = sum(pdf_obs) * shift_obs;
        norm_diff = 1.0 - norm
        
        h = plot(bins_obs, pdf_obs, 'LineWidth', 3);
        legend(h, sprintf('obs=%d', obs_id))
        set(gca, 'FontSize', 30);
        xlabel('random observable', 'Interpreter', 'latex');
        set(gca, 'FontSize', 30);
        ylabel('$PDF$', 'Interpreter', 'latex');
        hold all;
        
        title_str = sprintf('$Nc=%d$ \\quad $W=%d$ \\quad $U=%d$ \\quad $J=%d$ \\quad $tr=%d$', ...
            Nc, ...
            prm_W, ...
            prm_U, ...
            prm_J, ...
            tr_id);
        title(title_str, 'interpreter', 'latex', 'FontSize', 20)
        
    end
    
    propertyeditor(fig)
    
    
end



