clear all;

path_figures = 'C:/YandexDisk/Work/os_d/drafts/floquet/ospm/lambdas';

path_data = 'C:/Work/os_d/source/cpp/QJX/QJX';

cmap_type = 'rwb'; %% 'hot'

sys_id = 9;
task_id = 7;
prop_id = 1;
seed = 1;
mns = 1000000;

num_trajectories = 50;
num_tp_periods = 100;
num_obs_periods = 100;

nph = 1;
Delta = 0.10;
Omega = 1.0;
gamma = 0.01;

lpn_type = 0;
lpn_delta_s = log10(5.0e-1);
lpn_delta_f_high = log10(5.0e-1);
lpn_delta_f_low = log10(5.0e-1);

ampl_start = 0.6;
ampl_shift = 0.6;
ampl_num = 50;
ampls = linspace(ampl_start, ampl_start + (ampl_num - 1) * ampl_shift, ampl_num);

freq_start = 0.07;
freq_shift = 0.07;
freq_num = 50;
freqs = linspace(freq_start, freq_start + (freq_num - 1) * freq_shift, freq_num);

lambdas = zeros(ampl_num, freq_num);

for freq_id = 1:freq_num
    freq = freq_start + (freq_id - 1) * freq_shift
    for ampl_id = 1:ampl_num
        ampl = ampl_start + (ampl_id - 1) * ampl_shift;

        suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_nph(%d)_ampl(%0.4f)_freq(%0.4f)_D(%0.4f)_O(%0.4f)_gamma(%0.4f)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            seed, ...
            mns, ...
            nph, ...
            ampl, ...
            freq, ...
            Delta, ...
            Omega, ...
            gamma, ...
            lpn_type, ...
            lpn_delta_s, ...
            lpn_delta_f_high, ...
            lpn_delta_f_low);

        path = sprintf('%s/lambda_%s.txt', path_data, suffix);
        data = importdata(path);

        lambdas(ampl_id, freq_id) = lambdas(ampl_id, freq_id) + mean(data(num_trajectories / 2 + 1:end));
    end  
end

fig = figure;
hLine = imagesc(ampls, freqs, lambdas');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\omega$', 'Interpreter', 'latex');
if strcmp(cmap_type,'rwb') == 1
    h = colorbarpwn(min(lambdas(:)), max(lambdas(:)), 'level', 512);
else
    colormap hot;
    h = colorbar;
end
set(gca, 'FontSize', 30);
title(h, '$\lambda$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;

suffix = sprintf('traj(%d)_tp(%d)_obs(%d)_n(%d)_ampl(%0.4f_%0.4f_%d)_freq(%0.4f_%0.4f_%d)_D(%0.4f)_J(%0.4f)_gamma(%0.4f)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    num_trajectories, ...
    num_tp_periods, ...
    num_obs_periods, ...
    nph, ...
    ampl_start, ...
    ampl_shift, ...
    ampl_num, ...
    freq_start, ...
    freq_shift, ...
    freq_num, ...
    Delta, ...
    Omega, ...
    gamma, ...
    lpn_type, ...
    lpn_delta_s, ...
    lpn_delta_f_high, ...
    lpn_delta_f_low);

fn_fig = sprintf('%s/lambdas_%s', path_figures, suffix);
oqs_save_fig(fig, fn_fig);