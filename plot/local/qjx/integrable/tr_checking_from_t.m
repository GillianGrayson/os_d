clear all;

base_color = [0 0 1];
jump_color = [0 1 0];
renorm_color = [1 0 0];
linewidth = 2;

T = 1;

sys_id = 5;
task_id = 8;
prop_id = 0;

N = 7;
sys_size = 2^N;
seed = 1;
tau = 1;
k = -1;

lpn_type = -1;
delta_lpn = 1e-6;

deep_dump = 1;
deep_num_steps = 200;

num_target_traj = 5;
base_tr_id = 1;
var_tr_id = num_target_traj + base_tr_id;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(0_1000000)_N(%d)_seed(%d)_tau(%d)_k(%d)_T(%0.4f)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    N, ...
    seed, ...
    tau, ...
    k, ...
    T, ...
    lpn_type, ...
    log10(delta_lpn), ...
    log10(delta_lpn), ...
    log10(delta_lpn));

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(deep_dump == 1)
    num_periods = (num_dumps - 1) / (1 * deep_num_steps);
    dump_shift = 1 / (1 * deep_num_steps);
    dump_periods(1) = 0;
    for dump_id = 2:num_dumps
        dump_periods(dump_id) = dump_shift * (dump_id - 1);
    end
end
dump_periods = dump_periods * T;
num_periods = (size(dump_periods, 1) - 1) / deep_num_steps;
x_ticks = T * linspace(1, num_periods, num_periods)';
x_labels = {};
for i = 1 : num_periods
    x_labels{i} = sprintf('%d\\tau', i);
end

fn = sprintf('%s/norm_evo_%s.txt', data_path, suffix);
norm_evo_data = importdata(fn);
norm_evo = norm_evo_data(:, base_tr_id);

fn = sprintf('%s/jump_times_%d_%s.txt', data_path, base_tr_id - 1, suffix);
jump_times = importdata(fn);
jump_times = jump_times / T;
jump_times_after = jump_times + 10e-10;

fn = sprintf('%s/jump_norms_%d_%s.txt', data_path, base_tr_id - 1, suffix);
jump_norms = importdata(fn);
jump_norms_after = ones(size(jump_norms, 1), 1);

fn = sprintf('%s/jump_etas_%d_%s.txt', data_path, base_tr_id - 1, suffix);
jump_etas = importdata(fn);

all_times = vertcat(dump_periods, jump_times, jump_times_after);
all_norms = vertcat(norm_evo, jump_norms, jump_norms_after);

[all_times_sorted, order] = sort(all_times);
all_norms_sorted = zeros(size(all_times_sorted, 1), 1);

for t_id = 1:size(all_times_sorted, 1)
    all_norms_sorted(t_id) = all_norms(order(t_id));
end

fig = figure;
hLine = plot(all_times_sorted, all_norms_sorted);
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$||\psi||^2$', 'Interpreter', 'latex');
hold all;
hLine = scatter(jump_times, jump_norms);
propertyeditor(fig)

fn = sprintf('%s/phi_evo_%d_%s.txt', data_path, base_tr_id - 1, suffix);
phi_evo_data = importdata(fn);
phi_evo_base = zeros(sys_size, num_dumps);
norm_evo_base = zeros(num_dumps, 1);
for dump_id = 1:num_dumps
    s_id = (dump_id - 1) * sys_size + 1;
    f_id = dump_id * sys_size;
    phi_evo_base(:, dump_id) = phi_evo_data(s_id: f_id, 1) + 1i * phi_evo_data(s_id: f_id, 2);
    norm_evo_base(dump_id) = norm(phi_evo_base(:, dump_id)).^2;
end
fn = sprintf('%s/phi_evo_%d_%s.txt', data_path, var_tr_id - 1, suffix);
phi_evo_data = importdata(fn);
phi_evo_var = zeros(sys_size, num_dumps);
for dump_id = 1:num_dumps
    s_id = (dump_id - 1) * sys_size + 1;
    f_id = dump_id * sys_size;
    phi_evo_var(:, dump_id) = phi_evo_data(s_id: f_id, 1) + 1i * phi_evo_data(s_id: f_id, 2);
end
norm_diff = norm(norm_evo_base - norm_evo)

fn = sprintf('%s/spec_lpn_evo_%s.txt', data_path, suffix);
rnd_evo = importdata(fn);
fn = sprintf('%s/lambda_%s.txt', data_path, suffix);
lambda_data = importdata(fn);
fn = sprintf('%s/lambda_evo_%s.txt', data_path, suffix);
lambda_evo_data = importdata(fn);
lambda_evo = lambda_evo_data(:, num_target_traj + base_tr_id);
spec_evo_base = rnd_evo(:, 2 * (base_tr_id - 1) + 1);
spec_evo_var = rnd_evo(:, 2 * (base_tr_id - 1 + num_target_traj) + 1);
spec_evo_diff = abs(spec_evo_base - spec_evo_var);

fn = sprintf('%s/spec_mtx_%s.txt', data_path, suffix);
tmp = importdata(fn);
metric_mtx = zeros(sys_size);
for s_id_1 = 1:sys_size
    for s_id_2 = 1:sys_size
        index = (s_id_1 - 1) * sys_size + s_id_2;
        metric_mtx(s_id_1, s_id_2) = tmp(index, 1) + 1i * tmp(index, 2);
    end
end

spec_evo_base_check = zeros(num_dumps, 1);
for dump_id = 1:num_dumps
    curr_phi = phi_evo_base(:, dump_id) / sqrt(norm_evo_base(dump_id));
    tmp = curr_phi' * metric_mtx * curr_phi;
    spec_evo_base_check(dump_id) = tmp;
end

spec_evo_diff = norm(spec_evo_base_check - spec_evo_base)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = sprintf('%s/hamiltonian_%s.txt', data_path, suffix);
qjx_data = importdata(fn);
H = zeros(sys_size);
for s_id_1 = 1:sys_size
    for s_id_2 = 1:sys_size
        index = (s_id_1 - 1) * sys_size + s_id_2;
        H(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Dissipators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dissipators = {};
for diss_id = 1:N
    fn = sprintf('%s/dissipator_%d_%s.txt', data_path, diss_id - 1, suffix);
    qjx_data = importdata(fn);
    qjx_d = zeros(sys_size);
    for s_id_1 = 1:sys_size
        for s_id_2 = 1:sys_size
            index = (s_id_1 - 1) * sys_size + s_id_2;
            qjx_d(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
        end
    end
    dissipators{diss_id} = qjx_d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Effective Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_effective_matlab = H;
for diss_id = 1:N
    curr_diss = dissipators{diss_id};
    H_effective_matlab = H_effective_matlab - 0.5 * 1i * curr_diss' * curr_diss; 
end
H_effective_matlab = -1i * H_effective_matlab;
fn = sprintf('%s/hamiltonians_qj_0_%s.txt', data_path, suffix);
qjx_data = importdata(fn);
H_effective_qjx = zeros(sys_size);
for s_id_1 = 1:sys_size
    for s_id_2 = 1:sys_size
        index = (s_id_1 - 1) * sys_size + s_id_2;
        H_effective_qjx(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
    end
end
H_effective_diff = norm(H_effective_matlab - H_effective_qjx)

prop_mtx = expm(H_effective_matlab * 0.1);

init_state = complex(1/sqrt(sys_size) * ones(sys_size, 1),zeros(sys_size, 1));
next_state = prop_mtx * init_state;

% fn = sprintf('%s/random_obs_lpn_evo_1_%s.txt', data_path, suffix);
% rnd_obs_evo_data = importdata(fn);
% rnd_evo = zeros(size(rnd_obs_evo_data, 1), num_target_traj * 4);
% scan_template = repmat('(%e,%e)\t', 1, num_target_traj * 2);
% for line_id = 1:size(rnd_obs_evo_data, 1)
%     str = string(rnd_obs_evo_data(line_id));
%     data = sscanf(str, scan_template, num_target_traj * 4)';
%     rnd_evo(line_id, :) = data;
% end

fn = sprintf('%s/spec_lpn_evo_%s.txt', data_path, suffix);
rnd_evo = importdata(fn);

fn = sprintf('%s/lambda_%s.txt', data_path, suffix);
lambda_data = importdata(fn);

fn = sprintf('%s/lambda_evo_%s.txt', data_path, suffix);
lambda_evo_data = importdata(fn);

lambda_evo = lambda_evo_data(:, num_target_traj + base_tr_id);

spec_evo_base = rnd_evo(:, 2 * (base_tr_id - 1) + 1);
spec_evo_var = rnd_evo(:, 2 * (base_tr_id - 1 + num_target_traj) + 1);
spec_evo_diff = abs(spec_evo_base - spec_evo_var);

time_diff = dump_periods(2) - dump_periods(1);

fn = sprintf('%s/jump_times_%d_%s.txt', data_path, (base_tr_id - 1), suffix);
jump_data = importdata(fn);

interruptions = vertcat(jump_data, x_ticks);

fig = figure;
propertyeditor(fig);
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$\Delta(t)$', 'Interpreter', 'latex');
set(gca, 'YScale', 'log')
xtickangle(0)
xticks(x_ticks);
xticklabels(x_labels);
hold all;
curr_jump_id = 1;
prev_start = 1;
for p_id = 1:num_periods
    
    while(curr_jump_id <= size(jump_data, 1) && jump_data(curr_jump_id) < dump_periods(p_id * deep_num_steps))
        jump_time = jump_data(curr_jump_id);
        jump_time_id = floor(jump_time/time_diff);
        if(jump_time_id + 1 ~= (p_id-1) * deep_num_steps)
            plot_d = spec_evo_diff(prev_start:jump_time_id + 1);
            plot_t = dump_periods(prev_start:jump_time_id + 1);
            if (size(plot_t, 1) > 0)
                h = plot(plot_t, plot_d, 'LineWidth', linewidth, 'Color', base_color);
                flipud(h);
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                h = plot([jump_time jump_time], [spec_evo_diff(jump_time_id + 1) spec_evo_diff(jump_time_id + 2)], 'Color', jump_color, 'LineWidth', linewidth, 'LineStyle', ':');
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        else
            ololo = 1;
        end

        prev_start = jump_time_id + 2;
        curr_jump_id = curr_jump_id + 1;
    end
    
    plot_d = spec_evo_diff(prev_start:p_id * deep_num_steps);
    plot_t = dump_periods(prev_start:p_id * deep_num_steps);
    h = plot(plot_t, plot_d, 'LineWidth', linewidth, 'Color', base_color);
    flipud(h)
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    if (p_id < num_periods)
        h = plot([dump_periods(p_id * deep_num_steps) + 0.5 * time_diff, dump_periods(p_id * deep_num_steps) + 0.5 * time_diff], ...
            [spec_evo_diff(p_id * deep_num_steps), spec_evo_diff(p_id * deep_num_steps + 1)], ...
            'LineWidth', linewidth, 'Color', renorm_color, 'LineStyle', ':');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        prev_start = p_id * deep_num_steps + 1;
    end  
end

h = plot([dump_periods(1), dump_periods(end)], [delta_lpn, delta_lpn], 'Color', 'k', 'LineStyle', '--');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
box on;
hold all;
% yyaxis right;
% set(gca, 'FontSize', 30);
% ylabel('$L(t) = \frac{1}{t} \sum \ln{\frac{\Delta_f}{\Delta_s}}$', 'Interpreter', 'latex');
% h = plot(dump_periods, lambda_evo);
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% h = plot([dump_periods(1), dump_periods(end)], [0, 0], 'Color', 'k', 'LineStyle', '--');
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


lambda_avg = mean(lambda_data(num_target_traj+1:end));
title(sprintf('$\\langle \\lambda \\rangle_{%d \\mathrm{traj}} = %0.2e \\quad \\lambda=%0.2e$', num_target_traj, lambda_avg, lambda_data(num_target_traj + base_tr_id)), 'Interpreter', 'latex')

    

renorm_diffs_sum_log = sum(log(spec_evo_diff(200:200:end) / delta_lpn))
