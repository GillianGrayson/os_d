clear all;

T_part = 1;

tr_id = 1;

seed = 1;
mns = 1000000;

diss_type = 1;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 1.00; 
ps_drv_part_2 = 1.00; 
ps_drv_ampl = 1.75;
ps_prm_alpha = 5;
ps_prm_d = 10;
ps_prm_g = 10;
ps_diss_w = 0.05;
start_type = 0;
start_state = 49;

drv_T_1 = ps_drv_part_1 * T_part;
drv_T_2 = ps_drv_part_2 * T_part;
T = drv_T_1 + drv_T_2;

cd_dump_deep = 1;
cd_num_sub_steps = 50;

path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
    seed, ...
    mns, ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    diss_type, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ps_drv_ampl, ...
    ps_prm_alpha, ...
    ps_prm_d, ...
    ps_prm_g, ...
    start_type, ...
    start_state);
fn = sprintf('%s/periods_%s.txt', path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(cd_dump_deep == 1)
    num_periods = (num_dumps - 1) / (2 * cd_num_sub_steps);
    dump_shift_1 = drv_T_1 / cd_num_sub_steps;
    dump_shift_2 = drv_T_2 / cd_num_sub_steps;
    curr_id = 1;
    dump_periods(curr_id) = 0;
    for period_id = 1:num_periods
        for dump_id = 1:cd_num_sub_steps
            dump_periods(curr_id + 1) = dump_periods(curr_id) + dump_shift_1;
            curr_id = curr_id + 1;
        end
        
        for dump_id = 1:cd_num_sub_steps
            dump_periods(curr_id + 1) = dump_periods(curr_id) + dump_shift_2;
            curr_id = curr_id + 1;
        end
    end
    
    dump_periods = dump_periods / T;
end

fn = sprintf('%s/spec_2_evo_%s.txt', path, suffix);
evo_data = importdata(fn);
spec_3_evo = evo_data(:, 2 * tr_id - 1);

fn = sprintf('%s/mean_evo_%s.txt', path, suffix);
evo_data = importdata(fn);
mean_evo = evo_data(:, tr_id);

fig = figure;
yyaxis left;
hLine = plot(dump_periods, spec_3_evo, 'LineWidth', 2);
legend(hLine, sprintf('tr=%d', tr_id))
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$J_z = \frac{1}{2} \sum_{j=1}^{M} \sigma^z_j$', 'Interpreter', 'latex');
hold all;

yyaxis right;
hLine = plot(dump_periods, mean_evo, 'LineWidth', 2);
legend(hLine, sprintf('tr=%d', tr_id))
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
hold all;

propertyeditor(fig)



