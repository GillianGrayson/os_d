clear all;

T_part = 1;

tr_id = 1;

seed = 0;
mns = 1000000;

diss_type = 1;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 0.98; 
ps_drv_part_2 = 1.00; 
ps_drv_ampl = 3.2;
ps_prm_alpha = 5;
ps_prm_d = 1.0;
ps_prm_g = 1.0;
ps_diss_w = 0.05;
start_type = 0;
start_state = 0;

drv_T_1 = ps_drv_part_1 * T_part;
drv_T_2 = ps_drv_part_2 * T_part;
T = drv_T_1 + drv_T_2;

cd_dump_deep = 0;
cd_num_sub_steps = 100;

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



fn = sprintf('%s/spec_evo_%s.txt', path, suffix);
spec_evo_data = importdata(fn);

spec_evo = complex(spec_evo_data(:, 2 * tr_id - 1), spec_evo_data(:, 2 * tr_id));


fig = figure;

hLine = plot(dump_periods, real(spec_evo));
legend(hLine, sprintf('tr=%d', tr_id))
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$Re(\theta)$', 'Interpreter', 'latex');
hold all;

propertyeditor(fig)

fig = figure;

hLine = plot(dump_periods, imag(spec_evo));
legend(hLine, sprintf('tr=%d', tr_id))
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$Im(\theta)$', 'Interpreter', 'latex');
hold all;

propertyeditor(fig)


