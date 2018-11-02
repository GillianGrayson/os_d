import socket
import getpass
import collections
import os


def get_root_path():
    root_path = ''
    host_name = socket.gethostname()
    if host_name == 'MSI':
        root_path = 'D:/YandexDisk/Work/os_d/data'
    elif host_name == 'DESKTOP-K9VO2TI':
        root_path = 'E:/YandexDisk/Work/os_d/data'
    elif host_name == 'DESKTOP-4BEQ7MS':
        root_path = 'D:/Aaron/Bio/os_d/data'
    elif host_name == 'master' or host_name[0:4] == 'node':
        user = getpass.getuser()
        root_path = '/common/home/' + user + '/Work/os_d/data'
    return root_path


def get_data_path(config):
    path = get_root_path() + '/' \
           + get_run_path(config) + '/' \
           + get_details_path(config) + '/' \
           + get_params_path(config) + '/' \
           + get_random_path(config) + '/'

    if not os.path.exists(path):
        os.makedirs(path)

    return path


def get_run_path(config):
    run = config.run
    path = 'run_' \
           + str(run.quantum_system.value) + '_' \
           + str(run.process.value) + '_' \
           + str(run.propagation.value)
    return path


def get_details_path(config):
    details = config.details
    path = 'details_' \
           + str(details.step_metrics) + '_' \
           + str(details.num_periods_trans) + '_' \
           + str(details.num_periods_obser)
    return path


def get_random_path(config):
    random = config.random
    path = 'random_' \
           + str(random.seed) + '_' \
           + str(random.num_seeds)
    return path


def get_dict_path(name, data):
    ordered_data = collections.OrderedDict(sorted(data.items()))
    path = name + '_' + '_'.join(list(map(str, ordered_data.values())))
    return path


def get_params_path(config):
    params = config.params
    path = get_dict_path('size', params.size) + '/' \
           + get_dict_path('dissipation', params.dissipation) + '/' \
           + get_dict_path('driving', params.driving)+ '/' \
           + get_dict_path('model', params.model)+ '/' \
           + get_dict_path('init_cond', params.init_cond)
    return path

