from enum import Enum

class FSType(Enum):
    local = 0
    cluster = 1
    mpipks_sd = 2
    mpipks_mv = 3

def get_root(type):

    root = ''
    if type == FSType.cluster:
        root = '/common/home/yusipov_i/Work/os_d/data/qjx'
    elif type == FSType.mpipks_sd:
        root = '/data/biophys/denysov/yusipov/os_d/data/qjx'
    elif type == FSType.mpipks_mv:
        root = '/data/condmat/ivanchen/yusipov/os_d/qjx'
    elif type == FSType.local:
        root = 'E:/Work/os_d/data/qjx'

    return root

def get_dir(path, ss, type):
    root = get_root(type)
    dir = root + path + '/ss_' + str(ss)
    return dir