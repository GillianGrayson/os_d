from enum import Enum

class FSType(Enum):
    local_big = 0
    unn = 1
    mpipks = 2
    local_msi = 3

def get_root(type):

    root = ''
    if type is FSType.unn:
        root = '/common/home/yusipov_i/Work/os_d/data/jcs'
    elif type is FSType.mpipks:
        root = '/data/biophys/yusipov/os_d/data/jcs'
    elif type is FSType.local_big:
        root = 'E:/Work/os_d/data/jcs/mock'
    elif type is FSType.local_small:
        root = 'D:/Work/os_d/data/jcs/mock'

    return root

def get_dir(path, type):
    root = get_root(type)
    dir = root + path
    return dir