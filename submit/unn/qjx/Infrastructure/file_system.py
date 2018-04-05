from enum import Enum

class FSType(Enum):
    local = 0
    cluster = 1

def get_root(type):

    root = ''
    if type is FSType.cluster:
        root = '\\home\\yusipov_i\\Work\\os_d\\qjx'
    elif type is FSType.local:
        root = 'E:\\Work\\os_d\\data\\qjx'

    return root

def get_dir(path, ss, type):
    root = get_root(type)
    dir = root + path + '\\ss_' + str(ss)
    return dir