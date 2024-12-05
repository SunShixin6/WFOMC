from enum import Enum

from .IncrementalWFOMC_teacher import incremental_wfomc_new
from .IncrementalWFOMC_new_realdp import incremental_wfomc_new2
from .StandardWFOMC import standard_wfomc
from .FastWFOMC import fast_wfomc
from .IncrementalWFOMC import incremental_wfomc
from .RecursiveWFOMC import recursive_wfomc

__all__ = [
    "standard_wfomc",
    "fast_wfomc",  # 不用了解
    "incremental_wfomc",  # 不用了解
    "incremental_wfomc_new",  # 不用了解
    "incremental_wfomc_new2",
    "recursive_wfomc"
]


class Algo(Enum):
    STANDARD = 'standard'
    FAST = 'fast'
    FASTv2 = 'fastv2'
    INCREMENTAL = 'incremental'
    INCREMENTAL_NEW = 'incremental_new'
    INCREMENTAL_NEW2 = 'incremental_new2'
    RECURSIVE = 'recursive'
    DFT = "dft" #这里表示使用dft
    DFT_V = "dft_v"

    def __str__(self):
        return self.value
