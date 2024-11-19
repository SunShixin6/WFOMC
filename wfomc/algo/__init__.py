from enum import Enum

from .StandardWFOMC import standard_wfomc
from .FastWFOMC import fast_wfomc
from .IncrementalWFOMC import incremental_wfomc
from .RecursiveWFOMC import recursive_wfomc

__all__ = [
    "standard_wfomc",
    "fast_wfomc",  # 不用了解
    "incremental_wfomc",  # 不用了解
    "recursive_wfomc"
]


class Algo(Enum):
    STANDARD = 'standard'
    FAST = 'fast'
    FASTv2 = 'fastv2'
    INCREMENTAL = 'incremental'
    RECURSIVE = 'recursive'
    DFT = "dft" #这里表示使用dft
    DFT_V = "dft_v"

    def __str__(self):
        return self.value
