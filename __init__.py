
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the APAC package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

'''
Usage
Example

import sys
path_root = 'C:/Users/ksbyeon/Desktop'

sys.path.append(path_root)
import apac
import os

path_out = os.path.join(os.path.dirname(apac.__file__), 'test')
fs_path = os.path.join(os.path.dirname(apac.__file__), 'fsaverage_LR32k')
core_maker = apac.core.core(path_out)
core_maker.call(fs_path)

core_maker.initial_roi('L')
core_maker.def_ppcore('L')

core_maker.initial_roi('R')
core_maker.def_ppcore('R')

'''


from . import mesh
from . import util
from . import core
