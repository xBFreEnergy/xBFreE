
# __all__ = ['alamdcrd', 'amber_outputs', 'analyzer', 'API', 'app', 'calculation', 'commandlineparser', 'createinput',
#            'exceptions', 'infofile', 'input_parser', 'main', 'make_top', 'make_trajs',
#            'output_file', 'parm_setup', 'timer', 'utils', '__version__', '__mmpbsa_version__', '__ambertools_version__']

__author__ = "Mario S. Valdes Tresanco, Mario E. Valdes Tresanco, Pedro A. Valiente PhD and Ernesto Moreno PhD"
__license__ = "GPLv3"
__mmpbsa_author__ = "Jason Swails, Dwight McGee, and Bill Miller III"
__mmpbsa_version__ = "16.0"
__ambertools_version__ = "20"
__prog__ = 'xBFreE'

from xBFreE._version import get_versions
__version__ = get_versions()['version']
del get_versions