from __future__ import division, print_function, absolute_import

import pkg_resources

__version__ = pkg_resources.get_distribution('pulse_world').version

from .pyworld import *
