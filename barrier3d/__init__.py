import pkg_resources

from .barrier3d import Barrier3d
from .bmi import Barrier3dBmi
from .configuration import Barrier3dConfiguration
from .load_input import load_inputs

__version__ = pkg_resources.get_distribution("barrier3d").version
__all__ = ["Barrier3d", "Barrier3dBmi", "Barrier3dConfiguration", "load_inputs"]

del pkg_resources
