"""gueye : inversion semi-analytique de l'operateur de Laplace/Helmholtz discretise, 1D-3D, toutes conditions aux limites."""
from .kernel import TridiagInverse, SingularOperator
from .problem1d import Problem1D
from .wave1d import Wave1D
from .problem2d import Problem2D, Wave2D
from .problem3d import Problem3D, Wave3D
from .api import problem, wave, Solver, WaveSolver
from . import recipes
__version__ = "0.7.2"


def launch_gui(port=0, open_browser=True):
    """Ouvre l'interface graphique interactive (navigateur, serveur local) : voir gueye.gui."""
    from .gui import serve
    return serve(port, open_browser)
