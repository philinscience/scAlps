"""Small spatial data. Big landscapes."""

from .core import Terrain, terrain
from .datasets import demo
from .render import plot

__all__ = ["Terrain", "demo", "plot", "terrain"]
__version__ = "0.1.0"
