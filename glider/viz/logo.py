import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.patheffects
from matplotlib import transforms
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties

"""
Global variables for plotting logos
"""

fp = FontProperties(family="Arial", weight="bold") 
GLOBSCALE = 1.35
LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
            "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
            "A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
            "C" : TextPath((-0.366, 0), "C", size=1, prop=fp) }
COLOR_SCHEME = {'G': 'orange', 
                'A': 'red', 
                'C': 'blue', 
                'T': 'darkgreen'}

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)


def letterAt(letter, x, y, yscale=1, ax=None):
    text = LETTERS[letter]
    t = transforms.Affine2D().scale(1*GLOBSCALE, yscale*GLOBSCALE) + \
        transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p


def plot_logo(base_frequencies, ax=None):
    """
    base_frequencies = 
    """
    from glider.utils.motifs  import get_scores
    scores = get_scores(base_frequencies)

    if ax == None:
        fig, ax = plt.subplots(figsize=(10,3))

    x = 1
    maxi = 0
    for base_scores in scores:
        y = 0
        for base, score in zip("ACGT", base_scores):
            letterAt(base, x, y, score, ax)
            y += score
        x += 1
        maxi = max(maxi, y)

    plt.xticks(range(1,x))
    plt.xlim((0, x)) 
    plt.ylim((0, maxi)) 
    plt.tight_layout()      
    return ax
