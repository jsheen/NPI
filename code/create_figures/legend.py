#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 23:46:59 2021

@author: Justin Sheen

@description: script for plotting solely the legend of simulation figures.
"""

# Plot just the legend --------------------------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
from pathlib import Path
home = str(Path.home())
palette = dict(zip(['Control', 'Treatment'], ['orange', 'grey']))
handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]
plt.legend(handles=handles)
plt.gca().set_axis_off()
plt.tight_layout()
plt.savefig(home + '/NPI/code_output/figs/legend.png', dpi=300)