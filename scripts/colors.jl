SEISMIC_WIDE = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "seismic_wide",
    [(0, 0, .3, 1), (0, 0, 1, 1), (.6, .6, 1, 1), 
     (.9, .9, 1, 1), (1, 1, 1, 1), (1, .9, .9, 1),
     (1, 0.6, 0.6, 1), (1, 0, 0, 1), (.5, 0, 0, 1)],  # Extra white in middle from seismic
    N=250,
)
plt.register_cmap("seismic_wide", SEISMIC_WIDE)
SEISMIC_WIDE2 = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "seismic_wider",
    [(0, 0, .3, 1),
     (0, 0, .7, 1),
     (.1, .1, .9, 1),
     (.3, .3, .95, 1), 
     (.6, .6, 1, 1), 
     (.85, .85, 1, 1), 
     (.92, .92, 1, .9), 
     (.98, .98, 1, .8), 
     (1, 1, 1, .7), 
     (1, .98, .98, .8), 
     (1, .92, .92, .9),
     (1, .85, .85, 1),
     (1, 0.6, 0.6, 1), 
     (1, 0.3, 0.3, 1), 
     (.9, .1, .1, 1),
     (.7, 0, 0, 1),
     (.3, 0, 0, 1)],
    N=250,
)
plt.register_cmap("seismic_wider", SEISMIC_WIDE2)

SEISMIC_WIDE_Y = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "seismic_wide_y",
    [
     (0, 0, .7, 1),
     (0, 0, 1, 1),
     (.1, .3, 1, 1), 
     (.2, .6, .95, 1), 
     (.5, .8, .85, 1), 
     (.8, .95, .80, 1), 
     (.9, .92, .78, .9),
     (.95, .95, .75, .8),
     (.95, .92, .70, .9),
     (.95, .85, .65, 1),
     (.97, .65, .3, 1),
     (1, 0.3, 0.2, 1), 
     (1, 0.1, 0.05, 1), 
     (1, 0, 0, 1),
     (.7, 0, 0, 1),
    ],
    N=450,
)
plt.register_cmap("seismic_wide_y", SEISMIC_WIDE_Y)
