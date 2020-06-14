cb = pyimport("palettable.colorbrewer.diverging")
cmgood = cb.RdYlBu_11_r.mpl_colormap
rdylbl = cb.RdYlBu_11_r.mpl_colormap


matlab_colors(n=4) = [[0, 0.4470, 0.7410, 1],
                   [0.8500, 0.3250, 0.0980, 1],
                   [0.9290, 0.6940, 0.1250, 1],
                   [0.4940, 0.1840, 0.5560, 1]][1:n]


SEISMIC_WIDE = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "seismic_wide",
    [
        (0, 0, 0.3, 1),
        (0, 0, 1, 1),
        (0.6, 0.6, 1, 1),
        (0.9, 0.9, 1, 1),
        (1, 1, 1, 1),
        (1, 0.9, 0.9, 1),
        (1, 0.6, 0.6, 1),
        (1, 0, 0, 1),
        (0.5, 0, 0, 1),
    ],  # Extra white in middle from seismic
    N = 250,
)
plt.register_cmap("seismic_wide", SEISMIC_WIDE)
SEISMIC_WIDE2 = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "seismic_wider",
    [
        (0, 0, 0.3, 1),
        (0, 0, 0.7, 1),
        (0.1, 0.1, 0.9, 1),
        (0.3, 0.3, 0.95, 1),
        (0.6, 0.6, 1, 1),
        (0.85, 0.85, 1, 1),
        (0.92, 0.92, 1, 0.99),
        (0.98, 0.98, 1, 0.98),
        (1, 1, 1, 0.95),
        (1, 0.98, 0.98, 0.98),
        (1, 0.92, 0.92, 0.99),
        (1, 0.85, 0.85, 1),
        (1, 0.6, 0.6, 1),
        (1, 0.3, 0.3, 1),
        (0.9, 0.1, 0.1, 1),
        (0.7, 0, 0, 1),
        (0.3, 0, 0, 1),
    ],
    N = 250,
)
plt.register_cmap("seismic_wider", SEISMIC_WIDE2)

SEISMIC_WIDE_Y = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "seismic_wide_y",
    [
        (0, 0, 0.7, 1),
        (0, 0, 1, 1),
        (0.1, 0.3, 1, 1),
        (0.2, 0.6, 0.95, 1),
        (0.5, 0.8, 0.85, 1),
        (0.8, 0.95, 0.80, 1),
        (0.9, 0.92, 0.78, 1),
        (0.95, 0.95, 0.75, 1),
        (0.95, 0.92, 0.70, 1),
        (0.95, 0.85, 0.65, 1),
        (0.97, 0.65, 0.3, 1),
        (1, 0.3, 0.2, 1),
        (1, 0.1, 0.05, 1),
        (1, 0, 0, 1),
        (0.7, 0, 0, 1),
    ],
    N = 450,
)
plt.register_cmap("seismic_wide_y", SEISMIC_WIDE_Y)


DSY5 = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "discrete_seismic_y5",
    broadcast(
        t -> t ./ 256,
        [
            (5, 113, 176, 256),
            (146, 197, 222, 256),
            # (247, 247, 247, 245),  # To make red-white-blue
            (247, 247, 191, 256),  # To make red-yellow-blue
            (244, 165, 130, 256),
            (202, 0, 32, 256),
        ],
    ),
    N = 5,
)
plt.register_cmap("discrete_seismic_y5", DSY5)

DSY7 = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "discrete_seismic_y7",
    broadcast(
        t -> t ./ 256,
        [
            (69, 117, 199, 256),
            (145, 191, 219, 256),
            (224, 243, 248, 256),
            (255, 255, 191, 256),
            (254, 224, 144, 256),
            (252, 141, 89, 256),
            (215, 48, 39, 256),
        ],
    ),
    N = 7,
)
plt.register_cmap("discrete_seismic_y7", DSY7)


DSY7 = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "discrete_seismic_y7",
    broadcast(
        t -> t ./ 256,
        [
            (69, 117, 199, 256),
            (145, 191, 219, 256),
            (224, 243, 248, 256),
            (255, 255, 191, 256),
            (254, 224, 144, 256),
            (252, 141, 89, 256),
            (215, 48, 39, 256),
        ],
    ),
    N = 7,
)
plt.register_cmap("discrete_seismic_y7", DSY7)

NSY = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "seismic_narrow_y",
    broadcast(
        t -> t ./ 256,
        [
            (9, 37, 99, 256),
            (39, 97, 199, 256),
            (105, 121, 219, 256),
            (255, 255, 191, 255),
            (252, 111, 69, 256),
            (215, 38, 19, 256),
            (105, 18, 9, 256),
        ],
    ),
    N = 256,
)
plt.register_cmap("seismic_narrow_y", NSY)
