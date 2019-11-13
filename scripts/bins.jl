import MapImages
import MapImages: grid_extent, MapImage
import Sario
import CSV
import PyPlot; plt = PyPlot
using ImageFiltering
import Images: imresize
using PyCall; gps = pyimport("apertools.gps");
include("./map_plots.jl")

demrsc = Sario.load("dem_all.rsc")
water_prod = CSV.read("water_production.csv")
oil_prod = CSV.read("oil_production.csv")
water_inj = CSV.read("water_injection.csv")

waterinjbins = bin_vals(demrsc, water_inj)

waterprodbins_small = coarse_bin_vals(demrsc, water_prod, digits=1)
waterprodbins = imresize(waterprodbins_small, size(demrsc)...)
waterprodbins .*= (sum(bin_vals(demrsc, water_prod)) / sum(waterprodbins))

oilprodbins_small = coarse_bin_vals(demrsc, oil_prod, digits=1)
oilprodbins = imresize(oilprodbins_small, size(demrsc)...)
oilprodbins .*= (sum(bin_vals(demrsc, oil_prod)) / sum(oilprodbins))

netbins = waterinjbins .- waterprodbins .- oilprodbins

eqs = CSV.read("texnet_events.csv")
eq_bins = bin_vals(demrsc, eqs, loncol=:Longitude, latcol=:Latitude, valcol=:Magnitude)

fracs_all = CSV.read("FracDetails_Dates_TotalVolume.csv")
# frac_bins = bin_vals(demrsc, fracs, loncol=:mid_longitude, latcol=:mid_latitude, valcol=:total_volume_gals)
fracs = CSV.read("Frac_2015_simple.csv")
frac_bins = bin_vals(demrsc, fracs, loncol=:lon, latcol=:lat, valcol=:volume)


im78 = MapImage("../images_path78/velocities_2018_linear_max700.h5", "velos/1");
binsblur = imfilter(netbins, Kernel.gaussian(2))
eqsblur = imfilter(eq_bins, Kernel.gaussian(1))
fracsblur = imfilter(frac_bins, Kernel.gaussian(2))

fig, axes = plt.subplots(2, 2, sharex=true, sharey=true)
axim = axes[1, 1].imshow(im78, cmap="seismic", vmin=-20, vmax=20, extent=grid_extent(im78)); fig.colorbar(axim, ax=axes[1,1])
axim = axes[1, 2].imshow(binsblur, vmax=1e6, vmin=-1e6, cmap="seismic", extent=grid_extent(demrsc)); fig.colorbar(axim, ax=axes[1,2])
axim = axes[2, 1].imshow(fracsblur, vmax=2e7, cmap="Reds", extent=grid_extent(demrsc)); fig.colorbar(axim, ax=axes[2,1])
axim = axes[2, 2].imshow(eqsblur, vmax=5, cmap="Reds", extent=grid_extent(demrsc)); fig.colorbar(axim, ax=axes[2,2])


axes[1, 1].set_title("InSAR")
axes[1, 2].set_title("Water Inj. - Water Prod. - Oil Prod")
axes[2, 1].set_title("Fracs")
axes[2, 2].set_title("Eqs")
# axes[1, 2].imshow(netbins, vmax=1e6, vmin=-1e6, cmap="seismic")
