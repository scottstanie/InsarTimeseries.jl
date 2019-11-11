import MapImages: MapImage
import CSV

geolist78 = Sario.load_geolist_from_h5("/Users/scott/Documents/Learning/masters-thesis/zoom_images78/stack_max700_alpha100.h5", "stack/1");
dfmt = dateformat"y-m-d H:M:S.s"
stack78 = MapImage("/Users/scott/Documents/Learning/masters-thesis/zoom_images78/stack_max700_alpha100.h5", "stack/1");
eqs15 = CSV.read("TXAR_Q12events_2015.csv", dateformat=dfmt);
eqs17 = CSV.read("TXAR_Q12events_2017.csv", dateformat=dfmt);
eqs17texnet = CSV.read("texnet_events_2017.csv", dateformat=dfmt);

@time include("/Users/scott/repos/InsarTimeseries.jl/scripts/plotting.jl");

animate_imgs_vs_pts(stack78, eqs17texnet, geolist78; latcol=:Latitude, loncol=:Longitude, sizecol=:Magnitude, datecol=:date, alpha=.2, vm=6, c="r", outname="animation_texnet78.gif")
animate_imgs_vs_pts(stack78, [eqs15; eqs17];, geolist78; latcol=:Latitude, loncol=:Longitude, sizecol=:Magnitude, datecol=:date, alpha=.2, vm=6, c="r", outname="animation_TXAR78.gif")
