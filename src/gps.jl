""" Utilities for integrating GPS with InSAR maps

Links:

1. list of LLH for all gps stations: ftp://gneiss.nbmg.unr.edu/rapids/llh
Note: ^^ This file is stored in the `STATION_LLH_FILE`

2. Clickable names to explore: http://geodesy.unr.edu/NGLStationPages/GlobalStationList
3. Map of stations: http://geodesy.unr.edu/NGLStationPages/gpsnetmap/GPSNetMap.html

"""

import CSV

const GPS_BASE_URL = "http://geodesy.unr.edu/gps_timeseries/tenv3/NA12/{station}.NA12.tenv3"

const STATION_LLH_FILE = joinpath(@__DIR__, "station_llh_all.csv")
const START_YEAR = 2014  # So far I don't care about older data

"""Returns the config folder for the application.  The default behavior
is to return whatever is most appropriate for the operating system.

This is used to store gps timeseries data

the following folders could be returned:
Mac OS X:
  ``~/Library/Application Support/apertools``
Mac OS X (POSIX):
  ``~/.apertools``
Unix:
  ``~/.cache/apertools``
Unix (POSIX):
  ``~/.apertools``

force_posix: if this is set to `true` then on any POSIX system the
    folder will be stored in the home folder with a leading
    dot instead of the XDG config home or darwin's
    application support folder.

Base on: https://github.com/pallets/click/blob/master/click/utils.py#L368
"""
function get_cache_dir(force_posix=false, app_name="gps")
    if force_posix
        path = joinpath(expanduser("~/.$app_name"))
    end
    # if Sys.isapple()
    #     path = joinpth(expanduser("~/Library/Application Support"), app_name)
    # end
    path = joinpath(
        get(ENV, "XDG_CONFIG_HOME", expanduser("~/.cache")),
        app_name
    )
    !exists(path) && mkdir(path)
    return path
end


"Read in a DataFrame of gps stations with columns name, lat, lon, alt"
read_station_llas() = CSV.read(STATION_LLH_FILE)
# columns = ["name", "lat", "lon", "alt"]


"Return the (lon, lat) of a `station_name`"
function station_lonlat(station_name)
    df = read_station_llas()
    !(station_name in df.name)&& error("No station named $station_name found.")
    # TODO:
    # closest_names = difflib.get_close_matches(station_name, df.name, n=5)
    # error("No station named %s found. Closest: %s" % (station_name, closest_names))
    name, lat, lon, alt = df[df.name .== station_name, :][1, :]
    return lon, lat
end

station_latlon(station_name) = reverse(station_lonlat(station_name))

station_rowcol(name, demrsc) = MapImages.latlon_to_rowcol(station_latlon(name)..., demrsc)


# function station_distance(station_name1, station_name2)
#     lonlat1 = station_lonlat(station_name1)
#     lonlat2 = station_lonlat(station_name2)
#     return apertools.latlon.latlon_to_dist(lonlat1[::-1], lonlat2[::-1])
# end


# function stations_within_rsc(rsc_filename=nothing, demrsc=nothing, gps_filename=nothing)
#     if demrsc is nothing
#         if rsc_filename is nothing
#             error("Need demrsc or rsc_filename")
#         end
#         demrsc = apertools.sario.load(rsc_filename)
#     end
# 
#     df = read_station_llas(filename=gps_filename)
#     station_lon_lat_arr = df[["lon", "lat"]].values
#     contains_bools = [apertools.latlon.grid_contains(s, **demrsc) for s in station_lon_lat_arr]
#     return df[contains_bools][["name", "lon", "lat"]].values
# end


function download_station_data(station_name)
    station_name = station_name.upper()
    url = GPS_BASE_URL.format(station=station_name)
    response = requests.get(url)

    stationfile = split(url, "/")[-1]  # Just blah.NA12.tenv3
    filename = "$GPS_DIR/$stationfile"
    println("Saving to $filename")

    write(filename, response.text)
end


function load_station_enu(station, start_year=START_YEAR, end_year=nothing, to_cm=true)
    """Loads one gps station's ENU data since start_year until end_year
    as separate Series items

    Will scale to cm by default, and center the first data point
    to 0
    """
    station_df = load_station_data(station, to_cm=to_cm, start_year=start_year, end_year=end_year)

    start_val = station_df[["east", "north", "up"]].iloc[:10].mean()
    enu_zeroed = station_df[["east", "north", "up"]] - start_val
    dts = station_df["dt"]
    return dts, enu_zeroed
end

function load_station_data(station_name,
                      start_year=START_YEAR,
                      end_year=nothing,
                      download_if_missing=true,
                      to_cm=true)
    """Loads one gps station's ENU data since start_year until end_year as a dataframe

    Args:
        station_name (str): 4 Letter name of GPS station
            See http://geodesy.unr.edu/NGLStationPages/gpsnetmap/GPSNetMap.html for map
        start_year (int), default 2014, cutoff for beginning of GPS data
        end_year (int): default nothing, cut off for end of GPS data
        download_if_missing (bool): default true
    """
    station_name = station_name.upper()
    gps_data_file = os.path.join(GPS_DIR, "%s.NA12.tenv3" % station_name)
    if not os.path.exists(gps_data_file)
        logger.warning("%s does not exist.", gps_data_file)
        if download_if_missing
            logger.info("Downloading %s", station_name)
            download_station_data(station_name)
        end
    end

    df = pd.read_csv(gps_data_file, header=0, sep=r"\s+")
    clean_df = _clean_gps_df(df, start_year, end_year)
    if to_cm
        # logger.info("Converting %s GPS to cm" % station_name)
        clean_df[["east", "north", "up"]] = 100 * clean_df[["east", "north", "up"]]
    end
    return clean_df
end


function _clean_gps_df(df, start_year=nothing, end_year=nothing)
    df["dt"] = pd.to_datetime(df["YYMMMDD"], format="%y%b%d")

    if start_year
        df_ranged = df[df["dt"] >= datetime.datetime(start_year, 1, 1)]
    end
    if end_year
        df_ranged = df_ranged[df_ranged["dt"] <= datetime.datetime(end_year, 1, 1)]
    end
    df_enu = df_ranged[["dt", "__east(m)", "_north(m)", "____up(m)"]]
    df_enu = df_enu.rename(mapper=s -> replace(replace(s, "_" => ""), "(m)" => ""), axis="columns")
    df_enu.reset_index(inplace=true, drop=true)
    return df_enu
end


# function stations_within_image(image_ll=nothing, filename=nothing, mask_invalid=true, gps_filename=nothing)
#     """Given an image, find gps stations contained in area
# 
#     Must have an associated .rsc file, either from being a LatlonImage (passed
#     by image_ll), or have the file in same directory as `filename`
# 
#     Args
#         image_ll (LatlonImage) LatlonImage of area with data
#         filename (str) filename to load into a LatlonImage
#         mask_invalid (bool) Default true. if true, don't return stations
#             where the image value is NaN or exactly 0
# 
#     Returns
#         ndarray Nx3, with columns ["name", "lon", "lat"]
#     """
#     if image_ll is nothing
#         image_ll = apertools.latlon.LatlonImage(filename=filename)
#     end
# 
#     df = read_station_llas(filename=gps_filename)
#     station_lon_lat_arr = df[[:lon, :lat]].values
#     contains_bools = image_ll.contains(station_lon_lat_arr)
#     candidates = df[contains_bools][["name", "lon", "lat"]].values
#     good_stations = []
#     if mask_invalid
#         for name, lon, lat in candidates
#             val = image_ll[..., lat, lon]
#             if np.isnan(val)  # or val == 0 TODO with window 1 reference, it's 0
#                 continue
#             else
#                 good_stations.append([name, lon, lat])
#             end
#         end
#     else
#         good_stations = candidates
#     end
#     return good_stations
# end

# function save_station_points_kml(station_iter)
#     for name, lat, lon, alt in station_iter
#         apertools.kml.create_kml(
#             title=name,
#             desc="GPS station location",
#             lon_lat=(lon, lat),
#             kml_out="%s.kml" % name,
#             shape="point",
#         )
#     end
# end
