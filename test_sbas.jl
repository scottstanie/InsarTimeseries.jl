using InsarTimeseries
# igram_dir = "data/mask-test/"; ext=".unwflat"
igram_dir = "insar/tests/data/sbas_test/"; ext=".unw"

intlist = InsarTimeseries.sario.read_intlist(igram_dir);
geolist = InsarTimeseries.sario.read_geolist(igram_dir);
B = InsarTimeseries.build_B_matrix(geolist, intlist);
timediffs = InsarTimeseries.day_diffs(geolist);

unw_stack = InsarTimeseries.load_stack(directory=igram_dir, file_ext=ext);
unw_cols = InsarTimeseries.stack_to_cols(unw_stack);
