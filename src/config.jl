import TOML


# Just to identify the output for a run, use the current time
df = Dates.format(Dates.now(), DateFormat("YYYY-mm-ddTHH:MM:SS"))
_default_outfile() = "deformation_" * Dates.format(Dates.now(), df) * ".h5"

# TODO: figure out how to incorporate the argparse checking
# with this TOML parsing

DEFAULT_KEYS = Dict(
    ### Files used in solving:
    
    "unw_stack_file" => UNW_FILENAME,
    "input_dset" => STACK_FLAT_SHIFTED_DSET,

    "outfile" => _default_outfile(),
    # This will be the dataset in the output file
    "outdset" => "velos1",

    # name of .h5 file with mask data
    "mask_stack_file" => "masks.h5",

    "ignore_geo_file" => "",

    # If we want to add specific extra, list comma separated
    # TODO: implement or remove
    "ignore_geo_dates" => "",

    "max_temporal_baseline" => 500,

    ### Type of solver ###
    "constant_velocity" => true,
    # Instead of SBAS: just average all igrams by date in a stack
    "stack_average" => false,

    ### Extra solver options ### 
    # Use L1 norm for SBAS cost instead of L2 least squares
    "L1" => true,
    # (SBAS) Strength of Tikhonov regularization
    "alpha" => 0.0,
)

