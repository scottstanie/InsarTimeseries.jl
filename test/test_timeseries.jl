using Dates

import InsarTimeseries: _get_pairs

split1 = convert(Array{Date,1}, [])
split2 = [Date(2017, 1, 1)]
split3 = [Date(2016, 1, 1), Date(2017, 1, 1)]

@test _get_pairs(split1) == (nothing, nothing)
@test _get_pairs(split2) == [(nothing, 2017 - 01 - 01), (2017 - 01 - 01, nothing)]
@test _get_pairs(split3) == [
    (nothing, 2016 - 01 - 01),
    (2016 - 01 - 01, 2017 - 01 - 01),
    (2016 - 01 - 01, nothing),
]
