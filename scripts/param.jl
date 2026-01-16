# -------------------------------------------
# Parameter for the configuration of the runs
# -------------------------------------------

using Dates

timesubsetting = false

# Files and directories
rawdatadir = "../data/raw_data/"
dataprocdir = "../data/derived_data/"
datadir = "../data/derived_data/"
auxdatadir = "../data/aux_data/"
validdir = "../data/derived_data/validation/"
figdir = "../product/figures/"

mkpath.([auxdatadir, datadir, dataprocdir, validdir, figdir]);

# Data files
rawdatafile = joinpath(rawdatadir, "CPR_DINCAE_Data_290925.csv")

# Set domain extension
regiondict = Dict(
    "NortheastAtlantic" => (-30, 9., 42., 67.),
    "NorthAtlantic" => (-95, 27.5, 22.5, 79.0)
)
regionname = "NortheastAtlantic"

# Set grid extension 
# (to avoid boundary problems)
Δext = 5.0
Δlon = 1.0
Δlat = 1.0
domain = regiondict[regionname]
domaincompute = domain .+ (-Δext, Δext, -Δext, Δext)

# Set period of interest
yearmin = 1958
yearmax = 2022
fielddates_monthly = collect(Date(yearmin, 1, 15):Dates.Month(1):Date(yearmax, 12, 15));

# Grid for the interpolation
longrid = domaincompute[1]:Δlon:domaincompute[2]
latgrid = domaincompute[3]:Δlat:domaincompute[4]

# Grid for the validation data
Δlonvalid = 2.0
Δlatvalid = 2.0
longridvalid = domaincompute[1]:Δlonvalid:domaincompute[2]
latgridvalid = domaincompute[3]:Δlatvalid:domaincompute[4]