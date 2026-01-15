# 
dataprocdir = "../data/derived_data/"


# Set domain extension
domain = (-30, 12., 42., 70.)
domain = (-95, 27.5, 22.5, 79.0)

# Set grid extension 
# (to avoid boundary problems)
Δext = 5.0
Δlon = 1.0
Δlat = 1.0
domaincompute = domain .+ [-Δext, Δext, -Δext, Δext]

# Set period of interest
yearmin = 1958
yearmax = 2022
fielddates_monthly = collect(Date(yearmin, 1, 15):Dates.Month(1):Date(yearmax, 12, 15));


longrid = domaincompute[1]:Δlon:domaincompute[2]
latgrid = domaincompute[3]:Δlat:domaincompute[4]

# Create grid for the validation data
Δlonvalid = 2.0
Δlatvalid = 2.0
longridvalid = domaincompute[1]:Δlonvalid:domaincompute[2]
latgridvalid = domaincompute[3]:Δlatvalid:domaincompute[4]