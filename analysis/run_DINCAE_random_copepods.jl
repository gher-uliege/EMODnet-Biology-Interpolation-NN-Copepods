#   Run DINCAE analysis
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# 
#   Run the simplest analysis with DINCAE:
# 
#     •  no environmental variables
# 
#     •  few time periods
# 
#     •  all the parameters with default values.

using Pkg
Pkg.activate("../")
Pkg.instantiate()
include("../scripts/CopepodsNN.jl")
using CUDA
using cuDNN
using DINCAE
using Dates
using NCDatasets

domain = (-95, 27.5, 22.5, 79.0)
Δext = 10.0
Δlon = 1.0
Δlat = 1.0
domaincompute = domain .+ [-Δext, Δext, -Δext, Δext]
yearmin = 1958
yearmax = 2022
fielddates_monthly = collect(Date(yearmin, 1, 15):Dates.Month(1):Date(yearmax, 12, 15));
dataprocdir = "../data/derived_data/"


#   Set the correct variable types
#   ================================
# 
#   (depending if CUDA using GPU)

const F = Float32

# Test if CUDA is functional to use the GPU, otherwise the CPU is used.

if CUDA.functional()
    Atype = CuArray{F}
else
    @warn "No supported GPU found. We will use the CPU which is very slow. Please check https://developer.nvidia.com/cuda-gpus"
    Atype = Array{F}
end

#   Domain
#   ========
# 
#   (defined in the configuration file config.jl)

@show domain;
grid = (domaincompute[1]:Δlon:domaincompute[2], domaincompute[3]:Δlat:domaincompute[4])

#   Time period
#   –––––––––––––
# 
#   (also defined in config.jl (config.jl))
# 
#   Here we create a monthly climatology, i.e. we don't consider the year of the
#   measurement. These values will be written in the netCDF input file.

@show fielddates_monthly;

#   Read the data file
#   =======================

inputfile = joinpath(dataprocdir, "Small_copepods_DINCAE.nc")

# Create output directory
# =======================

outputdir = "../product/param_optim"
isdir(outputdir) ? @debug("already there") : mkpath(outputdir)

# Download environmental data files
# =================================

# bathyfile = joinpath(covariabledir_regrid, "gebco_30sec_16_interp_1deg.nc")
# sstfile = joinpath(covariabledir_regrid, "sst_hadley_1deg_ufill.nc")

# CPRDINCAE.download_check(bathyfile, covariable_urls["bathymetry_1deg"])
# CPRDINCAE.download_check(sstfile, covariable_urls["sst_hadley_1deg"])

# # Path of the auxiliary fielddates
# auxdata_files = [
#   (filename = bathyfile, varname = "bathymetry", errvarname = "bathymetry_error"),
#   (filename = sstfile, varname = "SST", errvarname = "SST_error")
# ]



# Main loop
# =========

dateformat = "yyyymmddTHHMMSS"
for iii = 1:40
    timestamp = Dates.format(Dates.now(), dateformat)
    paramdict = CopepodsNN.create_random_params()
    @info(" ")
    @info(paramdict)

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp).nc")
    paramfile = joinpath(outputdir, "smallcopepods_$(timestamp)_params.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(F, Atype, inputfile, "Small_copepods", grid, [outputfile]; 
    paramdict..., paramfile=paramfile) # , auxdata_files=auxdata_files)

end
