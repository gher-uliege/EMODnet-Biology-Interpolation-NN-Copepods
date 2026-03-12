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
using LinearAlgebra

include("../scripts/param.jl")
varname = "Small_copepods"

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
#   (also defined in param.jl)

@show length(fielddates_monthly);

#   Read the data file
#   =======================

datafile1 = joinpath(dataprocdir, "Small_copepods_DINCAE_NortheastAtlantic_main.nc")
datafile2 = joinpath(dataprocdir, "Large_copepods_DINCAE_NortheastAtlantic_main.nc")
outputfilevalid1 = joinpath(datadir, replace(basename(datafile1), "_main.nc" => "_valid.nc"))
outputfilevalid2 = joinpath(datadir, replace(basename(datafile2), "_main.nc" => "_valid.nc"))

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

# Read the validation field

ds = NCDataset(outputfilevalid1, "r")
field_valid = ds["Small_copepods_mean"][:, :, :]
close(ds);

# Main loop
# =========

dateformat = "yyyymmddTHHMMSS"
for iii = 1:2
    timestamp = Dates.format(Dates.now(), dateformat)
    paramdict = CopepodsNN.create_random_params()
    @info(" ")
    @info(paramdict)

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp).nc")
    paramfile = joinpath(outputdir, "smallcopepods_$(timestamp)_params.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(F, Atype, datafile1, "Small_copepods", grid, [outputfile]; 
    paramdict..., paramfile=paramfile) # , auxdata_files=auxdata_files)

    # Compute RMSE
    ds = NCDataset(outputfile, "r")
    field = ds[varname][:, :, :]
    error = ds[varname * "_error"][:, :, :]
    close(ds);

    fieldiff = collect(skipmissing(field - field_valid))
    rmse = norm(fieldiff) / sqrt(length(fieldiff))
    bias = mean(fieldiff);

    @info("RMSE: $(rmse); bias: $(bias)");

    # Write the value into the netCDF
    ds = NCDataset(outputfile, "a")
    ds.attrib["RMSE"] = rmse
    ds.attrib["bias"] = bias
    close(ds);

end
