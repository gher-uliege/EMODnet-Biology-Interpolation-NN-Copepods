#   Run DINCAE analysis using the best parameter combination
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# 
using Pkg
Pkg.activate("../")
Pkg.instantiate()
include("../scripts/CopepodsNN.jl")
include("../scripts/param.jl")
using CUDA
using cuDNN
using DINCAE
using Dates
using NCDatasets
using LinearAlgebra
using Statistics

# varname = "Small_copepods"
# expdirbest = "20260325T184717"

varname = "Large_copepods"
expdirbest = "20260408T073427"

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

datafile = joinpath(dataprocdir, "$(varname)_DINCAE_NortheastAtlantic.nc")
paramfile = joinpath("../product/param_optim/$(varname)", expdirbest, "$(varname)_$(expdirbest)_params.nc")
outputdir = "../product/"


outputfilevalid1 =
    joinpath(datadir, replace(basename(datafile), ".nc" => "_valid.nc"))
ds = NCDataset(outputfilevalid1, "r")
field_valid = ds["$(varname)_mean"][:, :, :]
close(ds);


# Read parameters from file
# =========================

function read_params(paramfile::AbstractString)
    ds = NCDataset(paramfile)
    

    paramdictopti = Dict(
        :learning_rate                 => ds.attrib["learning_rate"],
        :epochs                        => ds.attrib["epochs"],
        :laplacian_penalty             => ds.attrib["laplacian_penalty"]
    )

    close(ds)

    return paramdictopti
end

paramdictopti = read_params(paramfile)


# 1. Base run
# -----------

outputfile = joinpath(outputdir, "$(varname)_best.nc")
paramfileout = joinpath(outputdir, "$(varname)_best_params.nc")

theloss = DINCAE.reconstruct_points(
    F,
    Atype,
    datafile,
    "$(varname)",
    grid,
    [outputfile];
    paramdictopti...,
    paramfile = paramfileout
) # , auxdata_files=auxdata_files)

# Compute RMSE
nc1 = NCDataset(outputfile, "r")
field = nc1[varname][:, :, :]
error = nc1[varname*"_error"][:, :, :]
close(nc1);

fieldiff = collect(skipmissing(field - field_valid))
rmse = norm(fieldiff) / sqrt(length(fieldiff))
bias = mean(fieldiff);

@info("RMSE: $(rmse); bias: $(bias)");

# Write the value into the netCDF
nc2 = NCDataset(outputfile, "a")
nc2.attrib["RMSE"] = rmse
nc2.attrib["bias"] = bias
close(nc2);

# 2. Run with bathymetry
# ----------------------

auxdata_files = [(
    filename = joinpath(dataprocdir, "gebco_30sec_16_interp_1deg_time.nc"),
    varname = "bat",
    errvarname = "bat_error",
),]

outputfile = joinpath(outputdir, "$(varname)_best_bathy.nc")
paramfileout = joinpath(outputdir, "$(varname)_best_bathy_params.nc")

@info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

theloss = DINCAE.reconstruct_points(
    F,
    Atype,
    datafile,
    "$(varname)",
    grid,
    [outputfile];
    paramdictopti...,
    paramfile = paramfileout,
    auxdata_files = auxdata_files,
)

# Compute RMSE
nc1 = NCDataset(outputfile, "r")
field = nc1[varname][:, :, :]
error = nc1[varname*"_error"][:, :, :]
close(nc1);

fieldiff = collect(skipmissing(field - field_valid))
rmse = norm(fieldiff) / sqrt(length(fieldiff))
bias = mean(fieldiff);

@info("RMSE: $(rmse); bias: $(bias)");

# Write the value into the netCDF
nc2 = NCDataset(outputfile, "a")
nc2.attrib["RMSE"] = rmse
nc2.attrib["bias"] = bias
close(nc2);

# 3. Run with SST
# ---------------

auxdata_files = [(
    filename = joinpath(dataprocdir, "hadley_sst_1deg_time.nc"),
    varname = "SST",
    errvarname = "SST_error",
),]

outputfile = joinpath(outputdir, "$(varname)_best_sst.nc")
paramfileout = joinpath(outputdir, "$(varname)_best_sst_params.nc")

@info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

theloss = DINCAE.reconstruct_points(
    F,
    Atype,
    datafile,
    "$(varname)",
    grid,
    [outputfile];
    paramdictopti...,
    paramfile = paramfileout,
    auxdata_files = auxdata_files,
)

# Compute RMSE
nc1 = NCDataset(outputfile, "r")
field = nc1[varname][:, :, :]
error = nc1[varname*"_error"][:, :, :]
close(nc1);

fieldiff = collect(skipmissing(field - field_valid))
rmse = norm(fieldiff) / sqrt(length(fieldiff))
bias = mean(fieldiff);

@info("RMSE: $(rmse); bias: $(bias)");

# Write the value into the netCDF
nc2 = NCDataset(outputfile, "a")
nc2.attrib["RMSE"] = rmse
nc2.attrib["bias"] = bias
close(nc2);

# 4. Run with bathymetry and SST
# ------------------------------

auxdata_files = [
    (
        filename = joinpath(dataprocdir, "gebco_30sec_16_interp_1deg_time.nc"),
        varname = "bat",
        errvarname = "bat_error",
    ),
    (
        filename = joinpath(dataprocdir, "hadley_sst_1deg_time.nc"),
        varname = "SST",
        errvarname = "SST_error",
    ),
]

outputfile = joinpath(outputdir, "$(varname)_best_bathy_sst.nc")
paramfileout = joinpath(outputdir, "$(varname)_best_bathy_sst_params.nc")

@info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

theloss = DINCAE.reconstruct_points(
    F,
    Atype,
    datafile,
    "$(varname)",
    grid,
    [outputfile];
    paramdictopti...,
    paramfile = paramfileout,
    auxdata_files = auxdata_files,
)

# Compute RMSE
nc1 = NCDataset(outputfile, "r")
field = nc1[varname][:, :, :]
error = nc1[varname*"_error"][:, :, :]
close(nc1);

fieldiff = collect(skipmissing(field - field_valid))
rmse = norm(fieldiff) / sqrt(length(fieldiff))
bias = mean(fieldiff);

@info("RMSE: $(rmse); bias: $(bias)");

# Write the value into the netCDF
nc2 = NCDataset(outputfile, "a")
nc2.attrib["RMSE"] = rmse
nc2.attrib["bias"] = bias
close(nc2);

# 5. Run with bathymetry and SST and distance to coast
# ----------------------------------------------------

auxdata_files = [
    (
        filename = joinpath(dataprocdir, "gebco_30sec_16_interp_1deg_time.nc"),
        varname = "bat",
        errvarname = "bat_error",
    ),
    (
        filename = joinpath(dataprocdir, "hadley_sst_1deg_time.nc"),
        varname = "SST",
        errvarname = "SST_error",
    ),
    (
        filename = joinpath(dataprocdir, "dist_coast_interp_1deg_time.nc"),
        varname = "dist",
        errvarname = "dist_error",
    ),
]

outputfile = joinpath(outputdir, "$(varname)_best_bathy_sst_coast.nc")
paramfileout = joinpath(outputdir, "$(varname)_best_bathy_sst_coast_params.nc")

@info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

theloss = DINCAE.reconstruct_points(
    F,
    Atype,
    datafile,
    "$(varname)",
    grid,
    [outputfile];
    paramdictopti...,
    paramfile = paramfileout,
    auxdata_files = auxdata_files,
)

# Compute RMSE
nc1 = NCDataset(outputfile, "r")
field = nc1[varname][:, :, :]
error = nc1[varname*"_error"][:, :, :]
close(nc1);

fieldiff = collect(skipmissing(field - field_valid))
rmse = norm(fieldiff) / sqrt(length(fieldiff))
bias = mean(fieldiff);

@info("RMSE: $(rmse); bias: $(bias)");

# Write the value into the netCDF
nc2 = NCDataset(outputfile, "a")
nc2.attrib["RMSE"] = rmse
nc2.attrib["bias"] = bias
close(nc2);

# 6. Run with distance to coast
# ----------------------------------------------------

auxdata_files = [
    (
        filename = joinpath(dataprocdir, "dist_coast_interp_1deg_time.nc"),
        varname = "dist",
        errvarname = "dist_error",
    ),
]

outputfile = joinpath(outputdir, "$(varname)_best_coast.nc")
paramfileout = joinpath(outputdir, "$(varname)_best_coast_params.nc")

@info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

theloss = DINCAE.reconstruct_points(
    F,
    Atype,
    datafile,
    "$(varname)",
    grid,
    [outputfile];
    paramdictopti...,
    paramfile = paramfileout,
    auxdata_files = auxdata_files,
)

# Compute RMSE
nc1 = NCDataset(outputfile, "r")
field = nc1[varname][:, :, :]
error = nc1[varname*"_error"][:, :, :]
close(nc1);

fieldiff = collect(skipmissing(field - field_valid))
rmse = norm(fieldiff) / sqrt(length(fieldiff))
bias = mean(fieldiff);

@info("RMSE: $(rmse); bias: $(bias)");

# Write the value into the netCDF
nc2 = NCDataset(outputfile, "a")
nc2.attrib["RMSE"] = rmse
nc2.attrib["bias"] = bias
close(nc2);

