#   Run DINCAE analysis
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
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
outputfilevalid1 =
    joinpath(datadir, replace(basename(datafile1), "_main.nc" => "_valid.nc"))
outputfilevalid2 =
    joinpath(datadir, replace(basename(datafile2), "_main.nc" => "_valid.nc"))

# Create output directory
# =======================

outputbasedir = "../product/param_optim"
isdir(outputbasedir) ? @debug("already there") : mkpath(outputbasedir)

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

    outputdir = joinpath(outputbasedir, timestamp)
    mkpath(outputdir)

    # 1. Base run
    # -----------

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp).nc")
    paramfile = joinpath(outputdir, "smallcopepods_$(timestamp)_params.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(
        F,
        Atype,
        datafile1,
        "Small_copepods",
        grid,
        [outputfile];
        paramdict...,
        paramfile = paramfile,
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

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp)_bathy.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(
        F,
        Atype,
        datafile1,
        "Small_copepods",
        grid,
        [outputfile];
        paramdict...,
        paramfile = paramfile,
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

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp)_sst.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(
        F,
        Atype,
        datafile1,
        "Small_copepods",
        grid,
        [outputfile];
        paramdict...,
        paramfile = paramfile,
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

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp)_bathy_sst.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(
        F,
        Atype,
        datafile1,
        "Small_copepods",
        grid,
        [outputfile];
        paramdict...,
        paramfile = paramfile,
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

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp)_bathy_sst_coast.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(
        F,
        Atype,
        datafile1,
        "Small_copepods",
        grid,
        [outputfile];
        paramdict...,
        paramfile = paramfile,
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

    outputfile = joinpath(outputdir, "smallcopepods_$(timestamp)_coast.nc")
    @info("Writing results to file $(outputfile)\nand parameters to file $(paramfile)")

    theloss = DINCAE.reconstruct_points(
        F,
        Atype,
        datafile1,
        "Small_copepods",
        grid,
        [outputfile];
        paramdict...,
        paramfile = paramfile,
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


end
