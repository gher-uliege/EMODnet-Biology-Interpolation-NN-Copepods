#   Find the best results
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# The script performs a loop through all the results stored in a specified directory,
# reads the RMS computed by comparing the gridded field with the validation dataset,
# and report the best combination of parameters.
#
# Those parameters then have to be used with the full data files (i.e., without removing observations for validation)
# in order to create the final products.

using Pkg
Pkg.activate("../")
include("../scripts/CopepodsNN.jl")
include("../scripts/param.jl")
using NCDatasets
using Dates
using Glob
using DIVAnd
using Statistics
using JupyterFormatter
enable_autoformat()

#   Files and directories
#   =====================

# varname = "Small_copepods"
varname = "Large_copepods"

outputdir_optim = joinpath(outputdir, "param_optim/$(varname)")
isdir(outputdir_optim)
expdirlist = filter(x -> isdir(joinpath(outputdir_optim, x)), readdir(outputdir_optim));
@info("Working on $(length(expdirlist)) experiments")
@info("Files located in $(outputdir_optim)")

#   Find the best analysis
#   ======================

io = open(joinpath(outputdir, "results_$(varname).txt"), "w");
RMSlist = []
allfilelist = []
notfinished = 0
ndirproc = 0

for expdir in expdirlist
    @debug(expdir)
    resfilelist = filter(
        x -> isfile(joinpath(outputdir_optim, expdir, x)),
        readdir(joinpath(outputdir_optim, expdir)),
    )
    global ndirproc += 1
    if length(resfilelist) == 7
        filter!(xx -> .!occursin("param", xx), resfilelist)

        for resfile in resfilelist
            _, RMS = CopepodsNN.get_field(joinpath(outputdir_optim, expdir, resfile), varname)
            write(io, "$(resfile)\t $(RMS)\n")

            push!(RMSlist, RMS)
            push!(allfilelist, resfile)
        end
    else
        @debug("Experiment not completed")
        global notfinished += 1
    end
end
close(io);
@info(expdirlist[ndirproc]);

minRMS = minimum(RMSlist)
bestfile = allfilelist[argmin(RMSlist)];
@show(minRMS, bestfile)
expdirbest = split(split(bestfile, "_")[3], ".")[1];

#   Parameters of the best analysis
#   –––––––––––––––––––––––––––––––

bestfileparam = "../product/param_optim/$(varname)/$(expdirbest)/$(varname)_$(expdirbest)_params.nc"
ds = NCDataset(bestfileparam)
loss = ds["losses"][:]
@show(ds)
close(ds)