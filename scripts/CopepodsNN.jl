module CopepodsNN

using DelimitedFiles
using Dates
using NCDatasets
using DataStructures
using GridInterpolations
using Interpolations
using Statistics
using Random
using Distributions
using GeoArrays
using ImageFiltering
using GeoJSON
using PolygonOps
using Downloads
using GeoDatasets
using Makie 
using CairoMakie
using GeoMakie

juliaversion = "v1.12.4"

"""
    get_dox_url(ID)

Return the full URL to download one of the files stored on DoX,
the ULiège internal file storage service
"""
function get_dox_url(ID::AbstractString)
    return "https://dox.uliege.be/index.php/s/$(ID)/download"
end

covariable_urls = Dict(
    "sst_terra_1deg" => get_dox_url("BAbFOurNvbAlCMz"),
    "sst_reynolds_1deg" => get_dox_url("9CLtnsCSmciXdK1"),
    "sst_hadley_1deg" => get_dox_url("uq1BolcW3FIwAnf"),
    "dist2coast_1deg" => get_dox_url("ladbtZHwI9e0ZNM"),
    "bathymetry_1deg" => get_dox_url("fuTcDgXJMzx2xkv"),
    "bathymetry_1deg_19932020" => get_dox_url("McFJ1Eeatjb9ICg"),
    "gebco_30sec_16" => get_dox_url("U0pqyXhcQrXjEUX"),
    "gebco_30sec_08" => get_dox_url("wS6Y8P8NhIF60eG"),
    "gebco_30sec_04" => get_dox_url("RSwm4HPHImdZoQP"),
    "chloro_cmems" => get_dox_url("5qtJuEPedq6umSO") 
)

"""
    download_check(datafile, datafileURL)

Check if the file `datafile` exists, otherwise download it from DoX `datafileURL`.

# Example
```julia-repl
julia> download_check("gebco_30sec_16.nc", "https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download")
```
"""
function download_check(datafile::AbstractString, datafileURL::AbstractString)
    if isfile(datafile)
        @info("File $(datafile) already downloaded")
    else
        @info("Downloading file $(datafile)")
        Downloads.download(datafileURL, datafile)
    end
end

"""
    read_polygon_json(contourfile)

Read the coordinates as a list of tuples stored in the geoJSON file `contourfile`,
as downloaded from https://geojson.io

# Example
```julia-repl
julia> coordlist = read_polygon_json(contourfile)
```
"""
function read_polygon_json(contourfile::AbstractString)
    coordlist = []
    jsonbytes = read(contourfile);
    fc = GeoJSON.read(jsonbytes)
    for poly in fc
        coordinates = poly.geometry[1]
        push!(coordlist, coordinates)
    end
    return coordlist
end

"""
    edit_mask!(xi, yi, mask, coordinatelist)

Edit the land-sea mask (as read using `DIVAnd.load_bath`) by setting to zero (_land_ value)
the cells that fall inside the polygon(s) defined by `coordinatelist` (as read by using `read_polygon_json`).

The contour file can be downloaded from geojson.io. An example of such a file (`outsidemask.json`) is provided in 
the data directory.

# Example
```julia-repl
julia> coordinatelist = read_polygon_json(contourfile);
julia> edit_mask!(xi, yi, mask, coordinatelist)
```
"""
function edit_mask!(xi, yi, mask::BitMatrix, coordinatelist::Vector{Any})
    for (ii, xx) in enumerate(xi)
        for (jj, yy) in enumerate(yi)
            for coordinates in coordinatelist
                if PolygonOps.inpolygon((xx,yy), coordinates) != 0
                    mask[ii, jj] = 0.
                end
            end
        end
    end
    return nothing
end

"""
    read_copepods_csv(datafile, sorting=true)

Read the coordinates and observations from the original CSV file

# Examples
```julia-repl
julia> lon, lat, dates, Total_Copepods_Large, Total_Copepods_Small, SampleId = read_data_calanus(datafile)
julia> lon, lat, dates, Total_Copepods_Large, Total_Copepods_Small, SampleId = read_data_calanus(datafile, sorting=false)
```
"""
function read_copepods_csv(datafile::AbstractString; sorting::Bool=true)

    # Read the first line to get column names and indices
    column_names = split(readline(datafile), ",")
    index_longitude = findfirst(isequal("Longitude"), column_names)
    index_latitude = findfirst(isequal("Latitude"), column_names)
    index_year = findfirst(isequal("Year"), column_names)
    index_month = findfirst(isequal("Month"), column_names)
    index_day = findfirst(isequal("Day"), column_names)
    index_id = findfirst(isequal("SampleId"), column_names)
    index_coplarge = findfirst(isequal("Small_Copepods"), column_names)
    index_copsmall = findfirst(isequal("Large_Copepods"), column_names)
   
    # Read all the data in one go
    data = readdlm(datafile, ',', skipstart=1);

    # Read time variable first to sort chronologically
    year = Int.(data[:,index_year])
    month = Int.(data[:,index_month])
    day = Int.(data[:,index_day])
    dates = Date.(year, month, day)

    if sorting
        @info("Sorting observations chronologically")
        sortingindex = sortperm(dates)

        SampleId = String.(data[sortingindex,index_id])
        lon = Float64.(data[sortingindex,index_longitude])
        lat = Float64.(data[sortingindex,index_latitude])
        Total_Copepods_Large = Float64.(data[sortingindex, index_coplarge])
        Total_Copepods_Small = Float64.(data[sortingindex,index_copsmall])
        dates = dates[sortingindex]

    else
    
        SampleId = String.(data[:,1])
        lon = Float64.(data[:,index_longitude])
        lat = Float64.(data[:,index_latitude])
        Total_Copepods_Large = Float64.(data[:,index_coplarge])
        Total_Copepods_Small = Float64.(data[:,index_copsmall])

    end;
    
    return lon::Vector{Float64}, lat::Vector{Float64}, dates::Vector{Date}, Total_Copepods_Large::Vector{Float64}, Total_Copepods_Small::Vector{Float64}, 
    SampleId::Vector{String}
end

"""
    get_obsid(sampleid)

Create a list of unique identifiers (integers) based on the sample IDs 
(first column of the CSV file)

# Example
```julia-repl
julia> obsid = get_obsid(sampleid)
```
"""
function get_obsid(sampleid::Vector{String})
    SampleID_short = [split(ss, "-")[1] for ss in sampleid];
    SampleID_uni = unique(SampleID_short);
    obsid = ones(Int64, length(sampleid))
    tracknum = ones(Int64, length(SampleID_uni))
    for (iii, values) in enumerate(SampleID_uni)
        sel = findall(SampleID_short .== values)
        obsid[sel] = iii * ones(Int64, length(sel))
        tracknum[iii] = length(sel)
    end
    return obsid::Vector{Int64}, tracknum::Vector{Int64}
end

"""
    get_num_obs(dates, fielddates)

Count the number of observations (made at time `dates`) in each monthly period specified by `fielddates`.
`deltayear` indicate the number of years to be considered before and after the period of interest to be considered.

# Examples
```julia-repl
julia> datestest = [Date(2018, 1, 2), Date(2018, 1, 22), Date(2018, 2, 13), Date(2018, 3, 18), Date(2018, 3, 30)];
julia> fielddates = collect(Date(2018, 1, 15):Dates.Month(1):Date(2018, 3, 15))
julia> get_num_obs(datestest, fielddates)
3-element Vector{Int64}:
 2
 1
 2
```
"""
function get_num_obs(dates::Vector{Date}, fielddates::Vector{Date}, deltayear::Int64=0)

    # allocate
    numobs = zeros(Int64, length(fielddates))

    # Loop on the dates and select observations for the good period
    for (iii, ddd) in enumerate(fielddates)
        themonth = Dates.month(ddd)
        theyear = Dates.year(ddd)
        theyearmin = theyear - deltayear
        theyearmax = theyear + deltayear
        goodtime = findall( (Dates.month.(dates) .== themonth) .&
                            (Dates.year.(dates) .>= theyearmin) .& 
                            (Dates.year.(dates) .<= theyearmax))
        numobs[iii] = length(goodtime)
    end
    
    return numobs

end

"""
    create_random_params()

Return a dictionary storing the values of each parameters to be tested.     
The parameters not present in the dictory will be assigned the default value by `DINCAE`.     

The parameters are taken randomly within an interval that can be modified using the keyword parameters.

# Examples
```julia-repl
julia> parameters = create_random_params()
julia> parameters = create_random_params(epochs_min=1000, epochs_max=5000)
```
"""
function create_random_params(;epochs_min::Int64=500,epochs_max::Int64=5000, 
        learning_rate_min::Float64=1e-4, learning_rate_max::Float64=1e-2,
        laplacian_penalty_min::Float64=1e-6, laplacian_penalty_max::Float64=1e-2)
    
    paramdict = Dict()
    paramdict[:epochs] = rand(epochs_min:10:epochs_max)
    #paramdict["clip_grad"] = [1., 3., 5., 10.]
    #paramdict["ntime_win"] = [1, 3, 5, 7]
    #paramdict["batch_size"] = [8, 16, 32, 64, 128, 256, 512]
    paramdict[:learning_rate] = exp10(rand(Uniform(log10(learning_rate_min), log10(learning_rate_max))))
    #paramdict["jitter_std_pos"] = [(0.1,0.1), (0.5,0.5), (1., 1.), (2., 2.), (5., 5.)]
    #paramdict["enc_nfilter_internal"] = [[32,64,128], [32,64,128,256], [32,64,128,256,512], [32,64,128,256,512,1024]]
    #paramdict["regularization_L2_beta"] = [0.01, 0.1, 1.0, 0.]
    #paramdict["skipconnections"] = [2:5, 3:5, 4:5]
    paramdict[:laplacian_penalty] = exp10(rand(Uniform(log10(laplacian_penalty_min), log10(laplacian_penalty_max))))
    #paramdict[:laplacian_error_penalty] = [0.0001f0, 0.001f0, 0.01f0]

    return paramdict::Dict{Any, Any}
end


"""
    get_index_validation(obsid; percent_validation)

Create two lists of index, one for the observations to keep for the interpolation,
the other for the validation, keeping `percent_validation` percents of the available observations 
(default: 10% are converved)

# Examples
```julia-repl
julia> index4main, index4validation = get_index_validation(obsid)
julia> index4main, index4validation = get_index_validation(obsid, percent_validation=25)
```
"""
function get_index_validation(obsid::Vector{Int64}; percent_validation=10.)
    # Check if percent_validation has a valid value
    if percent_validation > 100 || percent_validation < 0
        @error("percent_validation should take a value between 0 and 100")
    end

    ndata = length(obsid)
    ntracks = length(unique(obsid))
    ntracks_valid = Int(ceil((percent_validation / 100) * ntracks))
    @debug(ntracks_valid);

    # Track numbers to use for the validation
    index4validation = Int[]
    for iii in randperm(MersenneTwister(123), ntracks)[1:ntracks_valid]
        append!(index4validation, findall(obsid .== iii))
    end
    index4main = setdiff(collect(1:ndata), index4validation)
    return index4main::Vector{Int64}, index4validation::Vector{Int64}
end

"""
    write_netcdf_CPR(ncfile, lon, lat, dates, values, varname, obsid, fielddates, valex)

Write the observations in a netCDF file.

# Inputs
- `ncfile` is the path to the netCDF file where the data will be written
- `lon`, `lat`, `dates` are the coordiates 
- `values` are the values of the observations
- `varname` is the name of the variable stored in `values` 
- `obsid` is the list of the observation identifiers
- `fielddates` is a vector of dates that indicates the times at which the field has to be reconstructed

# Example
```julia-repl
julia> write_netcdf_CPR("output.nc", lon, lat, dates, c_fin, "Calanus_finmarchicus", obsid, fielddates, 
    valex=-999.9)
```
"""
function write_netcdf_CPR(ncfile::AbstractString, lon::Vector{Float64}, lat::Vector{Float64}, 
    dates::Vector{Date}, observations::Vector{Float64}, varname::AbstractString, obsid::Vector{Int64}, 
    fielddates::Vector{Date}=Date[]; valex::Float64=-999.)

    ds = NCDatasets.Dataset(ncfile, "c", attrib = OrderedDict(
        "citation"                  => "Pierre Helaouët (Marine Biological Association of the United Kingdom) (2025): CPR Data request - DINCAE - 29/09/2025. The Archive for Marine Species and Habitats Data (DASSH). (Dataset). https://doi.org/10.17031/68da4a97650f1",
        "contact"                   => "Pierre Helaouet <pihe@MBA.ac.uk>",
        "creator_email"             => "ctroupin@uliege.be",
        "creator_name"              => "Charles Troupin",
        "creator_type"              => "University of Liège",
        "creator_url"               => "https://www.uliege.be",
        "date_created"              => Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"),
        "geospatial_lat_max"        => maximum(lat),
        "geospatial_lat_min"        => minimum(lat),
        "geospatial_lat_units"      => "degrees_north",
        "geospatial_lon_max"        => maximum(lon),
        "geospatial_lon_min"        => minimum(lon),
        "geospatial_lon_units"      => "degrees_east",
        "institution"               => "MBA, University of Liège",
        "metadata_character_set"    => "utf8",
        "metadata_language"         => "eng",
        "netcdf_version"            => "netCDF-4 classic model",
        "project"                   => "EMODnet Biology",
        "time_coverage_start"       => Dates.format(minimum(dates), dateformat"yyyy-mm-dd"),
        "time_coverage_end"         => Dates.format(maximum(dates), dateformat"yyyy-mm-dd"),
        "title"                     => "CPR data - $(varname)",
        "method"                    => "CopepodsNN.write_netcdf_CPR",
        "Julia_version"             => string(VERSION),
        "varname"                   => varname,
        "doi"                       => "10.17031/68da4a97650f1"
    ))

    # Dimensions

    ds.dim["obs"] = length(lon)
    ds.dim["timeinstances"] = length(fielddates)
    
    # Declare variables
    nctime = defVar(ds,"dtime", Float64, ("obs",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "actual_range"              => [Dates.format(minimum(dates), dateformat"yyyy-mm-dd"), Dates.format(maximum(dates), dateformat"yyyy-mm-dd")],
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "ioos_category"             => "Time",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))

    nclatitude = defVar(ds,"lat", Float64, ("obs",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Lat",
        "_FillValue"                => Float64(valex),
        "actual_range"              => Float64[minimum(lat), maximum(lat)],
        "axis"                      => "Y",
        "ioos_category"             => "Location",
        "latitude_reference_datum"  => "geographical coordinates, WGS84 projection",
        "long_name"                 => "Latitude",
        "missing_value"             => Float64(valex),
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
        "valid_max"                 => Float64(90.0),
        "valid_min"                 => Float64(-90.0),
    ))

    nclongitude = defVar(ds,"lon", Float64, ("obs",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Lon",
        "_FillValue"                => Float64(valex),
        "actual_range"              => Float64[minimum(lon), maximum(lon)],
        "axis"                      => "X",
        "ioos_category"             => "Location",
        "long_name"                 => "Longitude",
        "longitude_reference_datum" => "geographical coordinates, WGS84 projection",
        "missing_value"             => Float64(valex),
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
        "valid_max"                 => Float64(180.0),
        "valid_min"                 => Float64(-180.0),
    ))

    ncObservations = defVar(ds, varname, observations, ("obs",), attrib = OrderedDict(
        "_FillValue"                => Float64(valex),
        "actual_range"              => Float64[minimum(observations), maximum(observations)],
        "coordinates"               => "time",
        "long_name"                 => "Abundance of $(varname)",
        "missing_value"             => Float64(valex),
        "units"                     => "Individuals per cubic meter",
    ))

    ncid = defVar(ds,"id", Int64, ("obs",), attrib = OrderedDict(
        "long_name"                 => "Identifier",
    ))

    ncnumobs = defVar(ds, "numobs", Int64, ("timeinstances",), attrib = OrderedDict(
		"long_name"                 => "Number of observations in each time period",
        "sample_dimension"          => "obs" 
    ))

    ncdates = defVar(ds, "dates", Float64, ("timeinstances",), attrib = OrderedDict(
		"standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))

    # Define variables
    nclongitude[:] = lon
    nclatitude[:] = lat
    ncid[:] = obsid

    numobs = get_num_obs(dates, fielddates, 0)
    ncnumobs[:] = numobs
    
    if length(fielddates) > 0
        # Convert dates to days since date origin (1900-01-01)
        ncdates[:] = [(dd - Dates.Date(1900, 1, 1)).value for dd in fielddates]
    end
    
    nctime[:] = [(dd - Dates.Date(1900, 1, 1)).value for dd in dates]
    
    close(ds)
end

"""
    read_data_nc(datafile)

Read the coordinates, observations and identifiers from the netCDF file obtained 
with the function `write_netcdf_CPR`. It is assumed that there is one variable per file.
The variable name is defined as a global attribute in the netCDF file.

# Example
```julia-repl
julia> loncpr, latcpr, datescpr, observations, obsid, varname = 
read_data_nc("datafile.nc");
```
"""
function read_data_nc(datafile::AbstractString)
    NCDataset(datafile) do ds
        varname = ds.attrib["varname"]
        lon = coalesce.(NCDatasets.varbyattrib(ds, standard_name="longitude")[1][:], NaN)
        lat = coalesce.(NCDatasets.varbyattrib(ds, standard_name="latitude")[1][:], NaN)
        dates = Dates.Date.(NCDatasets.varbyattrib(ds, standard_name="time")[1][:])
        observations = coalesce.(ds[varname][:], NaN)
        obsid = ds["id"][:]
        
        return lon::Vector{Float64}, lat::Vector{Float64}, dates::Vector{Date}, 
        observations::Vector{Float64}, obsid::Vector{Int64}, varname::AbstractString
    end
end


"""
    get_data_mean_grid(loncpr, latcpr, datescpr, datavalues, longrid, latgrid, timegrid, valex=valex)
    
Create 3D arrays storing the mean, minimal and maximal values of all the observations on a regular grid 
defined by `longrid`, `latgrid` and `timegrid`. 

The positions (`loncpr`, `latcpr`), times (`datescpr`) and values (`datavalues`) of the observations 
are obtained with the function `read_cpr_data_nc`.

The values can then be written in a netCDF file with the function `write_netcdf_grid`.

# Example
```julia-repl
julia> data_mean, data_num = get_data_mean_grid(loncpr, latcpr, datescpr, datavalues, longrid, latgrid, timegrid)
```
"""
function get_data_mean_grid(loncpr::Vector{Float64}, latcpr::Vector{Float64}, datescpr::Vector{Date}, 
    datavalues::Vector{Float64}, longrid::StepRangeLen, latgrid::StepRangeLen, 
    timegrid::Vector{Date}; valex::Float64=-999.)

    # Compute resolution
    Δlon = longrid[2] - longrid[1]
    Δlat = latgrid[2] - latgrid[1]
    
    # Extract years and months from the observations
    monthcpr = Dates.month.(datescpr);
    yearcpr = Dates.year.(datescpr);

    # Allocate matrices (lon, lat, time)
    data_mean = valex .* ones(Union{Missing, Float64}, length(longrid), length(latgrid), length(timegrid))
    data_num = zeros(Int64, length(longrid), length(latgrid), length(timegrid))

    # Loop on time
    Threads.@threads for ttt=1:length(timegrid)#(ttt, ddates) in enumerate(fielddates_monthly)

        ddates = timegrid[ttt]

        # Select good values
        themonth = Dates.month(ddates)
        theyear = Dates.year(ddates)
        gooddates = findall( (monthcpr .== themonth) .& (yearcpr .== theyear))
        @debug(length(gooddates));

        lonsel = loncpr[gooddates]
        latsel = latcpr[gooddates]
        datavalues_dates = datavalues[gooddates]

        # Loop on space
        for (iii, llon) in enumerate(longrid[1:end])
            goodlon = findall( (lonsel .>= llon - Δlon/2) .& (lonsel .< llon + Δlon/2))

            for (jjj, llat) in enumerate(latgrid[1:end])
                goodlat = findall( (latsel .>= llat - Δlat/2) .& (latsel .< llat + Δlat/2))

                datasel = intersect(goodlon, goodlat)
                
                if length(datasel) > 0
                    data_mean[iii,jjj,ttt] = mean(datavalues_dates[datasel])
                    data_num[iii,jjj,ttt] = length(datasel)
                end
            end
        end
    end
    return data_mean::Array{Float64, 3}, data_num::Array{Int64, 3} 
end


"""
    write_netcdf_grid(ncfile, longrid, latgrid, dategrid, data_mean, data_std, data_min, data_max, 
    data_median, data_num, mask, valex)

Write the binned observations in a netCDF file. Inputs are produced by the function `grid_data_mean`.

# Inputs
- `ncfile` is the path to the netCDF file where the data will be written
- `longrid`, `latgrid` are vectors describing the spatial grid 
- `dategrid` is a vector of dates that indicates the times at which the field has to be reconstructed
- `c_fin_mean`, `c_helgo_mean`, ..., `c_helgo_min` are 3D arrays that stored the values on the grid (lon, lat and time)
- `mask` is the land-sea mask
- `valex` is the exclusion value
```
"""
function write_netcdf_grid(ncfile::String, longrid::Vector{Float64}, latgrid::Vector{Float64}, 
    dategrid::Vector{Date}, data_mean::Array{Float64, 3}, data_num::Array{Int64, 3}, mask::BitMatrix; varname::AbstractString="observations",valex::Float64=-999.)

    ds = NCDatasets.Dataset(ncfile, "c", attrib = OrderedDict(
        "citation"                  => "",
        "comment"                   => "Data binned on 2D grid",
        "contact"                   => "Charles Troupin <ctroupin@uliege.be>",
        "creator_email"             => "ctroupin@uliege.be",
        "creator_name"              => "Charles Troupin",
        "creator_type"              => "University of Liège",
        "creator_url"               => "https://www.uliege.be",
        "date_created"              => Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"),
        "geospatial_lat_max"        => maximum(latgrid),
        "geospatial_lat_min"        => minimum(latgrid),
        "geospatial_lat_units"      => "degrees_north",
        "geospatial_lon_max"        => maximum(longrid),
        "geospatial_lon_min"        => minimum(longrid),
        "geospatial_lon_units"      => "degrees_east",
        "institution"               => "MBA, University of Liège",
        "license"                   => " ",
        "metadata_character_set"    => "utf8",
        "metadata_language"         => "eng",
        "netcdf_version"            => "netCDF-4 classic model",
        "time_coverage_start"       => Dates.format(minimum(dategrid), dateformat"yyyy-mm-dd"),
        "time_coverage_end"         => Dates.format(maximum(dategrid), dateformat"yyyy-mm-dd"),
        "title"                     => "Binned values of observations of $(varname)",
        "method"                    => "write_netcdf_grid",
        "Julia_version"             => juliaversion,
        "varname"                   => varname
    ))

    # Dimensions

    ds.dim["lon"] = length(longrid)
    ds.dim["lat"] = length(latgrid)
    ds.dim["time"] = length(dategrid)

    # Declare variables
    nctime = defVar(ds, "time", Float64, ("time",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "actual_range"              => [Dates.format(minimum(dategrid), dateformat"yyyy-mm-dd"), 
    Dates.format(maximum(dategrid), dateformat"yyyy-mm-dd")],
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "ioos_category"             => "Time",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))

    defVar(ds,"lat", latgrid, ("lat",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Lat",
        "_FillValue"                => Float64(valex),
        "actual_range"              => Float64[minimum(latgrid), maximum(latgrid)],
        "axis"                      => "Y",
        "ioos_category"             => "Location",
        "latitude_reference_datum"  => "geographical coordinates, WGS84 projection",
        "long_name"                 => "Latitude",
        "missing_value"             => Float64(valex),
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
        "valid_max"                 => Float64(90.0),
        "valid_min"                 => Float64(-90.0),
    ))

    defVar(ds,"lon", longrid, ("lon",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Lon",
        "_FillValue"                => Float64(valex),
        "actual_range"              => Float64[minimum(longrid), maximum(longrid)],
        "axis"                      => "X",
        "ioos_category"             => "Location",
        "long_name"                 => "Longitude",
        "longitude_reference_datum" => "geographical coordinates, WGS84 projection",
        "missing_value"             => Float64(valex),
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
        "valid_max"                 => Float64(180.0),
        "valid_min"                 => Float64(-180.0),
    ))

    defVar(ds,"$(varname)_mean", data_mean, ("lon", "lat", "time"), attrib = OrderedDict(
        "_FillValue"                => valex,
        "actual_range"              => [minimum(data_mean), maximum(data_mean)],
        "coordinates"               => "time, lat, lon",
        "long_name"                 => "Mean abundance of $(varname)",
        "missing_value"             => valex,
        "units"                     => "Individuals per cubic meter",
    ))

    defVar(ds,"$(varname)_num", data_num, ("lon", "lat", "time"), attrib = OrderedDict(
        "actual_range"              => [0, maximum(data_num)],
        "coordinates"               => "time, lat, lon",
        "long_name"                 => "Number of $(varname) observations per grid cell",
        "units"                     => "Individuals per cubic meter",
    ))

    ncmask = defVar(ds,"mask", Int8, ("lon", "lat"), attrib = OrderedDict(
        "coordinates"               => "lat, lon",
        "long_name"                 => "mask (sea=1, land=0)",
        "units"                     => "1s"
    ))

    # Define time and mask
    nctime[:] = [(dd - Dates.Date(1900, 1, 1)).value for dd in dategrid]
    ncmask[:] = mask

    close(ds)
    return nothing
end

"""
    compute_data_num(loncpr, latcpr, datescpr, longrid, latgrid, timegrid, valex=valex)
    
Compute the number of observations in each grid cell defined by `longrid`, `latgrid` and `timegrid`. 

The values can then be written in a netCDF file with the function `write_netcdf_datanum_grid`.
"""
function compute_data_num(loncpr::Vector{Float64}, latcpr::Vector{Float64}, datescpr::Vector{Date}, 
    longrid::StepRangeLen, latgrid::StepRangeLen, timegrid::Vector{Date})

    Δlon = longrid[2] - longrid[1]
    Δlat = latgrid[2] - latgrid[1]
    
    # Extract years and months
    monthcpr = Dates.month.(datescpr);
    yearcpr = Dates.year.(datescpr);

    # Allocate matrices
    data_num = zeros(Int64, length(longrid), length(latgrid), length(timegrid))

    # Loop on time
    Threads.@threads for ttt=1:length(timegrid)#(ttt, ddates) in enumerate(fielddates_monthly)

        ddates = timegrid[ttt]
        # Select good values
        themonth = Dates.month(ddates)
        theyear = Dates.year(ddates)
        gooddates = findall( (monthcpr .== themonth) .& (yearcpr .== theyear))
        @debug(length(gooddates));

        lonsel = loncpr[gooddates]
        latsel = latcpr[gooddates]

        # Loop on space
        for (iii, llon) in enumerate(longrid[1:end])
            goodlon = findall( (lonsel .>= llon - Δlon/2) .& (lonsel .< llon + Δlon/2))

            for (jjj, llat) in enumerate(latgrid[1:end])
                goodlat = findall( (latsel .>= llat - Δlat/2) .& (latsel .< llat + Δlat/2))

                datasel = intersect(goodlon, goodlat)
                
                if length(datasel) > 0
                    data_num[iii,jjj,ttt] = length(datasel)
                end
            end
        end
    end
    return data_num::Array{Int64, 3}

end


"""
    write_netcdf_datanum_grid(ncfile, longrid, latgrid, dategrid, data_num, mask, valex)

Write the number of observations per grid cell (space and time) in a netCDF file. 
Inputs are produced by the function `compute_data_num`.

# Inputs
- `ncfile` is the path to the netCDF file where the data will be written
- `longrid`, `latgrid` are vectors describing the spatial grid 
- `dategrid` is a vector of dates that indicates the times at which the field has to be reconstructed
- `data_num` is a 3D array that stores the number of observations on the grid (lon, lat and time)
- `mask` is the land-sea mask
- `valex` is the exclusion value
```

# Example
```julia-repl
julia> write_netcdf_datanum_grid("Calanus_Finmarchicus_grid.nc", longrid, latgrid, 
dategrid, data_num, mask)
```
"""
function write_netcdf_datanum_grid(ncfile::String, longrid::Vector{Float64}, latgrid::Vector{Float64}, 
    dategrid::Vector{Date}, data_num::Array{Int64, 3}; valex::Float64=-999.)

    ds = NCDatasets.Dataset(ncfile, "c", attrib = OrderedDict(
        "citation"                  => "",
        "comment"                   => "Data binned on 2D grid",
        "contact"                   => "Charles Troupin <ctroupin@uliege.be>",
        "creator_email"             => "ctroupin@uliege.be",
        "creator_name"              => "Charles Troupin",
        "creator_type"              => "University of Liège",
        "creator_url"               => "https://www.uliege.be",
        "date_created"              => Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"),
        "geospatial_lat_max"        => maximum(latgrid),
        "geospatial_lat_min"        => minimum(latgrid),
        "geospatial_lat_units"      => "degrees_north",
        "geospatial_lon_max"        => maximum(longrid),
        "geospatial_lon_min"        => minimum(longrid),
        "geospatial_lon_units"      => "degrees_east",
        "institution"               => "MBA, University of Liège",
        "keywords"                  => "data, earth, Earth Science > Oceans, global, latitude, longitude, ocean, oceans, quality, research, science, sea, seawater, system, time, water",
        "keywords_vocabulary"       => "GCMD Science Keywords",
        "license"                   => " ",
        "metadata_character_set"    => "utf8",
        "metadata_language"         => "eng",
        "netcdf_version"            => "netCDF-4 classic model",
        "time_coverage_start"       => Dates.format(minimum(dategrid), dateformat"yyyy-mm-dd"),
        "time_coverage_end"         => Dates.format(maximum(dategrid), dateformat"yyyy-mm-dd"),
        "title"                     => "Number of observations per grid cell",
        "method"                    => "write_netcdf_datanum_grid",
        "Julia_version"             => juliaversion
    ))

    # Dimensions

    ds.dim["lon"] = length(longrid)
    ds.dim["lat"] = length(latgrid)
    ds.dim["time"] = length(dategrid)

    # Declare variables
    nctime = defVar(ds, "time", Float64, ("time",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "actual_range"              => [Dates.format(minimum(dategrid), dateformat"yyyy-mm-dd"), 
    Dates.format(maximum(dategrid), dateformat"yyyy-mm-dd")],
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "ioos_category"             => "Time",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))

    defVar(ds,"lat", latgrid, ("lat",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Lat",
        "_FillValue"                => Float64(valex),
        "actual_range"              => Float64[minimum(latgrid), maximum(latgrid)],
        "axis"                      => "Y",
        "ioos_category"             => "Location",
        "latitude_reference_datum"  => "geographical coordinates, WGS84 projection",
        "long_name"                 => "Latitude",
        "missing_value"             => Float64(valex),
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
        "valid_max"                 => Float64(90.0),
        "valid_min"                 => Float64(-90.0),
    ))

    defVar(ds, "lon", longrid, ("lon",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Lon",
        "_FillValue"                => Float64(valex),
        "actual_range"              => Float64[minimum(longrid), maximum(longrid)],
        "axis"                      => "X",
        "ioos_category"             => "Location",
        "long_name"                 => "Longitude",
        "longitude_reference_datum" => "geographical coordinates, WGS84 projection",
        "missing_value"             => Float64(valex),
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
        "valid_max"                 => Float64(180.0),
        "valid_min"                 => Float64(-180.0),
    ))

    defVar(ds, "data_num", data_num, ("lon", "lat", "time"), attrib = OrderedDict(
        "actual_range"              => [0, maximum(data_num)],
        "coordinates"               => "time, lat, lon",
        "long_name"                 => "Number of Small copepods observations per grid cell",
        "units"                     => "Individuals per cubic meter",
    ))


    # Define variables
    nctime[:] = [(dd - Dates.Date(1900, 1, 1)).value for dd in dategrid]

    close(ds)
    return nothing
end

"""
    create_geoaxis(fig)

Create a GeoAxis for the current figure `fig`
"""
function create_geoaxis(fig::Figure; proj="+proj=merc", title="", ii=1, jj=1)    
    ga = GeoAxis(fig[ii, jj],
    dest = proj,
    xgridstyle = :dashdot,
    xgridcolor = :grey,
    ygridstyle = :dashdot,
    ygridcolor = :grey,
    xticks = collect(-110:20:30),
    yticks = collect(10.0:10:80),
    title = title
    )

    return ga::GeoAxis
end

end # end of module 
