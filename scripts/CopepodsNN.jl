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

    end

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

end # end of module 
