"""
# Module VectorMaths

Data Analysis and Manipulation tools containing functions:

- `DFvector_mean` to return a DataFrame with common x data, adjusted y data,
  the mean and sum of the y data of two DataFrames
- `multivector_mean` to return a DataFrame with the common x data, adjusted y data,
  the mean and sum of the y data of any number of individaul DataFrames with an
  x and a y column either for their common x data range or over the whole x data
  range
"""
module VectorMaths


##################
###  PREAMBLE  ###
##################

# Load packages
using DataFrames
using Dierckx
using Statistics

export vector_maths


##########################
###  PUBLIC FUNCTIONS  ###
##########################


"""
    multivector_mean(vectors::Vector{DataFrames.DataFrame};
                     output::String="offset")

From a vector of DataFrames (`vectors`) with x values in the first column and
y values in the second column, return a single DataFrame with unified x data in
the first column `x`, the ydata in the following columns `y1` to `yn`, the mean
in column `mean` and the sum in column `sum`. The x data in each DataFrame must
be strictly monotonic. Moreover, there must be a range, where the x data of all
vectors overlap.

There are two options for the keyword argument `datarange`: `common` and `all`.
If `common` is chosen, data is only returned for the overlapping range of all
vectors. If the data have different step sizes, all x values are collected and
cubic splines are used to interpolate missing data.

If `all` is chosen, all vectors are filled with `NaN`s, if their range is shorter
than the maximum range of all vectors. `NaN`s are ignored for the calculation of
the `sum` and the `mean`.

Furthermore, for vectors of different lengths, the mean is either multiplied by
a scaling factor or an offset is added depending on the choice in the second
keyword argment `output`: `"offset"` or `"scale"`.
Scaling factors or offsets are applied before every delayed start of a vector and
after each premature end of a vector, respectively, to assure a contineous smooth
mean even if an outlier has ended. Offsets and means are calcualted as

    offset = Ø(y values of all vectors)-Ø(y values of remaining vectors)
    SF = Ø(y values of all vectors)/Ø(y values of remaining vectors)

This means, that the ends of each shorter vector will influence the mean outside
the range of this vector. Scaling factors should not be used, if the y data spans
a large range (over several orders of magnitudes) as otherwise a small offset at
the lower end of the data could lead to a large scaled offset at the upper end of
the y data. Scaling factors should be used, if you want to assure that the sign
in your y data does not change and the y data is in a reasonably small range.
"""
function vector_maths(vectors::DataFrames.DataFrame...;
  output::String="offset",
  xcols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[],
  ycols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[])

  # Define x and y data in each DataFrame
  xdata, ydata, xcols, ycols = compile_data(vectors, xcols, ycols)

  ### Check input data
  fail = test_monotonicity(vectors, xcols)
  if fail  return nothing  end

  ### Find common range of all vectors
  # Find start/end x value of each vector (DataFrame)
  mins = Float64[]; maxs = Float64[]
  for x in xdata
    push!(mins,x[1])
    push!(maxs,x[end])
  end
  # Find boundaries of common range in all datasets
  lowest_common  = maximum(mins)
  largest_common = minimum(maxs)

  ### Get all data points in each DataFrame within the common range
  common_xdata = get_xdata(xdata, lowest_common, largest_common)

  # Get cubic splines of all vectors
  spl = Dierckx.Spline1D[]
  for i = 1:length(xdata)
    push!(spl,Spline1D(xdata[i],ydata[i]))
  end


  # Assign y data from each DataFrame, use cubic spines for missing data
  common_ydata = get_ydata(common_xdata,xdata,ydata,mins,maxs,spl)

  # Generate output DataFrame with common x data, individual y data columns,
  # the mean and sum
  if output == "common"
    # Compile output DataFrame
    dfr = DataFrame(x = common_xdata)
    for i = 1:length(ydata)
      dfr[Symbol("y$i")] = common_ydata[:,i]
    end
    dfr[:mean] = vec(mean(common_ydata, dims=2))
    dfr[:sum]  = vec(sum(common_ydata, dims=2))
  else
    # Find all x and y data outside common range
    # Assign NaN's to missing data
    lower_xdata = get_xdata(xdata, -Inf, lowest_common, bounds="exclude")
    upper_xdata = get_xdata(xdata, largest_common, Inf, bounds="exclude")

    # Find indices in arrays, where some data columns end in the lower or upper data
    mi, Mi = find_boundaries(mins,maxs,lower_xdata,upper_xdata)

    lower_ydata = get_ydata(lower_xdata,xdata,ydata,mins,maxs,spl)
    upper_ydata = get_ydata(upper_xdata,xdata,ydata,mins,maxs,spl)
    return mi, Mi
  end
  #=
    SF = get_ScalingFactor(output, lower_xdata, common_xdata, upper_xdata,
                           lower_ydata, common_ydata, upper_ydata, mi, Mi)

    # Get mean and sum
    lower_mean, lower_sum = calc_MeanSum(lower_ydata,SF[1],output)
    upper_mean, upper_sum = calc_MeanSum(upper_ydata,SF[2],output)

    # Compile output DataFrame
    dfr = DataFrame(x = vcat(lower_xdata, common_xdata, upper_xdata))
    for i = 1:length(vectors)
      dfr[Symbol("y$i")] = vcat(lower_ydata[:,i], common_ydata[:,i], upper_ydata[:,i])
    end
    dfr[:mean] = vcat(lower_mean, vec(mean(common_ydata,2)), upper_mean)
    dfr[:sum]  = vcat(lower_sum, vec(sum(common_ydata,2)), upper_sum)
  end
  =#

  return dfr
end #function vector_maths


###########################
###  PRIVATE FUNCTIONS  ###
###########################


"""
    compile_data(vectors, xcols, ycols)

For each DataFrame in `vectors`, find the index of the x columns defined in `xcols`
or use first column by default and the index of the y columns defined in `ycols`
or use second column as default and return vectors with x and y data (as vectors)
and the revised x and y column indices.
"""
function compile_data(vectors, xcols, ycols)
  # Define x and y data in each DataFrame
  xcols = check_columns(vectors, xcols, 'x')
  ycols = check_columns(vectors, ycols, 'y')
  xdata = []; ydata = []

  # Loop over vectors
  for i = 1:length(vectors)
    # Filter mis-assigned column data
    if length(xcols[i]) > 1 && length(xcols[i]) ≠ length(ycols[i])
      println("\33[95mError! Different number of x and y columns definitions.")
      println("$i. DataFrame ignored.\33[0m")
      continue
    end
    # Retrieve x and y data from DataFrames
    if xcols[i] isa Vector
      for j in xcols[i]  push!(xdata, vectors[i][j])  end
      ydata = assign_ydata(ycols[i], ydata, vectors[i])
    elseif ycols[i] isa Vector
      for j = 1:length(ycols[i])  push!(xdata, vectors[i][xcols[i]])  end
      ydata = assign_ydata(ycols[i], ydata, vectors[i])
    else
      push!(xdata, vectors[i][xcols[i]])
      ydata = assign_ydata(ycols[i], ydata, vectors[i])
    end
  end #loop over vectors

  return xdata, ydata, xcols, ycols
end #function compile_data


"""
    assign_ydata(ycols, ydata, vector)

From the column index `ycols` in the current DataFrame `vector`, save the respective
column as vector to `ydata` and return the appended array.
"""
function assign_ydata(ycols, ydata, vector)
  if ycols isa Vector
    [push!(ydata, vector[j]) for j in ycols]
  else
    push!(ydata, vector[ycols])
  end
  return ydata
end


"""
    check_columns(vectors, cols, coltype::Char)

For every `DataFrame` in `vectors`, check that columns are assigned and assign
first column as x data and second column as ydata for missing assignments or
delete definitions in `cols` in excess of `DataFrames` in `vectors`.

`coltype` defines the data type (`'x'` or `'y'` data) in the columns of the `DataFrames`
in `vectors`.
"""
function check_columns(vectors, cols, coltype::Char)
  if isempty(cols) && coltype == 'x'
    cols = ones(Int64, length(vectors))
  elseif isempty(cols) && coltype == 'y'
    cols = 2 .* ones(Int64, length(vectors))
  elseif length(cols) > length(vectors)
    println("\033[95mWarning! More $coltype columns defined than vectors available.\33[0m")
    println("Last $(length(cols) - length(vectors)) columns in `$(coltype)cols` ignored.")
    deleteat!(cols, 1 + length(vectors):length(cols))
  elseif length(cols) < length(vectors)
    println("\033[95mWarning! Less x columns defined than vectors available.\33[0m")
    println("First column assumed to hold x data in last ",
      "$(length(vectors) - length(cols)) DataFrames.")
    if (typeof(cols[1]) == Symbol || typeof(cols[1]) == Vector{Symbol}) &&
      coltype == 'x'
      [push!(cols, names(vectors[i])[1]) for i = 1 + length(cols):length(vectors)]
    elseif (typeof(cols[1]) == Symbol || typeof(cols[1]) == Vector{Symbol}) &&
      coltype == 'y'
      [push!(cols, names(vectors[i])[2]) for i = 1 + length(cols):length(vectors)]
    elseif coltype == 'x'
      [push!(cols, 1) for i = 1:length(vectors) - length(cols)]
    elseif coltype == 'y'
      [push!(cols, 2) for i = 1:length(vectors) - length(cols)]
    end
  end #function check_column_defs

  return cols
end #function check_columns


"""
    test_monotonicity(vectors, x)

Test the the first column in each DataFrame in a vector of DataFrames (`vectors`)
is strictly monotonic. Issue a warning and stop the script, if not. Use the vector
with indices `x` to find the x column in each DataFrame.
"""
function test_monotonicity(vectors, x)
  fail = false
  for v = 1:length(vectors)
    if vectors[v][x[v]] != sort(unique(vectors[v][x[v]]))
      println("\033[95mError! X data in dataframe $v not strictly monotonic.")
      println("Script stopped.\033[0m")
      fail = true
      break
    end
  end

  return fail
end #function test_monotonicity

"""
    get_xdata(xdata,lb,ub;bounds::String="include")

From the `xdata`, extract the  common x data of all DataFrames
within the given bounds `lb` and `ub` and return it. The inclusion of the
boundaries can be set with the keyword argument `bounds`
(`"include"` (default)/`"exclude"`).
"""
function get_xdata(xdata,lb,ub;bounds::String="include")
  # Retrieve indices of all x datapoints within the specified boundaries
  xrange = []
  for x in xdata
    if bounds == "include"
      push!(xrange,findall(lb.≤x.≤ub))
    else bounds == "exclude"
      push!(xrange,findall(lb.<x.<ub))
    end
  end

  # Collect all possible datapoints from each vector and unify the data
  commonx = Float64[]
  for (i, x) in enumerate(xdata)
    commonx = vcat(commonx,x[xrange[i]])
  end
  commonx = unique(sort(commonx))
  # unifiedx = unique(sort(vcat(xdata...)))

  # Return a unified x data array
  return commonx#, unifiedx
end #function get_xdata


"""
    get_ydata(xdata,ydata,mins,maxs,spl)

From the `xdata` in the boundaries `mins`/`maxs` and the respective `ydata`
as well as splines `spl` of each dataset, return a matrix with unified
ydata using the selected xdata range and supplementing missing ydata
(due to different step sizes in the x datasets) within the common range with
estimates from the splines and fill missing data outside the common range with `NaN`s.
"""
function get_ydata(xrange,xdata,ydata,mins,maxs,spl)
  # Initialise output matrix
  all_y = Matrix{Float64}(undef,length(xrange),0)
  # Loop over all DataFrames
  for i = 1:length(xdata)
    # Initialise y data of current DataFrame
    yd = Float64[]
    # Loop over unified x data
    for x in xrange
      # Try to add value for current x value, if y is missing
      # interpolate with cubic spline, if inside the range of the current vector
      # or fill with NaN outside its range
      try push!(yd,ydata[i][findfirst(xdata.==x)])
      catch
        if mins[i] ≤ x ≤ maxs[i]
          push!(yd,spl[i](x))
        else
          push!(yd,NaN)
        end
      end
    end
    # Add a column with completed y data of current DataFrame to the output matrix
    all_y = hcat(all_y,yd)
  end

  # Return unified y data
  return all_y
end #function get_ydata


"""
get_ScalingFactor(output, lower_xdata, common_xdata, upper_xdata,
                  lower_ydata, common_ydata, upper_ydata, mi, Mi)

From the unified x and y data below, inside, and above the common range and
the indices of ending vectors outside the common range in the unified x data
(`mi` and `Mi`), calculate and return either offsets or scaling factors (as
defined in `output` with the keywords `"offset"` or `"scale"`) to allow
a smooth contineous mean at the delayed start or premature ends of vectors.
"""
function get_ScalingFactor(output, lower_xdata, common_xdata, upper_xdata,
                           lower_ydata, common_ydata, upper_ydata, mi, Mi)

  ### Get scaling factors for data below the common x range
  flow = ones(lower_xdata); fhigh = ones(upper_xdata)
  if isempty(mi) && length(lower_ydata) > 0
    # Calculate scaling factors, if no premature ending vectors are found
    # outside the common range
    if output == "offset"
      flow[1:end] = mean(common_ydata[1,:]) -
                    mean(common_ydata[1,:][!isnan.(lower_ydata[end,:])])
    else
      flow[1:end] = mean(common_ydata[1,:])/
                    mean(common_ydata[1,:][!isnan.(lower_ydata[end,:])])
    end
  elseif length(lower_ydata) > 0
    # Calculate scaling factors, if premature ending vectors outside the
    # common range are found
    y = lower_ydata[mi[1],:]
    # Scaling factors from common range to first ending vectors
    if output == "offset"
      flow[1:mi[1]-1] = mean(y[!isnan.(y)])-mean(y[!isnan.(lower_ydata[mi[1]-1,:])])
    else
      flow[1:mi[1]-1] = mean(y[!isnan.(y)])/mean(y[!isnan.(lower_ydata[mi[1]-1,:])])
    end
    # Scaling vectors for further ending vectors
    for i = 2:length(mi)
      y = lower_ydata[mi[i],:]
      if output == "offset"
        flow[mi[i-1]:mi[i]-1] = mean(y[!isnan.(y)]) -
                                mean(y[!isnan.(lower_ydata[mi[i]-1,:])])
      else
        flow[mi[i-1]:mi[i]-1] = mean(y[!isnan.(y)])/
                                mean(y[!isnan.(lower_ydata[mi[i]-1,:])])
      end
    end
    # Scaling factors for last prematurely ending vector to end of data range
    if output == "offset"
      flow[mi[end]:end] = mean(common_ydata[1,:]) -
                          mean(common_ydata[1,:][!isnan.(lower_ydata[end,:])])
    else
      flow[mi[end]:end] = mean(common_ydata[1,:])/
                          mean(common_ydata[1,:][!isnan.(lower_ydata[end,:])])
    end
  end

  ### Get scaling factors for data above the common x range
  if isempty(Mi) && length(upper_ydata) > 0
    # Calculate scaling factors, if no premature ending vectors are found
    # outside the common range
    if output == "offset"
      fhigh[1:end] = mean(common_ydata[end,:]) -
                     mean(common_ydata[end,:][!isnan.(upper_ydata[1,:])])
    else
      fhigh[1:end] = mean(common_ydata[end,:])/
                     mean(common_ydata[end,:][!isnan.(upper_ydata[1,:])])
    end
  elseif length(upper_ydata) > 0
    # Calculate scaling factors, if premature ending vectors outside the
    # common range are found
    # Scaling factors from common range to first ending vectors
    if output == "offset"
      fhigh[1:Mi[1]-1] = mean(common_ydata[end,:]) -
                         mean(common_ydata[end,:][!isnan.(upper_ydata[1,:])])
    else
      fhigh[1:Mi[1]-1] = mean(common_ydata[end,:])/
                         mean(common_ydata[end,:][!isnan.(upper_ydata[1,:])])
    end
    # Scaling vectors for further ending vectors
    for i = 2:length(Mi)
      y = upper_ydata[Mi[i],:]
      if output == "offset"
        fhigh[Mi[i-1]:Mi[i]-1] = mean(y[!isnan.(y)]) -
                                 mean(y[!isnan.(upper_ydata[Mi[i]+1,:])])
      else
        fhigh[Mi[i-1]:Mi[i]-1] = mean(y[!isnan.(y)])/
                                 mean(y[!isnan.(upper_ydata[Mi[i]+1,:])])
      end
    end
    # Scaling vectors for further ending vectors
    y = upper_ydata[Mi[end],:]
    if output == "offset"
      fhigh[Mi[end]:end] = mean(y[!isnan.(y)]) -
                           mean(y[!isnan.(upper_ydata[Mi[end]+1,:])])
    else
      fhigh[Mi[end]:end] = mean(y[!isnan.(y)])/
                           mean(y[!isnan.(upper_ydata[Mi[end]+1,:])])
    end
  end

  # Return scaling factors for mean below and above the common range
  return flow, fhigh
end #function get_ScalingFactor


"""
    find_boundaries(mins,maxs,lower_xdata,upper_xdata)

From the minima (`mins`) and maxima (`maxs`) in the x data of each DataFrame,
and the `lower_xdata` and `upper_xdata` below and above the common range, return
all indices in `lower_xdata` and `upper_xdata`, where vectors are prematurely
ending outside the common range.
"""
function find_boundaries(mins,maxs,lower_xdata,upper_xdata)
  # Initialise
  mi = Int64[]; Mi = Int64[]
  # Loop over extrema
  for m = 1:length(mins)
    # Save vectors starting later the absolute minimum,
    # but before common range
    if all(float(mins[m]).!=extrema(mins))
      push!(mi,findfirst(lower_xdata.==mins[m]))
    end
    # Save prematurely ending vectors outside common range
    if all(float(maxs[m]).!=extrema(maxs))
      push!(Mi,findfirst(upper_xdata.==maxs[m]))
    end
  end
  # Unify data for more than one ending vector at the same point
  mi = sort(unique(mi)); Mi = sort(unique(Mi))

  # Return indices
  return mi, Mi
end


"""
    calc_MeanSum(Ydata,SF,output)

From the unified `Ydata` of each DataFrame and the corresponding offsets or
scaling factors `SF` for the mean outside the common range (as set by
`output` with the keywords `"offset"` or `"scale"`), calculate the `mean`
and `sum` of all `Ydata`.
"""
function calc_MeanSum(Ydata,SF,output)
  Ymean = Float64[]; Ysum = Float64[]
  for i = 1:length(SF)
    if output == "offset"
      push!(Ymean,mean(Ydata[i,:][!isnan.(Ydata[i,:])])+SF[i])
    else
      push!(Ymean,mean(Ydata[i,:][!isnan.(Ydata[i,:])])⋅SF[i])
    end
    push!(Ysum,sum(Ydata[i,:][!isnan.(Ydata[i,:])]))
  end

  return Ymean, Ysum
end #function calc_MeanSum

end #module VectorMaths
