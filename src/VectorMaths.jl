"""
# Module VectorMaths

Data Analysis and Manipulation tools currently with function `vector_maths` to
calculate the mean and sum of selected vectors stored in DataFrames.
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
    vector_maths(vectors::DataFrames.DataFrame...; output::String="offset",
      xcols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[],
      ycols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[])

Select any number of DataFrames stored in vararg `vectors` and, from the x and y
column data in these DataFrames, calculate the mean and sum of all vectors.
A DataFrame is return with common `x` data as first column, ydata (`y1` ... `yn`)
in the following columns, where missing data outside the common range is filled
with `NaN`, and a column with the `mean` data and `sum`.

By default, the first column in each DataFrame is treated as x data and the second
column as y data. This can be changed by the kwargs `xcols` and `ycols`, which
specify the x and y columns in a DataFrame as follows:

`xcols`/`ycols` must be vectors holding either the position of the column in the
DataFrame as integer or the column name as `Symbol`. If you want to select more
than one column from one DataFrame, use vectors of integers or symbols as the
respective entry in `xcols`/`ycols`. You can either specify the same amount of
x and y columns in paired `xcols`/`ycols` entries or, if the y data in a DataFrame
shares the same x data, you can specify the x data with one integer or symbol and
the respective y data columns as vectors of integers or symbols.

To calculate the mean and sum, the following options exist, which can be specified
by the kwarg `output`:
- `"common"`: If vectors have different lengths, only for the overlapping range
  a mean and sum is calculated.
- `"offset"`: Mean and sum outside the common range are calculated by using an offset
  `offset + data` to avoid unneccessary steps in the output.
- `"scale"` (**default**): Mean and sum outside the common range are calculated
  by using a scaling factor `SF · data` to avoid unneccessary steps in the output.

    offset = Ø(y values of all vectors)-Ø(y values of remaining vectors)
    SF = Ø(y values of all vectors)/Ø(y values of remaining vectors)

This means, that the ends of each shorter vector will influence the mean outside
the range of this vector. Scaling factors should not be used, if the y data spans
a large range (over several orders of magnitudes) as otherwise a small offset at
the lower end of the data could lead to a large scaled offset at the upper end of
the y data. Scaling factors should be used, if you want to assure that the sign
in your y data does not change and the y data is in a reasonably small range.
"""
function vector_maths(vectors::DataFrames.DataFrame...; output::String="offset",
  xcols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[],
  ycols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[])

  # Retrieve x and y data in each DataFrame
  xdata, ydata = compile_data(vectors, xcols, ycols)

  ### Check for strictly monotonic rising x data
  test_monotonicity(xdata)

  ### Combine x data
  unifiedx, imin, imax = unify_xdata(xdata)

  # Get cubic splines of all vectors
  spl = dataspline(xdata, ydata)

  ### Assign y data from each DataFrame
  #   Use cubic spines for missing data in overlapping range and NaN's outside
  unifiedy = get_ydata(unifiedx,xdata,ydata,spl)


  if output == "common"
    # Generate output DataFrame with common x data, individual y data columns,
    # the mean and sum

    # Compile output DataFrame
    dfr = DataFrame(x = unifiedx[imin[end]:imax[1]])
    for i = 1:length(ydata)
      dfr[Symbol("y$i")] = unifiedy[imin[end]:imax[1],i]
    end
    # Add mean and sum of all vectors
    dfr[:mean] = vec(mean(unifiedy[imin[end]:imax[1],:], dims=2))
    dfr[:sum]  = vec(sum(unifiedy[imin[end]:imax[1],:], dims=2))
  else
    # Generate mean and sum of all data

    # Calculate scaling factors/offsets for data outside the common range
    SF = get_ScalingFactor(unifiedx, unifiedy, imin, imax, output)
    # Compile output DataFrame
    dfr = DataFrame(x = unifiedx)
    for i = 1:length(ydata)
      dfr[Symbol("y$i")] = unifiedy[:,i]
    end
    # Calculate mean and sum using the scaling factors/offsets
    if output == "scale"
      dfr[:mean] = SF.*[mean(filter(!isnan, unifiedy[i,:])) for i = 1:length(unifiedx)]
      dfr[:sum]  = SF.*[sum(filter(!isnan, unifiedy[i,:])) for i = 1:length(unifiedx)]
    elseif output == "offset"
      dfr[:mean] = SF.+[mean(filter(!isnan, unifiedy[i,:])) for i = 1:length(unifiedx)]
      dfr[:sum]  = SF.+[sum(filter(!isnan, unifiedy[i,:])) for i = 1:length(unifiedx)]
    end
  end
  # Return output DataFrame
  return dfr
end #function vector_maths


###########################
###  PRIVATE FUNCTIONS  ###
###########################


"""
    compile_data(vectors, xcols, ycols)

For each DataFrame in `vectors`, find the index of the x columns defined in `xcols`
or use first columns by default and the index of the y columns defined in `ycols`
or use second columns as default and return vectors with x and y data vectors.
"""
function compile_data(vectors, xcols, ycols)
  # Define x and y data in each DataFrame
  xcols = check_columns(vectors, xcols, 'x')
  ycols = check_columns(vectors, ycols, 'y')
  xdata = []; ydata = []

  # Loop over vectors
  for i = 1:length(vectors)
    # Filter mis-assigned column data
    if xcols[i] isa Vector && length(xcols[i]) > 1 && length(xcols[i]) ≠ length(ycols[i])
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

  return xdata, ydata
end #function compile_data


"""
    check_columns(vectors, cols, coltype::Char)

For every `DataFrame` in `vectors`, check that columns are assigned or assign
first column as x data and second column as ydata by default and
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
end #function assign_ydata


"""
    test_monotonicity(vectors, x)

Test that the `x` column in each DataFrame in a vector of DataFrames (`vectors`)
is strictly monotonic. Issue an error, if not. Use the vector with indices `x`
to find the x column in each DataFrame.
"""
function test_monotonicity(xdata)
  fail = false
  for (i, x) in enumerate(xdata)
    if x ≠ sort(unique(x))
      error("X data in dataframe $i not strictly monotonic.\n$x")
    end
  end
end #function test_monotonicity


"""
    unify_xdata(xdata)

Combine all `xdata` of every vector and return as vector with ascending `xdata`,
where every is listed once.
"""
function unify_xdata(xdata)
  # Initialise
  mins = Float64[]; maxs = Float64[]
  # Loop over all vectors and get min/max xdata
  for x in xdata
    push!(mins,x[1])
    push!(maxs,x[end])
  end
  # Sort and remove duplicates
  mins = sort(unique(mins))
  maxs = sort(unique(maxs))

  # Check for common range
  if maximum(mins) > minimum(maxs)
    error("No common x data range.")
  end
  # Return unified x data and indices for lower/upper bounds
  unifiedx = unique(sort(vcat(xdata...)))

  return unifiedx, [findfirst(unifiedx.==i) for i in mins],
    [findfirst(unifiedx.==i) for i in maxs]
end #function unify_xdata


"""
    dataspline(xdata, ydata)

Derive a cubic spline for the vectors with `xdata` and `ydata` and return it.
"""
function dataspline(xdata, ydata)
  spl = Dierckx.Spline1D[]
  for i = 1:length(xdata)
    push!(spl,Spline1D(xdata[i],ydata[i]))
  end

  return spl
end


"""
    get_ydata(xrange,xdata,ydata,spl)

From the `xrange` (unified x data stored in a vector), the original `xdata` stored
in an array of vectors and the respective `ydata` as well as cubic splines `spl`
of each dataset, return a matrix with unified ydata using the `xrange` and
supplementing missing ydata (due to different step sizes in the x datasets) within
the common range with estimates from the splines and fill missing data outside the
common range with `NaN`s.
"""
function get_ydata(xrange,xdata,ydata,spl)
  # Initialise output matrix
  unifiedy = Matrix{Float64}(undef,length(xrange),0)
  # Loop over all DataFrames
  for i = 1:length(ydata)
    # Initialise y data of current DataFrame
    yd = Float64[]
    # Loop over unified x data
    for x in xrange
      # Try to add value for current x value, if y is missing
      # interpolate with cubic spline, if inside the range of the current vector
      # or fill with NaN outside its range
      try push!(yd,ydata[i][findfirst(xdata.==x)])
      catch
        if xdata[i][1] ≤ x ≤ xdata[i][end]
          push!(yd,spl[i](x))
        else
          push!(yd,NaN)
        end
      end
    end
    # Add a column with completed y data of current DataFrame to the output matrix
    unifiedy = hcat(unifiedy,yd)
  end

  # Return unified y data
  return unifiedy
end #function get_ydata


"""
    get_ScalingFactor(x, y, mins, maxs, output)

From the unified `x` and `y` data and the indices of the start and end of vectors,
calculate and return either offsets or scaling factors (as defined in `output` with
the keywords `"offset"` or `"scale"`) to allow a smooth continuous mean at the
delayed start or premature ends of vectors.
"""
function get_ScalingFactor(x, y, mins, maxs, output)
  F = ones(length(x))
  if length(mins) ≥ 2
    for i = 1:length(mins)-1
      if output == "offset"
        F[mins[i]:mins[i+1]-1] .= mean(filter(!isnan,y[mins[i+1],:])) -
          mean(y[mins[i+1],broadcast(!,isnan.(y[mins[i+1]-1,:]))])
      elseif output == "scale"
        println(F)
        println(y)
        println(isnan.(y[mins[i+1]-1,:]))
        F[mins[i]:mins[i+1]-1] .= mean(filter(!isnan,y[mins[i+1],:])) /
          mean(y[mins[i+1],broadcast(!,isnan.(y[mins[i+1]-1,:]))])
      end
    end
  end
  if length(maxs) ≥ 2
    for i = 1:length(maxs)-1
      if output == "offset"
        F[maxs[i]+1:maxs[i+1]] .= mean(filter(!isnan, y[maxs[i],:])) -
          mean(y[maxs[i],:][broadcast(!,isnan.(y[maxs[i]+1,:]))])
      elseif output == "scale"
        F[maxs[i]+1:maxs[i+1]] .= mean(filter(!isnan, y[maxs[i],:])) /
          mean(y[maxs[i],:][broadcast(!,isnan.(y[maxs[i]+1,:]))])
      end
    end
  end

  return F
end #function get_ScalingFactor

end #module VectorMaths
