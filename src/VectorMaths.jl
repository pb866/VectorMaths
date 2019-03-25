"""
# Module VectorMaths

Data analysis and manipulation tools currently with function `vector_maths` to
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
    vector_maths(vectors::DataFrames.DataFrame...;
      output::String="offset", colnames::Union{String,Vector{Symbol}}="default",
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
  `offset + data` to avoid unnecessary steps in the output.
- `"scale"` (**default**): Mean and sum outside the common range are calculated
  by using a scaling factor `SF · data` to avoid unnecessary steps in the output.

    offset = Ø(y values of all vectors)-Ø(y values of remaining vectors)
    SF = Ø(y values of all vectors)/Ø(y values of remaining vectors)

This means, that the ends of each shorter vector will influence the mean outside
the range of this vector. Scaling factors should not be used, if the y data spans
a large range (over several orders of magnitudes) as otherwise a small offset at
the lower end of the data could lead to a large scaled offset at the upper end of
the y data. Scaling factors should be used, if you want to assure that the sign
in your y data does not change and the y data is in a reasonably small range.

Use kwarg `colnames` to define headers of the `output` DataFrame with the following
options:
- `"default"`: By default `colnames` is set to the keyword `"default"`, which assigns
  the following column names: `[:x, :y1, :y2,..., :mean, :sum]`.
- `"original"`: The original header of each DataFrame column are used. For the combined
  x data, the column name of the first column definition is used, duplicate header
  names are appended by `_1`, `_2`, ..., `:mean`, and `:sum` are used for the last
  2 columns.
- `Vector{Symbol}`: Directly assign column header names in a `Vector{Symbol}` defining
  the name of the combined x data in the first entry, followed by the names of the
  y columns, the mean, and sum: `[:x, :y1, :y2,..., :Ø, :∑]`
"""
function vector_maths(vectors::DataFrames.DataFrame...;
  output::String="offset", colnames::Union{String,Vector{Symbol}}="default",
  xcols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[],
  ycols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[])

  # Check for correct keywords in output
  if output != "offset" && output != "scale"
    error("Wrong keyword argument!\nUse \"offset\" or \"scale\" for `output`.")
  end

  # Retrieve x and y data in each DataFrame
  xdata, ydata, colnames = compile_data(vectors, xcols, ycols, colnames)

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
  # Assign column headers to output DataFrame
  names!(dfr, colnames)
  # Return output DataFrame
  return dfr
end #function vector_maths


###########################
###  PRIVATE FUNCTIONS  ###
###########################


"""
    compile_data(vectors, xcols, ycols, colnames)

For each DataFrame in `vectors`, find the x columns defined in `xcols` or use
first columns by default and the y columns defined in `ycols` or use second columns
as default and return vectors with x and y data vectors.
Return revised column headers according to the keyword `colnames`.
"""
function compile_data(vectors, xcols, ycols, colnames)

  # Define x and y data in each DataFrame and header names for output
  xcols = check_columns(vectors, xcols, 'x')
  ycols = check_columns(vectors, ycols, 'y')
  colheaders = get_headers(colnames, xcols, ycols)
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

  return xdata, ydata, colheaders
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
  coltype == 'x' ? column = 1 : column = 2
  if isempty(cols)
    # Use first column for x and 2nd for y data as default, if cols is empty
    cols = [names(v)[column] for v in vectors]
  else
    # If cols is specified, enforce Symbols over integers for naming of output data
    if typeof(cols) != Vector{Symbol}
      columns = []
      for (n, col) in enumerate(cols)
        if typeof(col) == Symbol
          push!(columns, col)
        elseif typeof(col) == Vector{Symbol}
          push!(columns, [c for c in col]...)
        elseif typeof(col) == Int
          push!(columns, names(vectors[n])[col])
        elseif typeof(col) == Vector{Int}
          push!(columns, [names(vectors[n])[c] for c in col])
        end
      end #loop over columns
      cols = columns
    end #non-Symbol columns
  end
  # Check and adjust correct length of column header definitions
  if length(cols) > length(vectors)
    println("\033[95mWarning! More $coltype columns defined than vectors available.\33[0m")
    println("Last $(length(cols) - length(vectors)) columns in `$(coltype)cols` ignored.")
    deleteat!(cols, 1 + length(vectors):length(cols))
  elseif length(cols) < length(vectors)
    println("\033[95mWarning! Less x columns defined than vectors available.\33[0m")
    println("$column. column assumed to hold $coltype data in last ",
      "$(length(vectors) - length(cols)) DataFrames.")
    push!(cols, [names(vectors[i])[column] for i = 1 + length(cols):length(vectors)]...)
  end

  return cols
end #function check_columns


"""
    get_headers(colnames, xcols, ycols)

Return a vector of symbols with the column headers according to the keyword `colnames`
from the headers stored for `xcols` and `ycols`.
"""
function get_headers(colnames, xcols, ycols)
  # Flatten arrays with column names or indices
  x = vcat(xcols...); y = vcat(ycols...)
  # Check correct assignment of colnames, if given as Vector of header names
  # and throw error or leave function
  if colnames isa Vector
    length(colnames) != length(y) + 3 ?
      error(join(["Wrong assignment of column headers.\n",
      "Use `Vector{Symbol}` with x column names as first entry, ",
      "followed by the y column names and names for the mean and sum."])) :
        return colnames
  end
  if colnames == "default"
    ### Define standard names as x, y1, y2, y3,...
    colheaders = [:x]
    push!(colheaders, [Symbol("y$i") for i = 1:length(y)]...)
    push!(colheaders, [:mean, :sum]...)
  elseif colnames == "original"
    ### Use original column names, append duplicate headers with "_1", "_2", ...
    # Use name of first x column for combined x data
    colheaders = [x[1]]
    # Loop over all y columns
    for (n, column) in enumerate(y)
      # Get header of current y column
      global col = column
      # Search for duplicate names in already defined headers
      global duplicate = match.(Regex("($col)(_[0-9])*\$"),string.(colheaders))
      # Remove already existing _#, where # is a number
      if !(all(isnothing.(duplicate)))
      c = match(Regex("(.*)_[0-9]\$"),string(col))

      if !isnothing(c)  col = c.captures[1] end
        # Search for already existing duplicates and the appended number last used
        dup = findall(broadcast(!,isnothing.([d for d in duplicate])))
        d = findlast(broadcast(!,isnothing.([duplicate[d].captures[2] for d in dup])))
        if !isnothing(d)
          # Save last appended number and the core header name
          global i = parse(Int, duplicate[dup[d]].captures[2][2:end])
          header = duplicate[dup[d]].captures[1]
        else
          # Set index to 0, if no other duplicates were previously found
          global i = 0
          # Save core header name
          header = duplicate[dup[end]].match
        end
        # Increase appendix number until a previously unused version is found
        while any(broadcast(!,isnothing.(duplicate)))
          # Increase appendix number
          global i += 1
          # Search for duplicates with that number
          global duplicate = match.(Regex("($header)(_$i)\$"),string.(colheaders))
          # Save new header name to be stored in colheaders
          global col = Symbol("$(header)_$i")
        end
      end
      # Save header in array of all column headers
      push!(colheaders,col)
    end
    push!(colheaders, [:mean, :sum]...)
  else
    error(join(["Wrong assignment of column headers.\n",
      "Use `Vector{Symbol}` with x column names as first entry ",
      "and y column names in the following.\nOr use keywords \"original\" ",
      "to assign the header names of the given DataFrames or ",
      "keyword \"default\" to assign default names `x, y1, ... yn'."]))
  end

  return colheaders
end #function get_headers

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
    test_monotonicity(xdata)

Test that all `x` column data in `xdata` is strictly monotonic ascending.
Issue an error, if not.
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
