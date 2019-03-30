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

export vector_maths,
       ColError, DataError

### Error handling

## Column Errors
# Define Error type
struct ColError <: Exception
  alarm::String
  info::String
end
# Format Error message
Base.showerror(io::IO, e::ColError) = print(io, typeof(e), "\n", e.alarm, "\n\033[0m", e.info)
# Define default alarm and info message
ColError(msg) = ColError("Wrong assignment of column headers.", msg)
ColError() = ColError("Wrong assignment of column headers.",
  string("Use `Vector{Symbol}` with x column names as first entry, ",
  "followed by the y column names and names for the mean and sum."))

## DataError
# Define Error type
struct DataError <: Exception
  msg::String
  data
end
# Format Error message
Base.showerror(io::IO, e::DataError) = print(io, typeof(e), "\n", e.msg, "\n\033[0m", e.data)
# Define default alarm and info message
DataError(msg) = DataError(msg, nothing)
DataError() = DataError("", nothing)

# Include functions for vector_maths
include("./vector_mean_and_sum.jl")


end #module VectorMaths
