VectorMaths
===========

Data analysis and manipulation tools currently with function `vector_maths` to
calculate the mean and sum of selected vectors stored in DataFrames.

Installation
------------

Install the Package by adding it to your environment going to the package
manager typing `]` in the julia prompt and in the package manager the following commands:

```julia
julia> ]
pkg> add https://github.com/pb866/filehandling.git
pkg> add https://github.com/pb866/VectorMaths.git
```

`VectorMaths` relies on the unregistered packages `filehandling`, which has to be
installed beforehand.


Usage
-----

### Function vector_maths
Function `vector_maths` calculates the mean and sum of any vectors stored in DataFrames.

```julia
vector_maths(vectors::DataFrames.DataFrame...;
  output::String="offset", colnames::Union{String,Vector{Symbol}}="default",
  xcols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[],
  ycols::Union{Vector{Int64},Vector{Symbol},Vector{Any}}=[])
```

Hand over any number of DataFrames stored in vararg `vectors` to `vector_maths`
and calculate the mean and sum of vectors from the stored x and y data in these
DataFrames.
A DataFrame is return with combined `x` data as first column, ydata (`y1`, `y2`, ...)
in the following columns, where missing data outside the common range is filled
with `NaN`s, and a column with the `mean` data and `sum` of these data.


#### Selecting columns

By default, the first column in each DataFrame is treated as x data and the second
column as y data. This can be changed with the kwargs `xcols` and `ycols`, which
specify the x and y columns in a DataFrame as follows:

`xcols`/`ycols` must be vectors holding either the position of the column in the
DataFrame as integer or the column name as `Symbol` (mixtures are also possible).
If you want to select more than one column from one DataFrame, use vectors of
integers or symbols as the respective entry in `xcols`/`ycols`. You can either
specify the same amount of x and y columns in paired `xcols`/`ycols` entries or,
if the y data in a DataFrame shares the same x data, you can specify the x data
with one integer or symbol and the respective y data columns as vector of integers
or symbols.


#### Selecting output format

To calculate the mean and sum, the following options exist, which can be specified
by the kwarg `output`:
- `"common"`: If vectors have different lengths, only for the overlapping range
a mean and sum is calculated and outputted.
- `"offset"`: Mean and sum outside the common range are calculated by using an offset
`offset + data` to avoid unnecessary steps in the output.
- `"scale"` (**default**): Mean and sum outside the common range are calculated
by using a scaling factor `SF · data` to avoid unnecessary steps in the output.

    offset = Ø(y values of all vectors) - Ø(y values of remaining vectors)
    SF = Ø(y values of all vectors) / Ø(y values of remaining vectors)

This means, that the ends of each shorter vector will influence the mean outside
the range of this vector. Scaling factors should not be used, if the y data spans
a large range (over several orders of magnitudes) as otherwise a small offset at
the lower end of the data could lead to a large scaled offset at the upper end of
the y data. Scaling factors should be used, if you want to assure that the sign
in your y data does not change and the y data is in a reasonably small range.

Use the kwarg `colnames` to define headers of the `output` DataFrame with the
following options:
- `"default"`: By default `colnames` is set to the keyword `"default"`, which assigns
the following column names: `[:x, :y1, :y2,..., :mean, :sum]`.
- `"original"`: The original header of each DataFrame column are used. For the combined
x data, the column name of the first column definition is used, duplicate header
names are appended by `_1`, `_2`, ..., `:mean`, and `:sum` are used for the last
2 columns.
- `Vector{Symbol}`: Directly assign column header names in a `Vector{Symbol}` defining
the name of the combined x data in the first entry, followed by the names of the
y columns, the mean, and sum: `[:x, :y1, :y2,..., :Ø, :∑]`.


Version history
===============

Version 0.1.1
-------------
- Improved error handling with new error types `ColError` and `DataError`
- Fix error handling of mis-assigned x and y columns
- Restructured source code

Version 0.1.0
-------------
- Improved version of `VectorMaths` in [auxdata](https://github.com/pb866/auxdata.git) with language updated to Julia 1.1.
- Combine functions `DFvector_mean` and `multivector_mean` from previous module
- Use vararg instead of vector of DataFrames for input DataFrames
- Updated kwargs with `xcols`/`ycols` to select x/y data columns, `output` to select data range and scaling of data outside the common range, and `colnames` to manipulate column headers
