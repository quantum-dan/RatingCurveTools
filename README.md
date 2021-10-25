# Rating Curve Tools

Some utility scripts for creating and applying rating curves.  RatingCurve.R generates rating curves fitted to data
and fast_ftable.py efficiently applies rating curves to a large volume of timeseries (or other input) data.

Note that this is NOT an actively-maintained project and is not consciously written to be generic; the scripts are developed as
necessary for my own use and tailored to my own use cases.  Use at your own risk.

## RatingCurve.R

Takes a table of flow (or other) data and fits rating curve functions to it.  The input table should have
the first two columns named `ID` (location identifier) and `x` (x variable, e.g. Q).  The remaining columns
can be an arbitrary number of arbitrarily-named output functions.

Given such a table, `fit.curves` fits rating curves to all the given functions and identifiers using the
1st through 10th powers of X (outputted as `X1`...`X10`), natural log `LnX`, square root `SqrtX`, and cube
root `Rt3X`.  It then outputs a table of the function coefficients for each ID and function, where the rating
curve function as a whole is the sum of the results of multiplying each coefficient by the corresponding
function of x.  It also outputs the NSE and RMSE of the fit and the MinQ and MaxQ, identifying, respectively,
the lowest value of x in the input data or the lowest value of x corresponding to a positive value of y,
and the highest value of x in the input data.  During the fitting, it can show fit plots to the user or save them
to a specified file `save.plot`.

The output table, written to a CSV, can be used directly as the input to fast_ftable.

## fast_ftable.py

The primary functionality is handled by `run`.

fast_ftable efficiently applies rating curves to large timeseries or other data sets.  It can accept arbitrary
metadata and preserve it into the output, as it ignores all columns except `total_flow`.

`run` takes an inputs tuple specifying the ID, flow file path, curve file path, and output file path, in that order.
The curve file should be in the format of the output from RatingCurve.R.  The ID identifies which set of functions
to apply, corresponding to the ID column in the curve table.  The flow file must be a CSV, but the format is irrelevant
except that it must contain a `total_flow` column.  All other columns are simply kept with the corresponding flow.

Note that the current implementation assumes that `total_flow` is in cubic meters per second and converts it to cubic feet
per second, due to my use case.  If you are using it under different assumptions, you should modify `parse_ts` to remove
or modify the conversion as appropriate.
