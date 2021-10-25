"""
Written beginning October 11, 2020 by Daniel Philippus for the Los Angeles
River Environmental Flows project at the Colorado School of Mines.  Revised
beginning October 7, 2021 for use with the improved rating curve setup.

All functionality is handled by the `run` function.

This program is an optimized re-implementation of ftable.R; it applies fitted
rating curve tables to flow timeseries data.  It is especially intended for use
with multiple scenarios.  The program is being rewritten using numpy because
the original R implementation was unacceptably slow for large volumes of data;
although it could be sufficiently faster just by using matrix operations in R,
it was decided to use numpy as well to hedge our bets.

Flow files are assumed to have a column total_flow; any remaining
columns are simply preserved.

The rating curve file should have the format
    ID,Function,nse,rmse,X1...X10,LnX,SqrtX,Rt3X
With X1 onwards describing coefficients.

1. Convert into a matrix with rows for x, x^2, sqrt(x), etc.
2. Multiply the coefficient vectors by the matrix.  Result: a matrix with rows
    corresponding to output variables and columns corresponding to flows.
3. Transpose matrix and add in the artefact columns.
4. Write CSV.

The final format will be <artefact columns>, <functions...>
"""

from numpy import array, transpose, real
# from glob import glob
from math import sqrt, log

BASE_PATH = r"C:\Users\dphilippus\Dropbox\LA-River Flows-SCCWRP\05-modeling\restoration"
CURVE_PATH = BASE_PATH + r"\RatingCurves\Curves"
FLOW_PATH = BASE_PATH + r"\Restoration_flows"
OUT_PATH = BASE_PATH + r"\RatingCurves\Outputs"

CURVE_COEFFS = ["X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10",
                "LnX", "SqrtX", "Rt3X"]
CURVE_IDC = "ID"
CURVE_FC = "Function"

CURVE_FNS = [
    lambda x: x,
    lambda x: x**2,
    lambda x: x**3,
    lambda x: x**4,
    lambda x: x**5,
    lambda x: x**6,
    lambda x: x**7,
    lambda x: x**8,
    lambda x: x**9,
    lambda x: x**10,
    lambda x: log(x if x > 0 else 0.001),
    lambda x: sqrt(x if x > 0 else 0),
    lambda x: x**(1/3)
    ]


def read_csv(path):
    with open(path, "r") as f:
        return [line.strip().split(",") for line in f]


def write_csv(x, path):
    with open(path, "w") as f:
        f.write("\n".join([",".join(map(str, r)) for r in x]))


def parse_ts(path):
    # Convert CSV into a flow vector and artefact columns
    raw = read_csv(path)
    flow_ix = raw[0].index("total_flow")
    flows = [float(r[flow_ix]) * 3.28**3 for r in raw[1:]]  # convert to cfs
    for ix in range(1, len(raw)):
        raw[ix][flow_ix] = flows[ix-1]  # Switch to cfs
    return (flows, raw)


def prepare_ftab(path):
    # Prepare function table for application
    raw = read_csv(path)
    header = raw[0]
    # Validate function table
    idx = header.index(CURVE_IDC)
    fnx = header.index(CURVE_FC)
    hix = header.index(CURVE_COEFFS[0])
    hjx = header.index(CURVE_COEFFS[-1])
    if not header[hix:hjx+1] == CURVE_COEFFS:
        raise ValueError("Wrong set of curve coefficients")
    # Prepare for application: dictionary of
    # {id: (function names, functions matrix)}
    # So that for a given ID the matrix can be multiplied by the inputs
    output = {}
    for row in raw[1:]:
        ident = row[idx]
        if ident in output:
            output[ident][0].append(row[fnx])
            output[ident][1].append(list(map(float, row[hix:hjx+1])))
        else:
            output[ident] = ([row[fnx]], [list(map(float, row[hix:hjx+1]))])
    return output


def apply_ftab(ftab, flows, raw, id):
    # ftab from prepare_ftab, flows and raw from parse_ts
    # id associated with ftab identifier
    # Apply function table to flows
    #
    # First: expand flows; this gives a matrix with rows
    # corresponding to input functions, which we can multiply by
    # the coefficient matrix
    exflows = array([
            list(map(lambda x: real(CURVE_FNS[ix](x)), flows))
            for ix in range(len(CURVE_FNS))
        ])
    # Next: coeff is a matrix with rows corresponding to output functions
    # and columns corresponding to input functions
    coeff = array(ftab[id][1])
    # Product gives flow predictions
    prediction = coeff @ exflows
    # Now just re-arrange things and label it
    fnames = ftab[id][0]
    vertical = transpose(prediction).tolist()
    named = [fnames] + vertical
    # Put raw back at the front to include all the labels etc
    return [(["ID"] if ix == 0 else [id]) + raw[ix] + named[ix]
            for ix in range(len(raw))]


def run(inputs, basepaths=True):
    # Run all input files
    # inputs: [(id, flow path, curve path, output path)]
    # if basepaths is true, prepend base paths to all three; otherwise, user
    # must provide full, suitable paths
    pathify = lambda bp: bp + "\\" if basepaths else ""
    for (id, fp, cp, op) in inputs:
        (flows, raw) = parse_ts(pathify(FLOW_PATH) + fp)
        ftab = prepare_ftab(pathify(CURVE_PATH) + cp)
        out = apply_ftab(ftab, flows, raw, id)
        write_csv(out, pathify(OUT_PATH) + op)


base_name_13m = "Sc1013.minuswidth\\%s_Sc1013_n035_z4_afp_lfc13.csv"
base_name_13p = "Sc1013.pluswidth\\%s_Sc1013_n035_z4_afp_lfc13.csv"
curve_name_13m = "Curve_Sc1013m.csv"
curve_name_13p = "Curve_Sc1013p.csv"
cn40 = "Curve_Sc40_n035_z4_afp_lfc13.csv"
cn11 = "Curve_Sc1011.csv"
cn_iter = "Curve_Sc1014_iterations.csv"
base_iter = "Sc1014_iter\\%s_Sc1014_n035_z4_afp_lfc13.csv"

inputs_iter = [
    ("F45B", "F45B_all_scenarios.csv", cn_iter,
     base_iter % "F45B"),
    ("F37B",
     "F37B_all_scenarios.csv", cn_iter,
     base_iter % "F37BHigh"),
    ("LA14",
     "LA14_all_scenarios.csv", cn_iter,
     base_iter % "LA14")
    ]

inputs13m = [
    ("F45B", "F45B_all_scenarios.csv", curve_name_13m,
     base_name_13m % "F45B_w443"),
    ("F37BHigh",
     "F37B_all_scenarios.csv", curve_name_13m,
     base_name_13m % "F37BHigh_w246"),
    ("LA14",
     "LA14_all_scenarios.csv", curve_name_13m,
     base_name_13m % "LA14_w920")
    # ("LA3",
    #  "LA3_all_scenarios_iterations.csv", curve_name, "LA3_w1084" + base_name)
    ]

inputs13p = [
    ("F45B", "F45B_all_scenarios.csv", curve_name_13p,
     base_name_13p % "F45B_w607"),
    ("F37BHigh",
     "F37B_all_scenarios.csv", curve_name_13p,
     base_name_13p % "F37BHigh_w410"),
    ("LA14",
     "LA14_all_scenarios.csv", curve_name_13p,
     base_name_13p % "LA14_w1248")
    ]

inputs11 = [
    ("F45B", "F45B_all_scenarios.csv", cn11,
     "Sc1011\\F45B_w525_Sc1011_n035_z4_afp8_lfc13.csv"),
    ("F37BHigh",
     "F37B_all_scenarios.csv", cn11,
     "Sc1011\\F37BHigh_w328_Sc1011_n035_z4_afp8_lfc13.csv")
    ]

inputs40 = [
    ("F45B", "F45B_all_scenarios.csv", cn40,
     "Sc40\\F45B_w39_Sc40_n035_z4_afp_lfc13.csv"),
    ("F37BHigh",
     "F37B_all_scenarios.csv", cn40,
     "Sc40\\F37BHigh_w39_Sc40_n035_z4_afp_lfc13.csv"),
    ("LA14",
     "LA14_all_scenarios.csv", cn40,
     "Sc40\\LA14_w1084_Sc40_n035_z4_afp_lfc13.csv")
    ]
