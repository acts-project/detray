# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from collections import namedtuple

#-------------------------------------------------------------------------------
# Common helpers for plotting measurement data
#-------------------------------------------------------------------------------

""" Filter the data in a data frame by a given prescription """
def filter_data(data, filter = lambda df: [], variables = []):
    dataColl = []

    # Get global data
    if len(filter(data)) == 0:
        for var in variables:
            dataColl.append(data[var].to_numpy())

    # Filtered data
    else:
        filtered = data.loc[filter]
        for var in variables:
            dataColl.append(filtered[var].to_numpy())

    return tuple(dataColl)

#-------------------------------------------------------------------------------
