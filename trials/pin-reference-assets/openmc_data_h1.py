"""openmc.data で H-1 の ENDF から NJOY を通して HDF5 を作る（ADR-0124 の確認。openmc_library.py と同じ道筋）。

使い方: python3 openmc_data_h1.py <ENDF の H-1 のファイル> <出力の HDF5>
"""

import sys

import openmc.data

data = openmc.data.IncidentNeutron.from_njoy(sys.argv[1], temperatures=[600.0])
data.export_to_hdf5(sys.argv[2], "w")
