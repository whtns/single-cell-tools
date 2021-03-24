import loompy

def cellnames_from_loom(loom_path):
  with loompy.connect(loom_path) as ds:
    print(ds.ca["CellID"])
