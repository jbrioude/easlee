# Compare two EASLEE grid outputs (last frame of `conc`): total mass and Pearson correlation.
# Used for the backend validations (threads vs ka_gpu) in docs/BENCHMARKS.md.
import sys, glob, numpy as np, netCDF4 as nc, os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
def load(d):
    f = sorted(glob.glob(f'{d}/grid_*.nc'))[-1]; c = nc.Dataset(f); v = c.variables['conc']
    x = np.ma.filled(v[-1], 0).astype(np.float64); return f, x
fa, a = load(sys.argv[1]); fb, b = load(sys.argv[2])
ma, mb = a.sum(), b.sum(); r = np.corrcoef(a.ravel(), b.ravel())[0, 1]
print(f"{os.path.basename(fa)} vs {os.path.basename(fb)}: mass {ma:.6e} vs {mb:.6e} ({abs(mb/ma-1):.2e} apart), correlation {r:.6f}")
