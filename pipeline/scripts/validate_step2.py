import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import BASE, OUT, IIRS_L, IIRS_S

NB = 247
mine = np.memmap(os.path.join(OUT, r"ch2\step2_thermal\ch2_iirs_step2_reflectance_thermcorr.img"),
                 dtype='<f4', mode='r', shape=(NB, IIRS_L, IIRS_S))
inc = np.load(os.path.join(OUT, r"ch2\aux\angle.npy"))
myTemp = np.fromfile(os.path.join(OUT, r"ch2\aux\ch2_iirs_Temp.img"), '<f4').reshape(IIRS_L, IIRS_S)

d = BASE + r"\ch2\ch2\ch2\ch2_iir_nri_20250804T0253259678_d_img_d18"
refIC = d + "_inc_corrRef"
refT  = d + "_Temp"

# Temperature
rT = np.fromfile(refT, '<f4').reshape(IIRS_L, IIRS_S)
dT = myTemp - rT
print(f"Temperature: mean|d|={np.nanmean(np.abs(dT)):.4f}K max|d|={np.nanmax(np.abs(dT)):.3f}K  "
      f"ref {rT.min():.1f}-{rT.max():.1f} mine {myTemp.min():.1f}-{myTemp.max():.1f}")

# _inc_corrRef == my corrRef / cos(inc)
ric = np.memmap(refIC, dtype='<f4', mode='r', shape=(NB, IIRS_L, IIRS_S))
cosi = np.cos(np.deg2rad(inc))
# compare a subset of bands to save memory
for b in [0, 60, 120, 180, 246]:
    myIC = mine[b] / cosi
    r = ric[b]
    m = np.isfinite(myIC) & np.isfinite(r) & (r != 0)
    dd = np.abs(myIC[m] - r[m])
    rel = dd / (np.abs(r[m]) + 1e-9)
    print(f"band {b:3d}: mean|d|={dd.mean():.3e} max|d|={dd.max():.3e} "
          f"medRel={np.median(rel)*100:.3f}%  ref~{np.nanmedian(r[m]):.4f}")
