import numpy as np, os
qub = r"D:\new project moon\ch2\ch2\ch2\ch2_iir_nri_20250804T0253259678_d_img_d18\data\raw\20250804\ch2_iir_nri_20250804T0253259678_d_img_d18.qub"
S,L,B = 250, 13569, 256
sz = os.path.getsize(qub)
print("file size:", sz, " expected:", S*L*B*2, " match:", sz==S*L*B*2)
# BSQ uint16 LE. memmap
arr = np.memmap(qub, dtype='<u2', mode='r', shape=(B,L,S))
# sample every 200th line to be fast
samp = arr[:, ::200, :].astype(np.float64)   # (B, ~68, S)
gmin, gmax = arr[:, ::400, :].min(), arr[:, ::400, :].max()
print("global sampled min/max:", gmin, gmax)
print("\nband  wl(approx)   min     mean      max     %zero")
wl0=712.3339889; dwl=(5009.686324-712.3339889)/255
for b in [0,7,30,60,80,100,130,160,190,220,246,253,255]:
    band = arr[b, ::100, :].astype(np.float64)
    z = 100.0*np.mean(band==0)
    print(f"{b:4d}  {wl0+b*dwl:8.1f}  {band.min():6.0f}  {band.mean():8.1f}  {band.max():6.0f}   {z:5.1f}")
# look for dark/reference rows (first & last lines) - IIRS often has dark frames
print("\nfirst 3 lines mean (band 100):", arr[100,:3,:].mean(axis=1))
print("last  3 lines mean (band 100):", arr[100,-3:,:].mean(axis=1))
