import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import BASE
import ch2_pipeline as P

inc = P.gen_angle()

# existing reference _angle (250 x 13569, float32 bsq, 1 band)
ref = os.path.join(BASE, r"ch2\ch2\ch2\ch2_iir_nri_20250804T0253259678_d_img_d18_angle")
if os.path.exists(ref):
    r = np.fromfile(ref, dtype='<f4').reshape(P.IIRS_L, P.IIRS_S)
    d = inc - r
    print(f"\nvs existing _angle:  mean|diff|={np.mean(np.abs(d)):.4f}  max|diff|={np.max(np.abs(d)):.4f}")
    print(f"  ref range {r.min():.3f}..{r.max():.3f}   mine {inc.min():.3f}..{inc.max():.3f}")
else:
    print("no existing _angle to compare")
