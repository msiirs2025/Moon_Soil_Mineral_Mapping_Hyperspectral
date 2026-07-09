import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import OUT, IIRS_L, IIRS_S
import ch2_seleno_destripe as SD

C2 = os.path.join(OUT, "ch2")
base = os.path.join(C2, "step5_destripe", "ch2_iirs_step5_destriped_selenoref")
# remove partial outputs
for f in glob.glob(base + ".*"):
    os.remove(f); print("removed partial", os.path.basename(f))

src = np.memmap(os.path.join(C2, "step5_destripe", "ch2_iirs_step5_destriped_sensor.img"),
                dtype='<f4', mode='r', shape=(SD.NB, IIRS_L, IIRS_S))
print("georeferencing destriped sensor cube ...")
SD._warp_cube(src, base,
              "Destriped + selenoreferenced reflectance (final product, Moon2000 eqc)")
print("DONE step5 georeference")
