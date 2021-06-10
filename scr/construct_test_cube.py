import numpy as np

import sys

num='123456'

T = np.zeros((512, 512, 180))
p = np.zeros((512, 512, 180))
d = np.zeros((512, 512, 180))

z = np.arange(180) * 10.0

T0l = 4500.0
#T0r = T0l
T0r = 6500.0

#b = 3.0
b = 0.0

Tl = b * z + T0l
Tr = b * z + T0r

d0l = 45000000
#d0r = d0l
d0r = 65000000

#a = 10000.0
a = 0.0

#dl = (a * z + d0l) * 1e-16
#dr = (a * z + d0r) * 1e-16
dl = 1.2
dr = 1.6

m = 1.27098 * 1.66e-24

k = 1.38e-16

pl = dl * k * Tl / m
pr = dr * k * Tr / m

T[:256, : , :] = Tl

p[:256, :, :] = pl

d[:256, :, :] = dl

T[256:, : , :] = Tr

p[256:, :, :] = pr

d[256:, :, :] = dr

#sys.exit()

T = T.flatten().astype(np.float32)
p = p.flatten().astype(np.float32)
d = d.flatten().astype(np.float32)

T.tofile('eosT.'     + num + '.bin')
p.tofile('eosP.'     + num + '.bin')
d.tofile('result_0.' + num + '.bin')
