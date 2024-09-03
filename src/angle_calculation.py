import numpy as np
import math

E = 1e+6

############  OXYGEN ################

for i in range(89):

    radian = (i+1)*math.pi / 180
    sin_4th = (math.sin(radian))**4

    # print(sin_4th)

    const = 3.4e-54


    k = const / sin_4th

    verovatnoca_za_ugao = k / (E**2)

    print(verovatnoca_za_ugao*1e59)

    # print(math.sin(radian))