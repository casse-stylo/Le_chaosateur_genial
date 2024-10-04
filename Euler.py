
import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

def Euler(wn, f, h, pot) :
    wn = wn + h*f(wn, pot)
    return wn 

