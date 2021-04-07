import pandas as pd
import matplotlib.pyplot as plt

Etot = pd.read_table('etot.dat', delim_whitespace=True)
Etot.columns=["Time","Etot"]
Etot.plot(x ='Time', y='Etot', kind = 'line')