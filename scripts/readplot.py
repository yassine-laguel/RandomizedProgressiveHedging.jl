import json
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib2tikz
import tkinter as tk
from tkinter import filedialog


root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename()

f = open(file_path)
data = json.load(f)
f.close()


results = data["simpleproblem"] # RESULTS PARTS
T = data["nstages"]
S = data["nscenarios"]
W = data["nworkers"]

fmin = results["progressivehedging"]["1"]["fopt"]

names = { "progressivehedging" : "Progressive Hedging" , "randomized_sync": "Randomized Progressive Hedging"  , "randomized_async": "Asynchronous Randomized Progressive Hedging" }
colors = { "progressivehedging" : "black" , "randomized_sync": "blue"  , "randomized_async": "red" }

########################################
### SUBOPTIMALITY vs CALLS
########################################
plt.figure()
for alg in {"progressivehedging","randomized_sync","randomized_async"}:
    RES = results[alg]["1"]
    F = [f- fmin for f in RES["functionalvalue"]]
    if alg == "progressivehedging":
        C = [S*i for i in range(len(F))]
    else:
        C = range(len(F))
    plt.plot(C,F,label=names[alg],color=colors[alg])
plt.ylabel("Suboptimality")
plt.xlabel("Number of scenarios treated")
plt.yscale('log', nonposy='clip')
plt.legend()
plt.savefig("./Figs/Subopt_Calls.png")
#matplotlib2tikz.save("./Tex/Subopt_Calls.tex")


########################################
### SUBOPTIMALITY vs TIME
########################################
plt.figure()
for alg in {"progressivehedging","randomized_sync","randomized_async"}:
    RES = results[alg]["1"]
    F =  [f- fmin for f in RES["functionalvalue"]]
    T = RES["time"]
    plt.plot(T,F,label=names[alg],color=colors[alg])
plt.ylabel("Suboptimality")
plt.xlabel("Time (s)")
plt.yscale('log', nonposy='clip')
plt.legend()
plt.savefig("./Figs/Subopt_Time.png")
#matplotlib2tikz.save("./Tex/Subopt_Time.tex")


########################################
### SUBOPTIMALITY/FILL vs CALLS
########################################
plt.figure()

fig,ax = plt.subplots()
for alg in {"progressivehedging","randomized_sync","randomized_async"}:
    Nruns = 2
    F = list()
    for i in range(Nruns):
        F.append([])
    
    lmin = np.inf
    for i in range(Nruns):
        RES = results[alg][str(i+1)]
        F[i] = [f- fmin for f in RES["functionalvalue"]]
        lmin = min(lmin,len(F[i]))
    
    Fmin = []
    Fmax = []
    for i in range(lmin):
        vmin = np.inf
        vmax = -np.inf
        for j in range(Nruns):
            vmin = min(vmin,F[j][i])
            vmax = max(vmax,F[j][i])
        Fmin.append(vmin)
        Fmax.append(vmax)
    

    if alg == "progressivehedging":
        C = [S*i for i in range(len(Fmin))]
    else:
        C = range(len(Fmin))
    ax.plot(C,Fmin,label=names[alg],color=colors[alg])
    ax.plot(C,Fmax,color=colors[alg])
    ax.fill_between(C,Fmin,Fmax,facecolor=colors[alg], alpha=0.5)
plt.ylabel("Suboptimality")
plt.xlabel("Number of scenarios treated")
plt.legend()
plt.yscale('log', nonposy='clip')
plt.savefig("./Figs/SuboptFill_Calls.png")
#matplotlib2tikz.save("./Tex/SuboptFill_Calls.tex")

