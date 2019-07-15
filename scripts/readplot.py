import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib2tikz
import tkinter as tk
from tkinter import filedialog
import os


root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename()

f = open(file_path)
data = json.load(f)
f.close()


# define the name of the directory to be created
path = "./results/"

try:  
    os.mkdir(path)
except OSError:  
    print ("Figures will be saved in %s " % os.path.realpath("./results/"))
else:  
    print ("Figures will be saved in %s (created)" % os.path.realpath("./results/"))



pbnames = data["problem_names"]

for pb in pbnames:

    results = data[pb] # RESULTS PARTS

    if "nstages" in data:
        T = data["nstages"]
        S = data["nscenarios"]
        # W = data["nworkers"]
    else:
        T = results["nstages"]
        S = results["nscenarios"]
        # W = results["nworkers"]


    if "algorithms" in results:
        algorithms = results["algorithms"]
    else:
        algorithms = results.keys()

    cmap = plt.get_cmap("nipy_spectral")
    colors_val = cmap(np.linspace(0, 1, len(algorithms)))

    names = {}
    colors = {}
    for (i,n) in enumerate(algorithms):
        names[n] = n
        colors[n] = colors_val[i]

    for alg in algorithms:
        if alg in results:
            fmin = results[alg]["1"]["fopt"]

            #RES = results[alg]["1"]
            #F = [f for f in RES["functionalvalue"]]
            #fmin = min(F)

    ########################################
    ### SUBOPTIMALITY vs CALLS 
    ########################################
    plt.figure()
    for alg in algorithms:
        if alg in results:
            RES = results[alg]["1"]
            step = results[alg]["1"]["logstep"]
            F = [abs((f - fmin)/fmin) for f in RES["functionalvalue"]]
            C = [step*i for i in range(len(F))]

            plt.plot(C,F,label=names[alg],color=colors[alg])

    plt.ylabel("Suboptimality")
    plt.xlabel("Number of scenarios treated")
    plt.yscale('log', nonposy='clip')
    plt.legend()
    plt.savefig(path+"Subopt_Calls"+pb+".png")
    matplotlib2tikz.save(path+"Subopt_Calls"+pb+".tex")


    ########################################
    ### SUBOPTIMALITY vs ALL
    ########################################
    plt.figure()
    for alg in algorithms:
        if alg in results:
            
            if type(results[alg]["seeds"]) == int:
                seeds = ['1']
                Nruns = 1
            else:
                seeds = [str(e) for e in results[alg]["seeds"]]
                Nruns = len(seeds)

            RES = results[alg]["1"]
            step = results[alg]["1"]["logstep"]
            F = [abs((f - fmin)/fmin) for f in RES["functionalvalue"]]
            C = [step*i for i in range(len(F))]

            plt.plot(C,F,label=names[alg],color=colors[alg])

            for run in range(2,Nruns):
                srun = seeds[run]
                RES = results[alg][srun]
                step = results[alg][srun]["logstep"]
                F = [abs((f - fmin)/fmin) for f in RES["functionalvalue"]]
                C = [step*i for i in range(len(F))]

                plt.plot(C,F,color=colors[alg])
    
    plt.ylabel("Suboptimality")
    plt.xlabel("Number of scenarios treated")
    plt.yscale('log', nonposy='clip')
    plt.legend()
    plt.savefig(path+"ALL_Subopt_Calls"+pb+".png")
    matplotlib2tikz.save(path+"ALL_Subopt_Calls"+pb+".tex")

    ########################################
    ### SUBOPTIMALITY vs TIME
    ########################################
    plt.figure()
    for alg in algorithms:
        if alg in results:
            RES = results[alg]["1"]
            F =  [(f- fmin)/fmin for f in RES["functionalvalue"]]
            T = RES["time"]

            plt.plot(T,F,label=names[alg],color=colors[alg])

    plt.ylabel("Suboptimality")
    plt.xlabel("Time (s)")
    plt.yscale('log', nonposy='clip')
    plt.legend()
    plt.savefig(path+"Subopt_Time"+pb+".png")
    matplotlib2tikz.save(path+"Subopt_Time"+pb+".tex")


    ########################################
    ### SUBOPTIMALITY vs TIME
    ########################################
    plt.figure()
    for alg in algorithms:
        if alg in results:

                        
            if type(results[alg]["seeds"]) == int:
                seeds = ['1']
                Nruns = 1
            else:
                seeds = [str(e) for e in results[alg]["seeds"]]
                Nruns = len(seeds)


            RES = results[alg]["1"]
            F =  [(f- fmin)/fmin for f in RES["functionalvalue"]]
            T = RES["time"]
            
            plt.plot(T,F,label=names[alg],color=colors[alg])


            for run in range(2,Nruns):
                srun = seeds[run]
                RES = results[alg][srun]
                step = results[alg][srun]["logstep"]
                F = [abs((f - fmin)/fmin) for f in RES["functionalvalue"]]
                T = RES["time"]

                plt.plot(T,F,color=colors[alg])

    plt.ylabel("Suboptimality")
    plt.xlabel("Time (s)")
    plt.yscale('log', nonposy='clip')
    plt.legend()
    plt.savefig(path+"ALL_Subopt_Time"+pb+".png")
    matplotlib2tikz.save(path+"ALL_Subopt_Time"+pb+".tex")


    ########################################
    ### SUBOPTIMALITY/FILL vs CALLS
    ########################################
    plt.figure()

    fig,ax = plt.subplots()
    for alg in algorithms:
        if alg in results:
            step = results[alg]["1"]["logstep"]

            if type(results[alg]["seeds"]) == int:
                seeds = ['1']
                Nruns = 1
            else:
                seeds = [str(e) for e in results[alg]["seeds"]]
                Nruns = len(seeds)
            
            F = list()
            for i in range(Nruns):
                F.append([])
            
            lmin = np.inf
            for i in range(Nruns):
                RES = results[alg][seeds[i]]
                F[i] = [(f- fmin)/fmin for f in RES["functionalvalue"]]
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
            
            C = [step*i for i in range(len(Fmin))]


            ax.plot(C,Fmin,label=names[alg],color=colors[alg])
            ax.plot(C,Fmax,color=colors[alg])
            ax.fill_between(C,Fmin,Fmax,facecolor=colors[alg], alpha=0.5)
    plt.ylabel("Suboptimality")
    plt.xlabel("Number of scenarios treated")
    plt.legend()
    plt.yscale('log', nonposy='clip')
    plt.savefig(path+"SuboptFill_Calls"+pb+".png")
    matplotlib2tikz.save(path+"SuboptFill_Calls"+pb+".tex")

