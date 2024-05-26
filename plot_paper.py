''' After having run 'proof_size_estimate.py', this python script computes the corresponding plots for Section 6.
To obtain Figure 3 (left) in submission, run the python script with the only function not commented out being 'plot_512_1024_AS_lin()'
To obtain Figure 3 (right) in submission, run the python script with the only function not commented out being 'plot_512_AS_lin_no_salt()'
To obtain Figure 4 in submission, first run the python script with the only function not commented out being 'plot_512_AS_2S_FS_lin()', then run it again with the only function not commented out being 'plot_1024_AS_2S_FS_lin()'.
'''

from pandas import *
import matplotlib.pyplot as plt
from matplotlib.ticker import  AutoMinorLocator
plt.style.use('tableau-colorblind10')

def kB(sizes_bits):
    conversion = lambda b: b / 8000
    sizes_kB = list(map(conversion, sizes_bits))
    return sizes_kB

def add_salts(num_sigs, sizes):
    add = lambda e: e[0] + (320 * e[1])
    return list(map(add, zip(sizes, num_sigs)))

#data = read_csv("estimates-lin-zoom.csv")
data = read_csv("estimates-lin-.csv")
# print(data.keys)
#with AS
heading ="Num-Sigs,Naive-512,Naive-1024,Falcon-512-2S,Falcon-512-AS,Falcon-1024-2S,Falcon-1024-AS"

num_sigs = data["Num-Sigs"].tolist()
naive_512 = data["Naive-512"].tolist()
naive_1024 = data["Naive-1024"].tolist()
f512_2S = data["Falcon-512-2S"].tolist()
f512_AS = data["Falcon-512-AS"].tolist()
f1024_2S = data["Falcon-1024-2S"].tolist()
f1024_AS = data["Falcon-1024-AS"].tolist()

data_mod = read_csv("estimates-mod-size.csv")
num_sigs_mod = data_mod["Num-Sigs"].tolist()
mod_512 = data_mod["Falcon-512-log-q"].tolist()
mod_1024 = data_mod["Falcon-1024-log-q"].tolist()

''' Naive signature sizes for optimized Falcon as described in [ETWY22]
 Numbers for skew factor gamma = 8 come from Table 1 in ETWY22, where we translate from bytes to bits; 
 The salt is still of size 40 bytes = 320 bits
'''

elips_512 = [8 * (410-40)* n for n in num_sigs] 
elips_1024= [8 * (780-40)* n for n in num_sigs]

fig, ax = plt.subplots()

def plot_512_1024_2S_lin(): 
    plt.title("Comparison With Trivial Aggregation")
    ax.plot(num_sigs, kB(add_salts(num_sigs, f512_2S)),  label="Our Aggregation of Falcon-512")
    ax.plot(num_sigs, kB(add_salts(num_sigs, f1024_2S)),  label="Our Aggregation of Falcon-1024", linestyle="dashed")
    plt.plot(num_sigs, kB(add_salts(num_sigs, naive_512)),  label="Trivial Aggregation of Falcon-512")
    plt.plot(num_sigs, kB(add_salts(num_sigs, naive_1024)),  label="Trivial Aggregation of Falcon-1024", linestyle="dashed")
    plt.ylabel("Size (in kB)")       
    plt.xlabel("Number of Signatures")
    ax.legend()
    ax.grid(axis="y")
    ratio = 0.5
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)
    plt.savefig("plot-512-1024-2S-lin.pdf")

def plot_512_1024_2S_lin_zoom(): 
    plt.title("Comparison With Trivial Aggregation")
    ax.plot(num_sigs, kB(add_salts(num_sigs, f512_2S)),  label="Our Aggregation of Falcon-512")
    ax.plot(num_sigs, kB(add_salts(num_sigs, f1024_2S)),  label="Our Aggregation of Falcon-1024", linestyle="dashed")
    plt.plot(num_sigs, kB(add_salts(num_sigs, naive_512)),  label="Trivial Aggregation of Falcon-512")
    plt.plot(num_sigs, kB(add_salts(num_sigs, naive_1024)),  label="Trivial Aggregation of Falcon-1024", linestyle="dashed")
    plt.ylabel("Size (in kB)")       
    plt.xlabel("Number of Signatures")
    ax.legend()
    ax.grid(axis="y")
    ratio = 0.5
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)
    plt.savefig("plot-512-1024-2S-lin_zoom.pdf")

def plot_512_2S_lin_no_salt(): 
    plt.title("Effect of Salt on Aggregate Signature")
    ax.plot(num_sigs, kB(add_salts(num_sigs, f512_2S)),  label="Our Aggregation of Falcon-512")
    ax.plot(num_sigs, kB(f512_2S),  label="Our Aggregation of Falcon-512 (No Salt)", linestyle="dashed")
    plt.plot(num_sigs, kB(add_salts(num_sigs, naive_512)),  label="Trivial Aggregation of Falcon-512")
    plt.plot(num_sigs, kB(naive_512),  label="Trivial Aggregation of Falcon-512 (No Salt)", linestyle="dashed")
    plt.ylabel("Size (in kB)")       
    plt.xlabel("Number of Signatures")
    ax.legend()
    ax.grid(axis="y")
    ratio = 0.5
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)
    plt.savefig("plot-512-2S-lin-no-salt.pdf")
    
def plot_512_AS_2S_lin(): 
    plt.title("Comparison With Other Challenge Sets (Falcon-512)")
    plt.plot(num_sigs, kB(add_salts(num_sigs, f512_AS)),  label="Almost-Fully-Splitting for Falcon-512")
    plt.plot(num_sigs, kB(add_salts(num_sigs, f512_2S)),  label="Two-Splitting for Falcon-512", linestyle="dashed")
    plt.ylabel("Size (in kB)")       
    plt.xlabel("Number of Signatures")
    ax.legend()
    ax.grid(axis="y")
    ratio = 0.5
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)
    plt.savefig("plot-512-AS-2S-lin.pdf")
    
def plot_1024_AS_2S_lin(): 
    plt.title("Comparison With Other Challenge Sets (Falcon-1024)")
    ax.plot(num_sigs, kB(add_salts(num_sigs, f1024_AS)),  label="Almost-Fully-Splitting for Falcon-1024")
    ax.plot(num_sigs, kB(add_salts(num_sigs, f1024_2S)),  label="Two-Splitting for Falcon-1024", linestyle="dashed")
    plt.ylabel("Size (in kB)")       
    plt.xlabel("Number of Signatures")
    ax.legend()
    ax.grid(axis="y")
    ratio = 0.5
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)
    plt.savefig("plot-1024-AS-2S-lin.pdf")
  

### Figure 3 (left) in submission
plot_512_1024_AS_lin() 

### zoom
#plot_512_1024_AS_lin_zoom() 

### Figure 3 (right) in submission
#plot_512_2S_lin_no_salt()

### Figure 4 (left) in submission
#plot_512_AS_2S_lin()

### Figure 4 (right) in submission
#plot_1024_AS_2S_lin()


