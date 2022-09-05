import matplotlib.pyplot as plt
import seaborn as sns

def ReadMatrix(B,data,start, size):
    for i in range(start,start+size):
        B.append(list(map(int, data[i].split())))

with open("GSA_d.txt") as file:
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]

[n,accuracy_factor,blocksize,N] = list(map(int, lines[0].split()))
c = float(lines[1])

B,B_red,U,sz = [[],[],[],[]]

t = 2

ReadMatrix(B,lines,t,n)
t += n +2 

ReadMatrix(B_red,lines,t,n)
t += n + 2

ReadMatrix(U,lines,t,n)
for l in U:
    for x in range(len(l)):
        l[x] = abs(l[x])    
t += n + 2

sz = list(map(float, lines[t].split()))

ratio = [(sz[i]/sz[i-1])**2 for i in range(1,len(sz))]

fig, axs = plt.subplots(nrows = 2,ncols=2)
axs[0][0].set_box_aspect(1)
axs[1][0].set_box_aspect(1)
axs[0][1].set_box_aspect(1)
fig.suptitle("[n,accuracy_factor,blocksize,N, c] = [" + str(n) + " " +  str(accuracy_factor) + " "+ str(blocksize) + " "+ str(N) + " "+ str(c) + " ]"  )
sns.heatmap(B, ax = axs[0][0]).set(title = "Basis")
sns.heatmap(B_red, ax = axs[0][1]).set(title = "Reduced basis")
sns.heatmap(U, ax = axs[1][0]).set(title = "Transformation Matrix")
sns.lineplot(range(len(ratio)),ratio, ax = axs[1][1]).set(title = "GSA test")
plt.tight_layout()
plt.savefig(str(n) + "_" +  str(accuracy_factor) + "_"+ str(blocksize) + "_"+ str(N) + "_"+ str(c) + ".png", dpi = 200)
# plt.show()