import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
levels = {
 "Patient level (12 studies)": {"nodes":{"AI":11,"Human":11,"Collab":2},
    "edges":{("AI","Human"):11,("AI","Collab"):2,("Human","Collab"):2}},
 "Eye level (29 studies)": {"nodes":{"AI":21,"Human":24,"Collab":8},
    "edges":{("AI","Human"):21,("AI","Collab"):5,("Human","Collab"):8}},
}
pos={"AI":(0,1),"Human":(-0.95,-0.62),"Collab":(0.95,-0.62)}
full={"AI":"Autonomous AI","Human":"Human readers","Collab":"Human–AI collaboration"}
colors={"AI":"#2166AC","Human":"#D6604D","Collab":"#1A9641"}
lblpos={"AI":(0,1.34,"center"),"Human":(-0.95,-0.95,"center"),"Collab":(0.95,-0.95,"center")}
fig,axes=plt.subplots(1,2,figsize=(11,5.4))
for ax,(title,g) in zip(axes,levels.items()):
    for (a,b),w in g["edges"].items():
        x=[pos[a][0],pos[b][0]]; y=[pos[a][1],pos[b][1]]
        ax.plot(x,y,'-',lw=1.5+w*0.7,color="#9aa0a6",zorder=1,solid_capstyle="round")
        mx,my=(x[0]+x[1])/2,(y[0]+y[1])/2
        ax.text(mx,my,str(w),fontsize=12,fontweight="bold",ha="center",va="center",
                bbox=dict(boxstyle="round,pad=0.28",fc="white",ec="#9aa0a6",lw=1.1),zorder=3)
    for n,(x,y) in pos.items():
        val=g["nodes"][n]
        ax.scatter([x],[y],s=2000+val*130,c=colors[n],zorder=2,edgecolors="white",linewidths=2)
        ax.text(x,y,n,fontsize=12,fontweight="bold",ha="center",va="center",color="white",zorder=4)
        lx,ly,ha=lblpos[n]
        ax.text(lx,ly,full[n],fontsize=10,ha=ha,va="center",color=colors[n],fontweight="bold")
    ax.set_title(title,fontsize=13,fontweight="bold",pad=10)
    ax.set_xlim(-1.75,1.8); ax.set_ylim(-1.4,1.65); ax.axis("off")
fig.suptitle("Network geometry: studies providing direct (within-study) head-to-head evidence per comparison",
             fontsize=12.5,y=1.0,fontweight="bold")
fig.text(0.5,0.015,"Edge label = number of studies contributing direct evidence for that comparison; node size reflects studies contributing each strategy. "
         "Every comparison carries within-study direct evidence (no purely indirect edges).",
         ha="center",fontsize=8.5,style="italic")
plt.tight_layout(rect=[0,0.045,1,0.95])
plt.savefig("results/Figure_Network_Geometry.png",dpi=300,bbox_inches="tight")
plt.savefig("results/Figure_Network_Geometry.svg",bbox_inches="tight")
print("saved")
