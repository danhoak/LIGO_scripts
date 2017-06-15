from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pylab
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import cm

matplotlib.rcParams.update({'savefig.dpi':250})

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

F_plus = np.zeros((len(u),len(v)))
F_cross = np.zeros((len(u),len(v)))
F_avg = np.zeros((len(u),len(v)))

x1 = np.zeros((len(u),len(v)))
y1 = np.zeros((len(u),len(v)))
z1 = np.zeros((len(u),len(v)))

x2 = np.zeros((len(u),len(v)))
y2 = np.zeros((len(u),len(v)))
z2 = np.zeros((len(u),len(v)))

x3 = np.zeros((len(u),len(v)))
y3 = np.zeros((len(u),len(v)))
z3 = np.zeros((len(u),len(v)))

for i in range(len(u)):
    for j in range(len(v)):
                    
        phi = u[i]
        theta = v[j]

        F_cross[i,j] = np.abs(-2 * np.cos(phi) * np.sin(2*theta))
        F_plus[i,j] = np.abs((np.cos(theta)**2 * np.cos(phi)**2 - np.sin(phi)**2) - (np.cos(theta)**2 * np.sin(phi)**2 - np.cos(phi)**2))

        #F_avg[i,j] = sqrt(F_plus[i,j]**2 + F_cross[i,j]**2)/2

        x1[i][j] = F_plus[i][j] * np.cos(phi) * np.sin(theta)
        y1[i][j] = F_plus[i][j] * np.sin(phi) * np.sin(theta)
        z1[i][j] = F_plus[i][j] * np.cos(theta)

        x2[i][j] = F_cross[i][j] * np.cos(theta) * np.sin(phi)
        y2[i][j] = F_cross[i][j] * np.sin(theta) * np.sin(phi)
        z2[i][j] = F_cross[i][j] * np.cos(phi)

        R = np.sqrt(F_plus[i,j]**2 + F_cross[i,j]**2)

        x3[i][j] = R * np.cos(theta) * np.sin(phi)
        y3[i][j] = R * np.sin(theta) * np.sin(phi)
        z3[i][j] = R * np.cos(phi)


fig = plt.figure(1)
ax = fig.gca(projection='3d')

plt.hold(True)

#ax.plot([0,2],[0,0],[-2,-2],'k--',linewidth=2,zorder=-1)
#ax.plot([0,0],[0,-2],[-2,-2],'k--',linewidth=2,zorder=-1)

#ax.arrow(0, 0, 2, 0, head_width=0.1, head_length=0.1, fc='k', ec='k',zorder=10)
#ax.arrow(0, 0, 0, -2, head_width=0.1, head_length=0.1, fc='k', ec='k',zorder=10)

ax.plot_surface(x1, y1, z1,  rstride=1, cstride=1, facecolors=cm.jet(F_plus/2),zorder=1)


#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')

ax.set_xlim(-2,2)
ax.set_ylim(-2,2)

"""
plt.tick_params(\
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.tick_params(\
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.tick_params(\
    axis='z',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.tick_params(axis='y',which='both',bottom='off',top='off',labelbottom='off')
"""

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

#pylab.title(r'$\mathtt{F^+(\Omega)}$',fontsize=22)
#pylab.title(r'${F_+(\Omega)}$',fontsize=20)

#ax.set_title('F+')

plt.savefig('F_plus.png',bbox_inches='tight')
#plt.savefig('F_plus.pdf',bbox_inches='tight')


fig = plt.figure(2)
ax = fig.gca(projection='3d')

plt.hold(True)

#ax.plot([0,2],[0,0],'k--',linewidth=2,zorder=0)
#ax.plot([0,0],[0,-2],'k--',linewidth=2,zorder=0)

ax.plot_surface(x2, y2, z2,  rstride=1, cstride=1, facecolors=cm.jet(F_cross/2))

#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')

ax.set_xlim(-2,2)
ax.set_ylim(-2,2)

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

#pylab.title(r'$\mathtt{F^x(\Omega)}$',fontsize=22)
#pylab.title(r'${F_{\times}(\Omega)}$',fontsize=20)

#ax.set_title('Fx')

plt.savefig('F_cross.png',bbox_inches='tight')
#plt.savefig('F_cross.pdf',bbox_inches='tight')





fig = plt.figure(2)
ax = fig.gca(projection='3d')

plt.hold(True)

#ax.plot([0,2],[0,0],'k--',linewidth=2,zorder=0)
#ax.plot([0,0],[0,-2],'k--',linewidth=2,zorder=0)

ax.plot_surface(x3, y3, z3,  rstride=1, cstride=1, facecolors=cm.jet(F_cross/2))

#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')

ax.set_xlim(-2,2)
ax.set_ylim(-2,2)

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

#pylab.title(r'$\mathtt{F^x(\Omega)}$',fontsize=22)
#pylab.title(r'${F_+^2 + F_{\times}^2}$',fontsize=20)

#ax.set_title('Fx')

plt.savefig('F_avg.png',bbox_inches='tight')
#plt.savefig('F_cross.pdf',bbox_inches='tight')
