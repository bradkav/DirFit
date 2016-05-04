import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
from scipy.stats import chi2, norm
from scipy.integrate import quad
from numpy.polynomial.chebyshev import chebval, chebfit
from scipy.special import iv,ive
import ternary
font = {'family' : 'sans-serif',
        'size'   : 16}

mpl.rc('font', **font)

#from scipy.interpolate import griddata
#from scipy.spatial import *

from numpy import pi

N = 3
N_terms = 4

chaindir = "chains/"
plotlabel = " "

data = np.loadtxt(chaindir+ "data.txt")


#------Options------
plot_ternary = 0
plot_f = 1
plot_2d = 0

#----Calculating 1d profile likelihoods and plots----
def profile_1d(data, likes, minval, maxval):
    nbin = 250
    proflike = np.zeros(nbin)
    binedges = np.linspace(minval, maxval, nbin+1)
    
    minL = np.min(likes)
    #Multinest spits out -2*lnL
    L = 0.5*(likes - minL)

    #Divide into bins
    Nsamps = likes.size    
    
    data1d = data
    
    for i in range(nbin):
        Li = L[np.where(data1d > binedges[i])]
        data1di = data1d[np.where(data1d > binedges[i])]
        Lii = Li[np.where(data1di < binedges[i+1])]
        
        if (Lii.size < 1):
            proflike[i] = 1e30
        else:
            proflike[i] = np.min(Lii)
        
    return binedges, proflike
     
def plot_profile_1d(data, likes, minval, maxval, xlabel='Parameter value'):
    bins, PL = profile_1d(data, likes, minval, maxval)

    pl.figure()
    pl.plot(bins[1:] - 0.5*(bins[1] - bins[0]), PL, 'k-', linewidth=2.0)
    pl.ylim([0,10])

    pl.axhline(0.5*chi2.ppf(0.683, 1), color='b', linestyle='--')
    pl.axhline(0.5*chi2.ppf(0.995, 1), color='b', linestyle='--')

    pl.text(2.7, 0.5*chi2.ppf(0.683, 1), r'$1\sigma$')
    pl.text(2.7, 0.5*chi2.ppf(0.995, 1), r'$2\sigma$')

    pl.xlabel(xlabel)
    pl.ylabel(r'$\Delta log \mathcal{L}$')
     
     
def param_intervals_1d(data, likes, paramstr='Parameter value'):  
    #Note that 'likes' is actually -2*lnL  
    cuts = data[np.where(0.5*(likes - np.min(likes)) < 0.5*chi2.ppf(0.683, 1))]
    plus = np.max(cuts) - data[np.argmin(likes)]
    minus = data[np.argmin(likes)] - np.min(cuts)
     
    #Print best fit...
    print "   " + paramstr + "    = (", '{:.3f}'.format(data[np.argmin(likes)]), " + ", '{:.3f}'.format(plus), " - ", '{:.3f}'.format(minus), ")"
     
#----Calculating 2d profile likelihood and plots----
def profile_2d(data, dim, minval, maxval, nbin=100):
    proflike = np.zeros((nbin,nbin))
    
    minL = np.min(data[:,1])
    L = 0.5*(data[:,1] - minL)
    
    binedges0 = np.linspace(minval[0], maxval[0], nbin+1)
    binedges1 = np.linspace(minval[1], maxval[1], nbin+1)
    
    for i in range(nbin):
        for j in range(nbin):
            inds = np.where((data[:, 2+dim[0]] > binedges0[i]) & (data[:, 2+dim[0]] < binedges0[i+1]) & (data[:, 2+dim[1]] > binedges1[j]) & (data[:, 2+dim[1]] < binedges1[j+1]))
            if (len(inds[0]) < 1):
                proflike[i,j] = 1e30
            else:
                proflike[i, j] = np.min(L[inds])
    
    points0 = binedges0[1:] - 0.5*(binedges0[1] - binedges0[0])
    points1 = binedges1[1:] - 0.5*(binedges1[1] - binedges1[0])
    
    return np.vstack((points0,points1)), proflike   


#---Velocity distribution stuff---

def f_poly(v, a):
    return np.exp(-chebval(2*v/1000.0 - 1, a))

#---------------------------------------------------
#Calculate the discretised 'f's
def calcf(v,j, v_lag=220, v_rms=156, f_s=0.167):
    #print v_lag, v_rms
    Norm = ((v**(-1))/(v_rms*v_lag*np.sqrt(2*pi)))
    R = v*v_lag/(v_rms**2)
    frac = (np.cos((j-1)*pi*1.0/N) - np.cos((j)*pi*1.0/N))
    d = ((v**2+v_lag**2)/(2*v_rms**2))
    #A = np.exp(-((v**2+v_lag**2)/(2*v_rms**2)))
    B = np.exp(R*np.cos(((j)-1)*pi*1.0/N) - d) - np.exp(R*np.cos((j)*pi*1.0/N) - d)
    #print -((v**2+v_lag**2)/(2*v_rms**2))
    #print R*np.cos(((j)-1)*pi*1.0/N)
    #print R*np.cos((j)*pi*1.0/N)
    #print "   "
    
    
    return (Norm/frac)*B*(1-f_s) + f_s*calcf2(v,j, v_lag=400, v_rms=10)

#Calculate the discretised 'f's
def calcf2(v,j, v_lag=400, v_rms=10):
    if (j != 2):
        return 0
    Norm = ((v**(-1))/(v_rms*v_lag*np.sqrt(2*pi)))
    R = v*v_lag/(v_rms**2)
    frac = (np.cos((j-1)*pi*1.0/N) - np.cos((j)*pi*1.0/N))
    d = ((v**2+v_lag**2)/(2*v_rms**2))
    
    #
    B = np.exp(+R - d) - np.exp(-R - d)
    #print -((v**2+v_lag**2)/(2*v_rms**2))
    #print R*np.cos(((j)-1)*pi*1.0/N)
    #print R*np.cos((j)*pi*1.0/N)
    #print "   "
    return (Norm/frac)*B

#--------
def generate_ternary_heatmap(PL, scale=1.0):
    from ternary.helpers import simplex_iterator
    d = dict()
    for (i,j,k) in simplex_iterator(scale):
        d[(i,j)] = PL[i-1,j-1]
    return d




#---Actual calculations---


print "------Parameter intervals--------"
param_intervals_1d(10**data[:,2], data[:,1], paramstr='m_x / GeV')
param_intervals_1d(10**(data[:,3]+39), data[:,1], paramstr='sigma_SD / 1e-39 cm^2')
#param_intervals_1d(data[:,4], data[:,1], paramstr='v_lag / km s^-1')
#param_intervals_1d(data[:,5], data[:,1], paramstr='sig_v / km s^-1')
print "---------------------------------"

#print chi2.ppf([0.683, 0.955, 0.99999943], 2)

#Let's do a little fit

c = np.zeros((N, N_terms))
vvals = np.linspace(1,1000, 100)
avals = 2.0*vvals/1000.0 - 1.0 
"""
for k in range(N):
    c[k,:] = chebfit(avals, -np.log(calcf(vvals, k+1)), N_terms-1)
#print c
for k in range(N):
    #print -np.sum(c[k,1:]*((-1.0)**np.arange(1,N_terms))) - np.log(3.5e-8)
    c[k,0] = -np.sum(c[k,1:]*((-1.0)**np.arange(1,N_terms))) - np.log(3.5e-8)
"""
livepts = np.loadtxt(chaindir + "dataphys_live.points")

bf = (livepts[np.argmax(livepts[:,-2]),:])
print "Max log-like:", np.max((livepts[:,-2]))
print "Best fit point:"
print bf[:-2]

angnorms = np.zeros(N)
for k in range(N):
    f = lambda v: v*v*calcf(v, k+1)
    angnorms[k] = 100*(np.cos(pi*(k)/N) - np.cos(pi*(k+1.0)/N))*quad(f, 0, 1000)[0]



if (plot_ternary):
    points = []
    points.append((angnorms[0], angnorms[1], angnorms[2]))
    print points

    #Ternary plot
    scale = 100
    fontsize = 14.0
    figure, tax = ternary.figure(scale=scale)

    #tax.scatter(points, color='Green', marker='^') 
    bins_2d, PL = profile_2d(data, [11,12], [0.0, 0.0], [1.0, 1.0], nbin=scale)
    #print PL
    d = generate_ternary_heatmap(np.clip(PL, 0, 50), scale=scale)
    tax.heatmap(d, style="h", cmap='Blues_r')
    tax.boundary(linewidth=2.0)
    tax.right_parallel_line(angnorms[0], linewidth=2., color='red', linestyle="--")
    tax.left_parallel_line(angnorms[1], linewidth=2., color='red', linestyle="--")
    tax.horizontal_line(angnorms[2], linewidth=2., color='red', linestyle="--")


    tax.left_axis_label(r"$\mathcal{N}_3$", fontsize=fontsize)
    tax.right_axis_label(r"$\mathcal{N}_2$", fontsize=fontsize)
    tax.bottom_axis_label(r"$\mathcal{N}_1$", fontsize=fontsize)


    tax.ticks(axis='lbr', linewidth=1, multiple=10)
    tax.clear_matplotlib_ticks()


    tax.savefig(chaindir + 'Ternary.pdf')
    tax.show()
    



#print bf

if (plot_f):
    #print np.arange(0, (N_terms-1)*N, (N_terms-1))
    #print bf
    #Add in the necessary zeros
    vp = np.insert(bf[2:(-2-N)], np.arange(0, (N_terms-1)*N, (N_terms-1)), 0)
    #vp = np.append(0, bf[2:-2])
    vp = vp.reshape(N,N_terms)
    #vp = c
    #vp[:,0] = 0

    
    
    #They seem to be inverted (k = 2, 3)
    for k in range(N):
        vp[k,0] = -np.sum(vp[k,1:]*((-1.0)**np.arange(1,N_terms)))
    
    

    Norm = 0
    for k in range(N):
        f = lambda x: x*x*f_poly(x, vp[k,:])
        Norm += (np.cos(pi*(k)/N) - np.cos(pi*(k+1.0)/N))*quad(f, 0, 1000)[0]

    for k in range(N):
        vp[k,0] += np.log(Norm)

    #print vp

    vlist = np.linspace(0,1000, 100)

    #print vp[2,:]

    fig, axarr = pl.subplots(N, sharex=True, figsize=(6,8)) 
   
    for k in range(N):
        #(10**(39+bf[1]))*
        axarr[k].plot(vlist, f_poly(vlist, c[k,:]), 'g:', linewidth=2.0)
        axarr[k].plot(vlist, f_poly(vlist, vp[k,:]), 'r--', linewidth=2.0)
        axarr[k].plot(vlist, calcf(vlist,k+1), 'b-', linewidth=2.0)
        axarr[k].set_yscale("log")
        #axarr[k].text("k = " + str(k+1))
        axarr[k].annotate("k = " + str(k+1), xy=(0.85, 0.85), xycoords='axes fraction')
        axarr[k].set_ylim(1e-13, 2e-7)
        axarr[k].xaxis.grid(True)
        axarr[k].yaxis.grid(True)

    axarr[0].set_title(plotlabel)
    axarr[1].set_ylabel(r'$f^{k}(v)$ km$^{-1}$ s')
    axarr[2].set_xlabel(r'$v$ / km s$^{-1}$')
    pl.subplots_adjust(hspace=None)
    pl.tight_layout()

    pl.savefig(chaindir + 'veldist.pdf')

    pl.show()

#Plot 1-D profile likelihood for mass
plot_profile_1d(10**data[:,2], data[:,1], 10.0, 1000.0, xlabel=r'$m_\chi$ / GeV')
#Plot 1-D profile likelihood for cross-section
plot_profile_1d(10**data[:,3], data[:,1], 1e-40, 1e-38, xlabel=r'$\sigma_{SD}$ / cm$^2$')
#pl.show()



livepts = np.loadtxt(chaindir + "dataphys_live.points")

fig = pl.figure(figsize=(16,7))
ax1 = fig.add_subplot(1,2,1)
#np.min(livepts[:,-2])
#chisq = 2.0*(-np.min(livepts[:,-2]) + livepts[:,-2])

sigma = 2
dof = 2
lim = chi2.ppf(2*norm.cdf(sigma)-1,dof)
#print lim

#The livepoints correspond to lnL
chisq = 2*(livepts[:,-2] - np.max(livepts[:,-2]))
valrange = np.where(chisq  > -1000)
pl.scatter(livepts[valrange,0], livepts[valrange,1], alpha=0.25, c=chisq[valrange])
pl.plot(livepts[np.argmax(chisq),0], livepts[np.argmax(chisq),1], 'g^', alpha=1.0, markersize=12)
pl.colorbar()
pl.axvline(np.log10(50.0), color='r', linestyle='--', linewidth=1.5)
pl.axhline(-39.0, color='r', linestyle='--', linewidth=1.5)
pl.xlim(1.0,3.0)
pl.ylim(-40.0, -38.0)
pl.title('Live samples')

ax2 = fig.add_subplot(1,2,2)
#np.min(data[:,1])
pl.scatter(data[:,2], data[:,3], alpha=0.25, c=0.5*(-data[:,1]))
pl.colorbar()
pl.axvline(np.log10(50.0), color='r', linestyle='--', linewidth=1.5)
pl.axhline(-39.0, color='r', linestyle='--', linewidth=1.5)
pl.xlim(1.0,3.0)
pl.ylim(-40.0, -38.0)
pl.title('Discarded samples')
pl.plot(data[np.argmin(data[:,1]),2], data[np.argmin(data[:,1]),3], 'g^', alpha=1.0, markersize=12)



#Plot 2-D contours for mass and cross-section
if (plot_2d):
    fig = pl.figure()
    ax = fig.add_subplot(1,1,1)

    lims = 2.0*chi2.ppf([0.683, 0.955, 0.99999943], 2)
    #print lims
    bins2d, PL2d = profile_2d(data, [0,1], [1,-40], [3,-37])

    CF = ax.contour(10**bins2d[0,:], 10**bins2d[1,:], PL2d.T, lims, linewidths=2.0, colors='k', alpha=0.7)
    #ax.plot(10**data[:,2], 10**data[:,3], 'k.', alpha=0.1)
    ax.plot(10**data[np.argmin(data[:,1]),2], 10**data[np.argmin(data[:,1]),3], 'g^', alpha=1.0)

    ax.axvline(50, color='r', linestyle='--')
    ax.axhline(1e-39, color='r', linestyle='--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([10, 100, 1000])
    ax.set_yticks([1e-40, 1e-39, 1e-38, 1e-37])
    fmt = {}
    strs = [r'$1\sigma$',r'$2\sigma$', r'$5\sigma$']
    for l,s in zip( CF.levels, strs ):
        fmt[l] = s


    #ax.plot(10**data[np.argmin(L),2], 10**data[np.argmin(L),3], 'g^', markersize=10)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.set_ylim([1e-40, 1e-37])

    # Label every other level using strings
    pl.clabel(CF,CF.levels,inline=True,fmt=fmt,fontsize=12)

    pl.xlabel(r'$m_\chi$ / GeV', size=18)
    pl.ylabel(r'$\sigma_p^{SD}$ / cm$^2$', size=18)
    pl.annotate(plotlabel, xy=(0.05, 0.95), xycoords='axes fraction')
    #pl.xlabel(r'$v_\mathrm{lag}$ / km s$^{-1}$', size=18)
    #pl.xlabel(r'$\sigma_v$ / km s$^{-1}$', size=18)
    pl.tight_layout()

    pl.savefig(chaindir + 'mx-sigma.pdf')

pl.show()