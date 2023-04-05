import numpy as np
import math
#import matplotlib.pyplot as plt
#import pylab as pyl
#import textwrap
import sys

#def set_eps( ncolumns = 1, aspect = None, eps=True ):
#    font = {'family' : 'serif',
#            'weight' : 'normal',
#            'size'   : 12 * ncolumns}
#    
#    plt.rc( 'font', **font )
#    plt.rc( 'text', usetex = True )
#    if eps: plt.rc( 'ps', usedistiller = 'xpdf' )
#    
#    plt.rc( 'xtick', labelsize = 10 * ncolumns )
#    plt.rc( 'ytick', labelsize = 10 * ncolumns )#
#
#    fig = plt.figure( figsize = aspect )
#    fsz = fig.get_size_inches()
#    fsz = [fsz[0], fsz[1]]

#    scale = 12.0 / fsz[0]
#    fig.set_size_inches( [scale*fsz[0], scale*fsz[1]] )
#    return fig


#def write_eps( filename, fig, dry_run = False, axes_on = True, eps=True ):
#    if eps:
#        fileext = "eps"
#    else:
#        fileext = "pdf"
#    if not axes_on:
#        ax_main = fig.add_subplot( 111 )
#        ax_main.set_axis_off()
#    if dry_run:
#        pyl.show()
#    else:
#        pyl.savefig( "{}.{}".format( filename, fileext ), transparent = True, bbox_inches='tight' )


#def write_tex( filename, caption = None, label = None, ncolumns = 1, dry_run = False ):
#    if ncolumns == 1:
#        rebel = "*"
#    else:
#        rebel = ""
#    filename2 = "./figures/{}".format( filename.split( "/" )[-1] )
#    tex = textwrap.dedent( r"""
#    \begin{{figure{0}}}[h!]
#      \centerline{{\includegraphics[width={1}\textwidth]{{{2}}}}}
#      \vfill
#      \caption{{
#        {3}
#      }}
#      \label{{{4}}}
#    \end{{figure{0}}}""".format( rebel, 1.00 / ncolumns, filename, caption, label ) )[1:]
#
#    if dry_run:
#        print(tex)
#    else:
#        open( "{}.tex".format( filename ), 'w' ).write( tex )


def parse_data_files(fname,perline=10):
    with open(fname, "r") as f:
        line = f.readline().split()
        nmin = int(line[0])
        entries = int(line[1])
        data = np.zeros((nmin*entries))
        nlines = math.ceil(nmin*entries/perline)
        for idline,line in enumerate(f):
            test_new_min = idline%nlines
            angles = [float(x) for x in line.split()]
            if (test_new_min<nlines):
                data[test_new_min*10:(test_new_min+1)*10] = angles
            else:
                data[test_new_min*10:nmin*entries] = angles
        return np.reshape(data,(nmin,entries))

def histogram_angles(data,weights,nbin,step):
    hist = np.zeros(nbin)
    nmin, nentries = data.shape
    for imin in range(nmin):
        for entry in range(nentries):
            ibin = math.floor(data[imin][entry]/step)
            hist[ibin] += weights[imin]
    return hist
            
    
system = sys.argv[1]
nentries = 10 # number of entries per line
kT = 0.616
nbin = 360
low = 0.0
high = 180.0
plot=False

step = (high-low)/nbin

pep_file = "pitch_pep.dat"
pro_file = "pitch_methylene_pro.dat"
hyp_file = "pitch_methylene_hyp.dat"

pro_data = parse_data_files(pro_file)
hyp_data = parse_data_files(hyp_file)
pep_data = parse_data_files(pep_file)

energies = np.genfromtxt("min.data", dtype=float, usecols=(0,))

energies = energies - np.min(energies)
weights = np.exp(-energies/kT)

#ang_range = list(np.arange(low, high, step, dtype=float))
xrange = np.arange(low+step/2, high+step/2, step, dtype=float)
#ang_range = np.asarray(ang_range+[high])


hist_pep  = histogram_angles(pep_data,weights,nbin,step)
hist_pro  = histogram_angles(pro_data,weights,nbin,step)
hist_hyp  = histogram_angles(hyp_data,weights,nbin,step)

np.save("hist_pep_"+system+".npy", hist_pep)
np.save("hist_pro_"+system+".npy", hist_pro)
np.save("hist_hyp_"+system+".npy", hist_hyp)

weights_dummy = np.zeros(pep_data.shape)
for i in range(pep_data.shape[1]):
    weights_dummy[:,i] = weights
    
print("Peptide - average: ",np.average(pep_data, weights=weights_dummy), " std: ", np.std(pep_data))

weights_dummy = np.zeros(pro_data.shape)
for i in range(pro_data.shape[1]):
    weights_dummy[:,i] = weights
    
print("Proline - average: ",np.average(pro_data, weights=weights_dummy), " std: ", np.std(pro_data))

weights_dummy = np.zeros(hyp_data.shape)
for i in range(hyp_data.shape[1]):
    weights_dummy[:,i] = weights
    
print("Hydroxyproline - average: ",np.average(hyp_data, weights=weights_dummy), " std: ", np.std(hyp_data))

#if plot:
#    nc = 2 # number of columns, controls font size
#    aspect = plt.figaspect(12.0/21.0) # aspect ratio
#    fig1=set_eps(nc , aspect)
#    
#    ax1 = fig1.add_subplot(3,1,1)
#    ax2 = fig1.add_subplot(3,1,2)
#    ax3 = fig1.add_subplot(3,1,3)
#    
#    ax1.bar(xrange,hist_pep/np.sum(hist_pep))
#    ax2.bar(xrange,hist_pro/np.sum(hist_pro))
#    ax3.bar(xrange,hist_hyp/np.sum(hist_hyp))
#    
#    plt.show()
