import os
import math
import numpy as np
import subprocess
from collections import OrderedDict

pars = OrderedDict([])
pars['poly_K'] = 1.982e-6
pars['poly_gamma'] = 2.75
pars['rho0'] = 10e14
pars['DR_SAVE'] = 0.01e5
pars['NDR_SAVE'] = 3000

pars_t = OrderedDict([])
pars_t['poly_K'] = 'd'
pars_t['poly_gamma'] = 'd'
pars_t['rho0'] = 'd'
pars_t['DR_SAVE'] = 'd'
pars_t['NDR_SAVE'] = 'i'

def float_to_string(flvar):
    return '{0:1.14e}'.format(flvar).replace('e','d')

def run_pars(pars,pars_t,klabel):
    f = open('hse.par','w')
    for k,v in pars.iteritems():
        if pars_t[k] == 'i':
            f.write(k + '          ' + str(v) + '\n')
        elif pars_t[k] == 'd':
            f.write(k + '          ' + float_to_string(v) + '\n')
    f.close()
    subprocess.call('./inthse', shell=True)

    f = open('hse_profile.dat','r')
    lfin = ''
    for l in f:
        lfin = l
    lsfin = lfin.split()

    data = {'mass':float(lsfin[0]), 'pres':float(lsfin[1]), 'radius':float(lsfin[-1])} 
    subprocess.call('mv hse.par hse_R-'+klabel+'.par', shell=True)
    subprocess.call('mv hse_profile.dat hse_profile_R-'+klabel+'.dat', shell=True)
    return data

j = 0
N = 1000
q = int(math.ceil(np.log10(N)))

def knum_to_str(kn):
    s = '{0:0>' + str(q) + '}'
    return s.format(kn)

GramsPerMsun = 1.988435e33

mr_data = OrderedDict([('mass',np.zeros(N)), ('pres',np.zeros(N)), ('radius', np.zeros(N)), ('rho0', np.zeros(N))])

# To reproduce approximately the observed NS masses, use:
# rho_lower = 3.0e14
# rho_upper = 9.0e14

for rho in np.linspace(3.0e14, 18.0e14, num=N, endpoint=True):
    pars['rho0'] = rho
    d = run_pars(pars,pars_t,knum_to_str(j))
    mr_data['mass'][j] = d['mass']
    mr_data['pres'][j] = d['pres']
    mr_data['radius'][j] = d['radius']
    mr_data['rho0'][j] = rho

    if (abs(d['pres']) > 1.0):
        print 'Error: pressure nonzero for j='+str(j)
        print 'Tried rho0 = ' + float_to_string(rho)
        exit()
    
    j += 1

mr_data['mass_msun'] = mr_data['mass']/GramsPerMsun

f = open('mr_data.dat','w')
ws = '       '
f.write('j' + 3*ws)
for k in mr_data.keys():
    f.write(k + 3*ws)
f.write('\n')
for i in xrange(0,N):
    f.write(str(i) + ws)
    for k,v in mr_data.iteritems():
        f.write(float_to_string(v[i]) + ws)
    f.write('\n')
f.close()
    
