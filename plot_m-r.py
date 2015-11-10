import numpy as np
import matplotlib.pyplot as plt

f = open('mr_data.dat','r')
d = {}
h = f.readline().strip().split()
for k in h:
    d[k] = []
for l in f:
    ls = l.strip().split()
    for k,v in zip(h,ls):
        d[k].append(float(v.replace('d','e')))
f.close()

for k,v in d.iteritems():
    d[k] = np.array(v)

plt.figure()
fig = plt.gcf()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(d['radius']/1.0e5,d['mass_msun'],color='black')

ax.annotate('$\\mathrm{\\rho_c = 18 \\times 10^{14} ~ g/cm^3}$', xy=(13.7,3.85),
           xytext=(13.2,3.5), arrowprops=dict(facecolor='black'))
ax.annotate('$\\mathrm{\\rho_c = 3 \\times 10^{14} ~ g/cm^3}$', xy=(13.42,0.8),
           xytext=(13.42,1.3), arrowprops=dict(facecolor='black'))
ax.annotate('$\\mathrm{\\rho_c = 5 \\times 10^{14} ~ g/cm^3}$', xy=(14.82,1.67),
           xytext=(14.5,0.7), arrowprops=dict(facecolor='black'))

x_ax_range = np.linspace(13.0,15.5,num=2,endpoint=True)

yup = 2.74*np.ones_like(x_ax_range)
ydn = 1.00*np.ones_like(x_ax_range)

ax.plot(x_ax_range,ydn,x_ax_range,yup,color='blue',alpha=0.2)
ax.fill_between(x_ax_range,ydn,yup,where=yup>ydn,facecolor='blue',alpha=0.2)

plt.text(13.15,2.25,'Observational NS mass range')
plt.text(13.15,2.00,'(cf. Lattimer, ARNPP, Vol. 62, 485, 2012)')

plt.xlabel('Radius (km)')
plt.ylabel('Mass $\\mathrm{(M_\\odot)}$')
#plt.tight_layout()
plt.savefig('m-r.pdf')

plt.clf()

plt.figure()
fig = plt.gcf()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(d['rho0']/1.0e14,d['mass_msun'],color='black')

x_ax_range = np.linspace(2,18,num=2,endpoint=True)

yup = 2.74*np.ones_like(x_ax_range)
ydn = 1.00*np.ones_like(x_ax_range)

ax.plot(x_ax_range,ydn,x_ax_range,yup,color='blue',alpha=0.2)
ax.fill_between(x_ax_range,ydn,yup,where=yup>ydn,facecolor='blue',alpha=0.2)

plt.text(7,2.25,'Observational NS mass range')
plt.text(7,2.00,'(cf. Lattimer, ARNPP, Vol. 62, 485, 2012)')

plt.xlabel('Central Density $\\mathrm{(\\times 10^{14}~g/cm^3)}$')
plt.ylabel('Mass $\\mathrm{(M_\\odot)}$')
#plt.tight_layout()
plt.savefig('m-rhoc.pdf')
