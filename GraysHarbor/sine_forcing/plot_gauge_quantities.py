
from pylab import *
import clawpack.pyclaw.gauges as gauges

figure(400, figsize=(10,6))
clf()
colors = ['k','b','r','m','g']

outdir = '_output'

for k,gaugeno in enumerate([0,1102,1187]):
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    t = gauge.t / 3600.   # convert to hours
    q = gauge.q
    eta = q[3,:]
    plot(t, eta, colors[k], label='Gauge %s' % gaugeno)

    m2 = int(floor(0.75*len(eta)))
    eta2 = eta[m2:]  # last part of eta signal
    etamax2 = eta2.max()
    etamin2 = eta2.min()
    t2 = t[m2:]
    jtmax = argmax(eta2) 
    tshift = (t2[jtmax] - 39.)*3600.
    
    print('At gauge %i, etamin2 = %.3f, etamax2 = %.3f at tshift = %.1f s' \
            % (gaugeno,etamin2,etamax2,tshift))

# forcing:
tperiod = 12
eta = 1.*sin(2*pi*t/tperiod)
plot(t, eta, 'k--', label='Ocean forcing')

#plot(t, 1.068*ones(t.shape), 'b--', label='MHW 1102')
#plot(t, -1.068*ones(t.shape), 'b--')
#plot(t, 1.21*ones(t.shape), 'r--', label='MHW 1187')
#plot(t, -1.21*ones(t.shape), 'r--')

legend()
xlabel('hours')
ylabel('Surface relative to MTL (m)')
grid(True)
title('Sine wave forcing in ocean and resulting GeoClaw gauge results')

if 1:
    xticks(arange(0,t[-1]+0.1,12))
    xlim(0,t[-1]+0.1)

if 1:
    fname = 'GaugeComparison.png'
    savefig(fname, bbox_inches='tight')
    print('Created %s' % fname)

