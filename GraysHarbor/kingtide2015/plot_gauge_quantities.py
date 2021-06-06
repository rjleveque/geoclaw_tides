
from pylab import *
import clawpack.pyclaw.gauges as gauges
from noaa_gauges import fetch_gauges

figure(400, figsize=(13,5))
clf()
colors = ['b','r','m','g']

outdir = '_output'
#outdir = '_gauges_mtl_rivers_4levels'

t_ocean,etap_ocean,t_westport,etap_westport,t_aberdeen,etap_aberdeen = \
        fetch_gauges()

plot(t_ocean/3600., etap_ocean, color=[.5,.5,.5],
     linestyle='--', label='Ocean forcing')

for k,gaugeno in enumerate([1102,1187]):
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    t = gauge.t / 3600.   # convert to hours
    q = gauge.q
    eta = q[3,:]
    if gaugeno==1102:
        plot(t, eta, colors[k], label='GeoClaw 1102 - Westport')
    elif gaugeno==1187:
        plot(t, eta, colors[k], label='GeoClaw 1187 - Aberdeen')
    else:
        plot(t, eta, colors[k], label='GeoClaw Gauge %s' % gaugeno)

# forcing:
if 0:
    tperiod = 12
    #tperiod = 12 * 3600.
    eta = 1.*sin(2*pi*t/tperiod)
    plot(t, eta, 'k--', label='Ocean forcing')

plot(t_westport/3600., etap_westport, 'c', label='NOAA 1102 - Westport')
plot(t_aberdeen/3600., etap_aberdeen, 'm', label='NOAA 1187 - Aberdeen')

xlim(0,96)
ylim(-2.5,2.5)
legend(loc='upper left', fontsize=8)
xlabel('hours')
ylabel('Surface relative to MTL (m)')
title('Grays Harbor tides, Dec. 21-26, 2015')
grid(True)

xlim(0,t[-1])
if 1:
    # ticks every 12 hours:
    xticks(arange(0,t[-1]+0.1,12))


if 1:
    fname = 'GaugeComparison.png'
    savefig(fname, bbox_inches='tight')
    print('Created %s' % fname)

