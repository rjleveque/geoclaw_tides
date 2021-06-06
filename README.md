# geoclaw_tides

Some examples of forcing tidal variation in a GeoClaw simulation.

So far there is only one example `GraysHarbor` with a test of the tidal 
forcing alone, by comparing tide gauges in Aberdeen and Westport.

After cloning this repository, open `GraysHarbor/README.html` for 
instructions and links to some sample results.

To be added: An example modeling a tsunami on top of this tide. And other 
examples...

Developed using [GeoClaw](http://www.geoclaw.org) from Clawpack Version 5.8.0.
Most things should work fine with other versions too, but note that 
the specification of topofiles in setrun.py changed in this version, see
the [release notes](http://www.clawpack.org/release_5_8_0.html).
