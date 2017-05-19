#!/usr/bin/env python
# coding: utf8

from pyrocko import gf, util, trace, io, moment_tensor
import numpy as np

engine = gf.LocalEngine(store_superdirs=['/home/zmaw/u300596/Documents/bachelor/gf'])
#distances = [100., 200., 300.]
moment = moment_tensor.magnitude_to_moment(-1.0)

source = gf.MTSource(
    time = util.str_to_time('2016-11-15 13:00:00.000'),
    depth= 80.0,
    mdd=-1*moment,
    mnd=-0.22*moment,
    med=0.69*moment,
    stf=gf.BoxcarSTF(duration=0.3))
#Anlegen einer DoubleCouple-Quelle

stations = (('FLBP1',-204.737835,-214.134002),('FLBP2',311.087988,10.094441),('FLBP3',-233.243096,364.604067), ('3171',996.004499,676.632386), \
            ('3172',489.160972,641.192668),('3G13',741.441641,341.971320))#,('CENT',-56.841044,174.671605),('GEO1',-71.522444,172.246905), \
 # ('GEO2',-47.986544,184.375505),('GEO3',-51.319044,163.431305))
components = ('HHZ'),('HHE'),('HHN')
# dies hier ist eine list comprehension gemaess a = [x+1 for x in range(0,4,1)] nur dass x+1 hier durch gf.Target() ersetzt ist
targets = [gf.Target(
  quantity='displacement',
  codes=('','%s' %stats[0],'','%s' %comps),	#runde Klammern, deshalb Tupel
  north_shift=stats[1],
  east_shift=stats[2],
  depth=0.,
  store_id='test_highsr_gf')
  for comps in components
  for stats in stations]

#Angelegt werden 10 Empfänger, die jeweils um 100 km getrennt in nörliche Richtung von der Quelle mit den angegeben Parametern entfernt sind

response = engine.process(source, targets)
#Erzeugung der synthetischen Spuren

traces = []
for resp_source, resp_target, resp_trace in response.iter_results():	#hier wird ueber die in response enthaltenen Eintraege iteriert und diese den Variablen resp_source, resp_target und resp_trace zugewiesen, lediglich resp_trace wird dann weiter verwand
  traces.append(resp_trace)
  #erzeugte Spuren werden in der Liste traces zusammengefasst - anscheinend stehen in response die Variablen iter_results
  
trace.snuffle(traces)
#Anzeigen der Spuren mit Snuffler
datname = 'synth_seism_HHZ.mseed'

io.save(traces, datname)
#Abspeichern der Liste traces in einer mseed-Datei
