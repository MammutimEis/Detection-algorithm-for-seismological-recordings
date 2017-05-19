# Detection-algorithm-for-seismological-recordings

## File description
#### triggertoolbox.py
The trigger toolbox Provides functions for detection algorithms. The module is imported by the other programs
#### plotCharacteristcFunction.py 	
Plots a seismogram and the corresponding characteristic function of the recursive STA/LTA.
#### plotTrigger.py
Applies the recursive STA/LTA trigger to the data and plots the recursive STA/LTA trigger on top of the seismograms.
#### syntheticEventStatisic.py
A trigger statistic with synthetic events. Three methods are implemented which apply bandpass filter with different corner frequencies to the data. An event statistic evaluates the performance of the trigger. 
#### totalGroundMovement.py 
The only file which do not import the trigger toolbox. Calculates the total ground movement, applies a recursive STA/LTA trigger to it and provides a trigger statistic as well. The code needs more comments and a better structure
#### synthetics_network.py
The file is only used for the generation of synthetic events
