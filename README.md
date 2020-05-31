# sc_simulate
Use Splatter to simulate single cell counts matrices and Minnow to simulate dscRNA-seq reads.

Main file: splatter_simulate_channel_patients.R

File of edited functions: multi_batch_functions.R

Edit parameters at top of splatter_simulate_channel_patients.R file and in splatSimulate_multi_batches() function. 

Currently can load 3 batch effects: pool (Batch effect in Splat algorithm), patients and channel and one group effect (larger effect size and differential expression between groups).

I've used the in-built batch effect to model pools

Batch numbers:
Can't set cell number directly must be done using batches, e.g. for 3 pools of 100 cells each `batchCells = c(100,100,100)`.

I've used .facLoc and .facScale parameters for batch effect, usage found here:
- https://github.com/Oshlack/splatter/blob/master/R/splat-simulate.R#L345
- https://github.com/Oshlack/splatter/blob/master/R/splat-simulate.R#L823

Quite complex relationship. Mean and sd of log normal distribution. I've just increased both for batch effects I think are bigger, not sure if this is right.

I've written nPatient and nChannels factors in. Somewhat *not* intuitively, its coded so nPatient is number of patients per pool- i.e. 14 and nChannels total- i.e. 70.

I have multiplied batch effects together for pool, channel and patient. 
As a default I have made patient > pool > channel for .facLoc and .facScale parameters

Group properties:
I have coded groups (healthy, mild, severe, sepsis,...) as groups. Groups have a larger effect on the variance of the count matrix than batch effects. As they have DE genes.

`group.prob = c(0.2,0.2,0.2,0.2,0.2)`, implies each group is equally likely, but can change it if certain groups are more likely than others.

`de.prob = c(0.1, 0.1, 0.1, 0.2, 0.2)`, is likelihood for DE for each group. 0.1 is the programme's default for each group. I have changed the last two values to 0.2








