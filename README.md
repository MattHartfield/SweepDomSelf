README FOR DOM SELF SWEEP

SLiM scripts used to simulate neutral diversity around a sweeping beneficial allele, while considering arbitrary levels of self-fertilisation and dominance. As used in the study "Selective sweeps under dominance and self-fertilisation". SLiM is available from https://messerlab.org/slim/.

There are two scripts. "StdVarSweep_Hard.slim" is for simulating a 'hard' sweep, where a derived mutations has a selective advantage from it's appearance. "StdVarSweepFromNeut.slim" simulates sweeps arising from standing variation. These programs are meant to be run from the command line. For each program one first has to run 'mkdir Mutations/ BIPops/ Seeds/' to generate directories to store outputs. In each script the haplotype size is set to 100,000 basepairs, and the homozygote selection coefficient of the beneficial allele equals 0.05.

To run "StdVarSweep_Hard.slim", use:

slim -seed <SEED> -d N=<POP SIZE> -d Theta=<MUT RATE> -d R=<REC RATE> -d sps=<OUTPUT HAPLOTYPES> -d numsamp=<NUMBER OF SAMPLES> -d h=<DOMINANCE> -d sfrate=<SELFING FREQUENCY> -d simID=<ID> StdVarSweep_Hard.slim

Where:
- seed is the starting seed for the random number generator
- N is the population size
- Theta = 4Nu is the population-level mutation rate
- R = 2Nr is the population-level recombination rate across the entire haplotype
- sps is the number of haplotypes to sample from the final population
- numsamp is the number of times 'sps' haplotypes are sampled per simulation run
- h is the dominance value
- sfrate is the frequency of self-fertilisation
- simID is the simulation ID, used to label outputs

To run "StdVarSweepFromNeut.slim", use:

slim -seed <SEED> -d N=<POP SIZE> -d Theta=<MUT RATE> -d R=<REC RATE> -d x0=<STARTING FREQUNCY> -d sps=<OUTPUT HAPLOTYPES> -d numsamp=<NUMBER OF SAMPLES> -d h=<DOMINANCE> -d sfrate=<SELFING FREQUENCY> -d simID=<ID> ./StdVarSweepFromNeut_V6.slim

Parameters are the same as in "StdVarSweep_Hard.slim", except 'x0' is the frequency at which the derived allele becomes selected for.

The 'Seeds' folder contains a list of random number seeds used in each simulation. 'BIPops' outputs the state of the population after the burn-in. 'Mutations' outputs diversity (in MS format) for 'sps' haplotype samples after the sweep has fixed.

