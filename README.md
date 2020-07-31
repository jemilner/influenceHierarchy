# influenceHierarchy

The code to fit the model from ... This model captures interactions between social animals through their movement behaviours.

## Fitting the Model

Almost all of the below information regards the gfit-setup-******.R files. These files load in data, set up the model, parameters etc for their respective analyses:
 - gfit-setup-sim.R - simulation analysis in Section 4.1
 - gfit-setup-baboon.R - baboon analyses in Sections 4.2 and 6
 - gfit-setup-sim-mass.R and gfit-setup-sim-mass-parallel.R - the reliability testing in Section 5

gfit-iter.R contains the MCMC loop. This ony needs to be edited when changing the number of BM states in the model (see below).

All of the other files contain functions used throughout the MCMC and contain comments explaining their function.

### Data Format

Data is loaded into the model through an RData file containing four columns:

1. the x coordinate
2. the y coordinate
3. time
4. animal ID

Currently, the code is setup for the case where the animal IDs are 1:total.animals.

### Number of Animals

The number of animals in the model is dictated by the number of distinct animal IDs in the data and is stored in:

```
total.animals <- length(animal.id)
```

The only section of code the user needs to be aware of with regards to the number of animals in the data/model is in gfit-iter.R:

```
cat(file = file.states1, all.data[, col.state, 1], "\n", append = TRUE)
cat(file = file.states2, all.data[, col.state, 2], "\n", append = TRUE)
cat(file = file.states3, all.data[, col.state, 3], "\n", append = TRUE)
cat(file = file.states4, all.data[, col.state, 4], "\n", append = TRUE)
cat(file = file.states5, all.data[, col.state, 5], "\n", append = TRUE)
```

These lines output the state estimations for each animal (5 animals in this example) and must correspond to the number of animals in the data. The filename variables are set in the relevant gfit-setup-******.R.

### Number of States

The number of Brownian motion (BM) states in the model can easily be amended in the code. The current code does not allow for Ornstein-Uhlenbeck (OU) states for leading/independent animals (see Section 2.3).
```
#dictates how many BM states are included
num.bm <- 2
```

A number of initial and tuning parameters for rho are required to correspond to the number of BM states. In this example there are 2 BM states to represent 'slow' and 'fast' movement.
```
#rho initialisation
rho <- c(1, 1)

#rho tuning parameters
prop.slow <- 0.02
prop.fast <- 0.10
```

The sampling of the new rho parameters in gfit-iter.R then needs to be updated to represent the new number of BM states:
```
#a vector for the new rho samples
new.rho <- c(NA, NA)
new.rho[1] <- rho[1] + rnorm(1, 0, prop.slow)
new.rho[2] <- rho[2] + rnorm(1, 0, prop.fast)
```

Note, as metioned in Section 3.2, to keep state labelling consistant, rho_1 < ... < rho_n for n BM states. The following line (which checks if the new movement parameter samples are in line with our boundaries - Section 3.2) in gfit-iter.R will also need updating with appropriate conditions for the number of BM states:
```
if(new.alpha > 0 & new.sigma > 0 & all(new.rho > 0) & new.sigma < new.rho[2] & new.rho[2] > new.rho[1])
```

In particular, new.sigma < new.rho[n] and new.rho[1] < ... < new.rho[n].

### Tuning Parameters

The model contains the following tuning parameters:
```
#for the parameters of the OU states
prop.alpha <- 0.02
prop.sigma <- 0.02
#for the parameters of the BM states
prop.slow <- 0.02
prop.fast <- 0.10

#partial.var is used to resample locations that are missing in the first observation
#this basically just creates a random walk around the previous sample
partial.var <- diag(1) * 0.1

#lambda_max (see Section 3.1)
max.rate <- 0.2

#number of MCMC iterations
iter <- 10000

#the lower and upper bounds of the trajectory length that is updated
traj.min <- 3
traj.max <- 12
#how many trajectory updates are performed each MCMC iteration
traj.iter <- ceiling(num.obs / ((traj.min + traj.max) / 2)) * total.animals
```

### Initial Parameter Values

Initial values are required for the movement parameters:
```
alpha <- 0.7
sigma <- 1
rho <- c(1, 1)
```

Initial values are also required for the behaviour parameters - the behaviour states and the transition rates. In all our analyses in the paper, the transition rates are initialised uniformly:
```
init.lambda <- max.rate / num.states
```
whilst the behaviour states are mostly random. Each animal is assigned a random state (either BM or subordinate) at each observation. However, if that state is subordinate and it caused a cyclic hierarchy, they are randomly assigned a BM state.

In the case of partial observations (Section 3.1), we need to initialise the locations of the missing animals. To do so, in our analyses, we simply linearly interpolated the missing location from the previous and next recorded observations.

Finally, we need to initialise the switching data that facilitates the initial behaviour states. In short, the initial switching times are randomly sampled in the intervals that contain a switch as dictated by the behaviour states assigned at the observations.

### Priors

The prior Dirichlet distribution for the transitions rates is parameterised with shape1 and shape2. shape1 corresponds to switching states, shape2 corresponds to staying in the current state.
```
shape1 <- shape2 <- 1 / num.states
```

## Reproduction

In all of the analyses undertaken for this paper, partial observations, behaviour states and the corresponding switching times are initialised as described above.

### Simulation Results (Section 4.1)

Data: 
```
linear-group5-seq-two-speed-missing-data.RData
```

Seeds:
```
#main analysis
set.seed(1107)

#secondary analysis
set.seed(500)
```

Number of states:
```
num.bm <- 2
```

Tuning Parameters:
```
prop.alpha <- 0.02
prop.sigma <- 0.02
prop.slow <- 0.02
prop.fast <- 0.10

partial.var <- diag(1) * 0.1

max.rate <- 0.2

iter <- 1400000

traj.min <- 3
traj.max <- 12
traj.iter <- ceiling(num.obs / ((traj.min + traj.max) / 2)) * total.animals
```

Initial Parameter Values:
```
#main analysis (mcmc.run <- 1)
alpha <- 0.7
sigma <- 1
rho <- c(1, 1)

#secondary analysis (mcmc.run <- 2)
alpha <- 2
sigma <- 0.5
rho <- c(2, 4)

#initial transition rates
init.lambda <- max.rate / num.states
```

Priors:

We used an uninformative Dirichlet prior
```
shape1 <- shape2 <- 1 / num.states
```

### Baboon Results (Section 4.2)

Data: see Data Availability section in the paper

Seeds:
```
#main analysis
set.seed(1107)

#secondary analysis
set.seed(500)
```

Number of states:
```
num.bm <- 2
```

Tuning Parameters:
```
prop.alpha <- 0.0005
prop.sigma <- 0.005
prop.slow <- 0.001
prop.fast <- 0.01

partial.var <- diag(1) * 2

max.rate <- 0.2

iter <- 500000

traj.min <- 3
traj.max <- 40
traj.iter <- ceiling(num.obs / ((traj.min + traj.max) / 2)) * total.animals
```

Initial Parameter Values:
```
#main analysis (mcmc.run <- 1)
alpha <- 0.1
sigma <- 1
rho <- c(1, 1)

#secondary analysis (mcmc.run <- 2)
alpha <- 0.05
sigma <- 0.5
rho <- c(0.5, 0.5)

#initial transition rates
init.lambda <- max.rate / num.states
```

Priors:

We used an uninformative Dirichlet prior
```
shape1 <- shape2 <- 1 / num.states
```

### Reliability (Section 5)

There are three lines of code to consider when running the code relating to Section 5. All are in gfit-setup-mass-parallel.R:

```
#define the number of cores to be used
num.cores <- 2

#which parameter set you are running
param.set <<- 1

#replace 1:4 with the 'seeds' you are running
parLapply(clust, 1:4, main)
```

Due to the large number of analyses, the reliability testing (or "mass testing" as it called throughout the code/data), is setup to run in parallel over a number of cores. The data for this is named as such:

linear-set-[i]-seed-[j]-data.RData

where [i] is the parameter set and [j] is the seed used to create the data (ranging from 1-100 for each parameter set). gfit-setup-sim-mass.R then controls the data load, tuning parameters etc depending on which parameter set and seeds have been set.

### Data Thinning (Section 6)

The data.freq variable controls whether to fit the model to the full baboon dataset or a thinned dataset (with 1 = full dataset, 5 = thin by 5, 20 = thin by 20)
```
#options: 1, 5, or 20
data.freq <- 1
```

Data: see Data Availability section in the paper

Seeds:
```
set.seed(1107)
```

Number of states:
```
num.bm <- 2
```

Tuning Parameters:
```
partial.var <- diag(1) * 2

#thin-by-5 (data.freq <- 5)
prop.alpha <- 0.001
prop.sigma <- 0.025
prop.slow <- 0.005
prop.fast <- 0.050

max.rate <- 0.2

iter <- 500000

traj.min <- 3
traj.max <- 12
traj.iter <- ceiling(num.obs / ((traj.min + traj.max) / 2)) * total.animals

#thin-by-20 (data.freq <- 20)
prop.alpha <- 0.002
prop.sigma <- 0.05 
prop.slow <- 0.01
prop.fast <- 0.1

max.rate <- 0.05

iter <- 500000

traj.min <- 3
traj.max <- 6
traj.iter <- ceiling(num.obs / ((traj.min + traj.max) / 2)) * total.animals
```

Initial Parameter Values:
```
alpha <- 0.1
sigma <- 1
rho <- c(1, 1)

#initial transition rates
init.lambda <- max.rate / num.states
```

Priors:

We used an uninformative Dirichlet prior
```
shape1 <- shape2 <- 1 / num.states
```

## Output

The following files are outputted, where [n] = 1,...,total.animals and [analysis-name] is set by:
```
dat <- "sim1"
```

 - group-[n]states-[analysis-name].txt - these contain the current state estimates for each animal at each observation and switching time
 - group-acc-[analysis-name].txt - contains the acceptance rates of the movement parameters and trajectory updates since the last output
 - group-alpha-[analysis-name].txt - contains the alpha parameter samples
 - group-obs[analysis-name].txt - contains the indexes of the augmented data (observations and switching times) that relate to observed data
 - group-rates[analysis-name].txt - contains the transition rate samples
 - group-rho-[analysis-name].txt - contains the rho parameter samples
 - group-sigma-[analysis-name].txt - contains the sigma parameter samples
 - group-times-[analysis-name].txt - contains the times of the observations and switching times to correspond with the state outputs in group-[n]states-[analysis-name].txt

## General Guidance

As alluded to in Section 4.2, the number of BM states included in the model is an important factor in allowing the inference to mix well. In particular, too few BM states may cause the subordinate states to act as a pseudo-BM state to capture different speeds of movement found in the data. Therefore, it may be pertinent to establish how many BM states is appropriate for the analysis at hand before fitting the model to the data.
