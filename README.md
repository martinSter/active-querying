# Active Querying Approach to Epidemic Source Detection on Contact Networks

Our implementation is heavily **based** on fast C code for temporal networks originally written by Petter Holme (thanks!). His code can be found here: https://github.com/pholme/tsir.

Our code can be compiled and run as follows:
1. Create a directory *o* for the object files.
2. Type `make` to compile with gcc.
3. Then run the code with `./run data/malawi.csv 0.3 0.01 24 0 13865926330996091410`.

In step 3, we first provide the network file, then the infection and recovery probability, and the duration of the epidemic. Finally, with the 0 we indicate that it is an undirected network and the last parameter is a random seed for the random number generator. Note that the seed can be generated by running the Python script rnd.py.

In the header file run.h you can change the following parameters:

* Minimal outbreak size that is required for the experiments (`MINOUTBREAK`).
* The factor how the cutoff for the prior on the duration is computed (`NUMT`). We define a prior on the interval `[0, NUMT * T]`, where `T` is the true duration.
* Number of Monte-Carlo simulations (`NSIM`).
* Number of experiments (`NEXP`).
* Number of queries (`NQUERIES`).

The example network provided in the directory *data* is a network of contacts between residents of a village in Malawi. The data has been introduced by [Ozella et al. (2021)](https://epjdatascience.springeropen.com/articles/10.1140/epjds/s13688-021-00302-w) and can be found [here](http://www.sociopatterns.org/datasets/contact-patterns-in-a-village-in-rural-malawi/). Note that we modfied the original data to get a less granular network.