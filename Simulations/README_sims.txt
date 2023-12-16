This folder contains functions which can be called to simulate community assembly and first-step evolution, producing output files such as those in the Raw_Results folder.

COMMUNITY ASSEMBLY FUNCTIONS

These functions generate assembled communities with mutations and output parameters describing the resulting communities into a "Results" folder.

single_step_mut: Generate an assembled community at equilibrium, introduce a single mtuation, and allow the community to re-equilibrate. Allows for knock-out mutants, knock-in mutants, and knock-out mutants without extinctions (i.e., allowing for "rescued species").

single_step_mut_alt: Alternative community assembly procedures used in Figs S2 and S4. In particular, this function is used to generate communities with Dirichlet-distributed resource usage (as opposed to binary), or communities subject to the consumer resource model in Tikhonov & Monasson 2017.

single_step_mut_Rp: Alternative community assembly procedure used in Fig S6 (left side). The parent species can use a different number of resources on average compared to the rest of the community. Because this means that the parent species is no longer chosen with a probability proportional to its abundance, this also outputs "weights" indicating the relative likelihood of a given run.

HELPER FUNCTIONS

These are functions internally called by the functions above.

findEquilibrium, findEquilibriumTikhonov: Functions for calculating the community equilibrium according to an optimization procedure, based on the consumer resource model discussed in the main text or in Tikhonov & Monasson 2017.

generateOrganisms, generateOrganismsDirichlet, generateOrganismsTikhonov: Functions for randomly sampling a set of organisms from a common statistical distribution. We use (1) a sampling procedure with binary resource usage, (2) Dirichlet-distributed resource usage, and (3) binary resource usage with energy costs ("pure fitnesses") drawn according to the procedure in Tikhonov & Monasson 2017.
