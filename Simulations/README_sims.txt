This folder contains functions which can be called to simulate community assembly and first-step evolution, producing output files such as those in the Raw_Results folder.

COMMUNITY ASSEMBLY FUNCTIONS

These functions generate assembled communities with mutations and output parameters describing the resulting communities into a "Results" folder.

comm_assembly: Generate an assembled community at equilibrium without mutation. This function allows either for a community to be generated either with explicit community assembly parameters (S and epsilon), or with implicit community assembly parameters by setting a target S* (astar = S*/R, phi = S*/S).

single_step_mut: Generate an assembled community at equilibrium, introduce a single mtuation, and allow the community to re-equilibrate. Allows for knock-out mutants, knock-in mutants, and knock-out mutants without extinctions (i.e., allowing for "rescued species").

single_step_mut_alt: Alternative community assembly procedures used in the top half of Figs S3 and S8. In particular, this function is used to generate communities with Dirichlet-distributed resource usage (as opposed to binary), or communities subject to the consumer resource model in Tikhonov & Monasson 2017.

single_step_mut_Rp: Alternative community assembly procedure used in Fig S9A. The parent species can use a different number of resources on average compared to the rest of the community. Because this means that the parent species is no longer chosen with a probability proportional to its abundance, this also outputs "weights" indicating the relative likelihood of a given run.

single_step_mut_gaussian_uptakes: Alternative community assembly procedure used in Fig S3C. This generates communities with Dirichlet-distributed resource usage, but with a large shape parameter so that resource usage is approximately Gaussian.

single_step_mut_no_tradeoffs: Alternative community assembly procedure used in Fig S8C-D. This generates communities where the sampling process has no metabolic tradeoffs, i.e. where using additional resources does not correlate with using each individual resource less efficiently.

single_step_mut_wide_R0_variation: Alternative community assembly procedure used in Fig S9B-D. This generates communities where the number of resources used by each organism is uniformly distributed.

single_step_mut_wide_resource_variation: Alternative community assembly procedure used in Fig S7B-C. This generates communities where the external resource supply is exponentially distributed, a wider distribution than used elsewhere in the text.

HELPER FUNCTIONS

These are functions internally called by the functions above.

findEquilibrium, findEquilibriumTikhonov: Functions for calculating the community equilibrium according to an optimization procedure, based on the consumer resource model discussed in the main text or in Tikhonov & Monasson 2017.

generateOrganisms, generateOrganismsDirichlet, generateOrganismsNoTradeoffs, generateOrganismsTikhonov, generateOrganismsWideR0: Functions for randomly sampling a set of organisms from a common statistical distribution. We use (1) a sampling procedure with binary resource usage, (2) Dirichlet-distributed resource usage, (3) binary resource usage without metabolic tradeoffs, (4) binary resource usage with energy costs ("pure fitnesses") drawn according to the procedure in Tikhonov & Monasson 2017, and (5) binary resource usage where the number of resources used by each species follows a uniform distribution.
