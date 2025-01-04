# EPBayes

## Simulation
Simulation scripts for 1D NPMLE-based empirical partially Bayes (1D-NPMLE), 2D NPMLE-based empirical partially Bayes (2D-NPMLE), Welch Approximation (Welch), Pooled t-test (t-test), Beherens-Fisher test (B-F), Equal-variance NPMLE-based empirical partially Bayes (EV-NPMLE).

## Real Dataset
1. Applying all algorithms (1D-NPMLE, 2D-NPMLE, Welch, t-test, B-F, EV-NPMLE) on Hedenfalk
2. Incorporating weight information estimated by voom and voomBygroup with all algorithms on PBMC1

## Data and Graphs
1. Distributions (G(\sigma^2_1/\sigma^2_2), H(\sigma^2_1, \sigma^2_2)) estimated based on real datasets (Hedenfalk, PBMC1) and visualization
2. Comparison between all algorithms concerning Power and FDR for (n1 $in$ (3,4,5,6,7), n2 = (n1: n1+6)) in equal-variance scenarios and unqeual-vairance scenarios.
    * Visualizations of distributions are generated via visualizer.Rmd and unequal_visualizer.Rmd
    * Separated Data contains simulation results for each of the algorithms
    * Full Data contains the combined simulation results processed via Data_organizer.R
