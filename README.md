`Warning: This project is still under development and is not fully mature`

This project has been packaged.
Draft documentation website on [gh-pages](https://nve.github.io/FlomKart/)

# Installation

Install the package with the following code:

``` r
library(devtools)
install_github("NVE/FlomKart")
```

# FlomKart
Flood frequency analysis for the FlomKart project. The output data of the repo is fed into the FlomKart_ShinyApp visualization tool.
This is still work in progress!

# Data
Flood records and station metadata (catchment descriptors) can be found in the `/data` folder.

# Methodology
Only station with more than 30 years of data are used.
We calculate fitted parameters for 5 distirbutions:
- Generalized Logistics (GL)
- Generalized extrem value (GEV)
- Pearson III
- Gumbel
- Gamma

and 4 parameter estimation methods:
- Ordinary moments (mom)
- Linear moments (Lmom)
- Maximum likelihood estimation (mle)
- Bayesian inference (bayes)

# Filing issues

Please try to follow those guidelines for filing issues:

One issue for one purpose. Don't add more than one bug, feature request, documentation request, question, etc.. on to the same issue.

- If you've found a bug, thanks for reporting!
- If you've a request of some kind, e.g., feature request or documentation request, it'd be much appreciated if you could add **[Request]** at the beginning of the title. This helps us to prioritise easily without having to go through the entire issue.
- If you need support, e.g., installation issues or upgrade issues, please add **[Support]** at the beginning of the title. This helps us to easily identify the most common support issues, and provide solutions in a separate page.
- If you have a general question, add **[Question]** at the beginning of the title.
- If you've an issue that doesn't fall into any of these categories, then it'd be helpful if you could add **[Misc]** to your title.
