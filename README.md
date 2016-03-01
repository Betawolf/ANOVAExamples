# ANOVA Examples

This `R` code generates the LaTeX markup for a filled-out example ANOVA.

Usage is simple.

```{R}
testdata <- c(16, 11, 34, 31, 22, 30, 24, 17, 36, 31, 24, 24)

# data3d takes the rows*cols*examples, and optionally an existing dataset
# (as here) or else a function to use for generation (default is rnorm). 
testmatrix <- data3d(2, 3, 2, testdata)

workedanova(testmatrix)
# ... output ...
```
