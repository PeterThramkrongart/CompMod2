Computational Modeling - Assignment 2
================
Riccardo Fusaroli
29/01/2019

## In this assignment we learn how to assess rates from a binomial distribution, using the case of assessing your teachers’ knowledge of CogSci

N.B. there is a second part at the bottom for next week.

### First part

You want to assess your teachers’ knowledge of cognitive science. “These
guys are a bunch of drama(turgist) queens, mindless philosophers,
chattering communication people and Russian spies. Do they really know
CogSci?”, you think.

To keep things simple (your teachers should not be faced with too
complicated things): - You created a pool of equally challenging
questions on CogSci - Each question can be answered correctly or not (we
don’t allow partially correct answers, to make our life simpler). -
Knowledge of CogSci can be measured on a scale from 0 (negative
knowledge, all answers wrong) through 0.5 (random chance) to 1 (awesome
CogSci superpowers)

This is the data: - Riccardo: 3 correct answers out of 6 questions -
Kristian: 2 correct answers out of 2 questions (then he gets bored) -
Josh: 160 correct answers out of 198 questions (Josh never gets bored) -
Mikkel: 66 correct answers out of 132 questions

Questions:

1.  What’s Riccardo’s estimated knowledge of CogSci? What is the
    probability he knows more than chance (0.5) \[try figuring this out.
    if you can’t peek into chapters 3.1 and 3.2 and/or the slides\]?

<!-- end list -->

  - First implement a grid approximation (hint check paragraph 2.4.1\!)
    with a uniform prior, calculate the posterior and plot the results
  - Then implement a quadratic approximation (hint check paragraph
    2.4.2\!).
  - N.B. for the rest of the exercise just keep using the grid
    approximation (we’ll move to quadratic approximations in two
    classes)

## Question 1 - Ricardo’s knowledge

``` r
# loading packages
library("rstan")
```

    ## Warning: package 'rstan' was built under R version 3.6.2

    ## Loading required package: StanHeaders

    ## Warning: package 'StanHeaders' was built under R version 3.6.2

    ## Loading required package: ggplot2

    ## rstan (Version 2.19.3, GitRev: 2e1f913d3ca3)

    ## For execution on a local, multicore CPU with excess RAM we recommend calling
    ## options(mc.cores = parallel::detectCores()).
    ## To avoid recompilation of unchanged Stan programs, we recommend calling
    ## rstan_options(auto_write = TRUE)

    ## For improved execution time, we recommend calling
    ## Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
    ## although this causes Stan to throw an error on a few processors.

``` r
pacman::p_load(rethinking, ggplot2, dplyr, tidyverse, tidyr, gridExtra, tidybayes)

### Riccardo's estimated knowledge

##grid aproximation##

# define grid
p_grid <- seq(from=0 , to=1 , length.out=1000)

## Priors
uni_prior <- dunif(p_grid, 0,1)
stan_prior <- dnorm(p_grid, mean = 0.8, sd = 0.2)

# compute likelihood at each value in grid
ric_likelihood <- dbinom(3 , size=6 , prob=p_grid )

# compute product of likelihood and prior
uns_ric_posterior <- ric_likelihood * uni_prior
ric_posterior <- uns_ric_posterior/sum(uns_ric_posterior)


# drawing samples
ric_samples <- sample( p_grid , prob=ric_posterior , size=1e4 , replace=TRUE )

#posterior
sum( ric_posterior[ p_grid > 0.5 ] ) # 50 % probability that Ricardo's knowledge above 0.5
```

    ## [1] 0.5

``` r
p_grid[which.max(ric_posterior)]
```

    ## [1] 0.4994995

``` r
## quadratic approximation
ric.qa <- rethinking::map(
  alist(
    w ~ dbinom(6,p) , # binomial likelihood
    p ~ dunif(0,1)), # uniform prior
    data=list(w=3) )

# display summary of quadratic approximation
precis( ric.qa )
```

    ##   Mean StdDev 5.5% 94.5%
    ## p  0.5    0.2 0.17  0.83

``` r
#curve( dnorm( x , 0.5 , 0.2 ) , lty=2 , add=TRUE )
```

## Question 2

2.  Estimate all the teachers’ knowledge of CogSci. Who’s best? Use grid
    approximation. Comment on the posteriors of Riccardo and Mikkel. 2a.
    Produce plots of the prior, and posterior for each teacher.

<!-- end list -->

``` r
## posterior
calc_teacher_knowledge <- function(n_correct, n_question, prior, length_out = 10000){
  p_grid <- seq(0, 1, length.out= length(prior))
  
  likelihood <- dbinom(n_correct, size = n_question, prob=p_grid)
  
  unstd.posterior <- prior*likelihood
  
  posterior <- unstd.posterior/sum(unstd.posterior)
  
  sampl <- sample(p_grid, prob = posterior, size=1e4 , replace=TRUE)
  
  max_posterior <- p_grid[which.max(posterior)]
  
  HPDI <- HPDI(sampl,prob=0.8) # if bell-shaped, prob doesn't matter
  SD <- sd(sampl)
  
  return(list(teacher_posterior = posterior, 
              likelihood = likelihood,
              max_posterior = max_posterior, 
              teacher_HPDI = HPDI,
              sample = sampl,
              SD=SD))
}

## plotting

pretty_plot <- function(p_grid, prior, posterior, title = " "){
  # define data
  d <- tibble(p_grid = p_grid, 
              prior = prior/length(p_grid),
              posterior = posterior)

  # make to long format
  d <- d %>% 
    pivot_longer(cols = c("prior", "posterior"), names_to = "name", values_to = "value")
  
  # make a 
  p <- ggplot(d, aes(x = p_grid, y = value, color = name)) + 
    geom_line() + 
    labs(x = "CogSci Knowledge (p)", y = "Density", title = title) + 
    theme_bw() + 
    ggplot2::theme(panel.background = element_rect(fill = "White"),
                   panel.border = element_blank()) +
    scale_colour_brewer(palette = "Dark2", direction = 1)
  return(p)
}

#function for giving summary stats
summary_stats <- function(teacher){
  print(teacher$teacher_HPDI)
  print(teacher$max_posterior)
  print(teacher$SD)
}

## Riccardo
ric <- calc_teacher_knowledge(3, 6, uni_prior)
ric_plot <- pretty_plot(p_grid, uni_prior, ric$teacher_posterior, title = "Riccardo's estimated knowledge")
summary_stats(ric)
```

    ##      |0.8      0.8| 
    ## 0.2702703 0.7117117 
    ## [1] 0.4994995
    ## [1] 0.1656042

``` r
## Kristian
kris <- calc_teacher_knowledge (2, 2, uni_prior)
kris_plot <- pretty_plot(p_grid, uni_prior, kris$teacher_posterior, title = "Kristian's estimated knowledge" )
summary_stats(kris)
```

    ##      |0.8      0.8| 
    ## 0.5945946 1.0000000 
    ## [1] 1
    ## [1] 0.191587

``` r
## Josh
josh<- calc_teacher_knowledge(160, 198, uni_prior)
josh_plot <- pretty_plot(p_grid, uni_prior, josh$teacher_posterior, title = "Josh's estimated knowledge" )
summary_stats(josh)
```

    ##      |0.8      0.8| 
    ## 0.7727728 0.8428428 
    ## [1] 0.8078078
    ## [1] 0.02792883

``` r
## Mikkel
mik <- calc_teacher_knowledge(66,132, uni_prior)
mik_plot <- pretty_plot(p_grid, uni_prior, mik$teacher_posterior, title = "Mikkel's estimated knowledge" )
summary_stats(mik)
```

    ##      |0.8      0.8| 
    ## 0.4424424 0.5535536 
    ## [1] 0.4994995
    ## [1] 0.0434567

``` r
## plotting knowledge of everyone at one
grid.arrange(ric_plot, kris_plot, josh_plot, mik_plot, nrow=2)
```

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
sum(kris$sample<josh$sample)/1e4*100 #josh is smarter than Kristian
```

    ## [1] 51.86

## Question 3

3.  Change the prior. Given your teachers have all CogSci jobs, you
    should start with a higher appreciation of their knowledge: the
    prior is a normal distribution with a mean of 0.8 and a standard
    deviation of 0.2. Do the results change (and if so how)? 3a. Produce
    plots of the prior and posterior for each teacher.

<!-- end list -->

``` r
## Riccardo
ric3 <- calc_teacher_knowledge (3, 6, stan_prior)
ric3_plot <- pretty_plot(p_grid, stan_prior, ric3$teacher_posterior, title = "Riccardo's estimated knowledge" )
summary_stats(ric3)
```

    ##      |0.8      0.8| 
    ## 0.4744745 0.7987988 
    ## [1] 0.6466466
    ## [1] 0.1260536

``` r
## Kristian
kris3 <- calc_teacher_knowledge (2, 2, stan_prior)
kris3_plot <- pretty_plot(p_grid, stan_prior, kris3$teacher_posterior, title = "Kristian's estimated knowledge" )
summary_stats(kris3)
```

    ##      |0.8      0.8| 
    ## 0.6946947 1.0000000 
    ## [1] 0.8898899
    ## [1] 0.1323245

``` r
## Josh
josh3 <- calc_teacher_knowledge(160, 198, stan_prior)
josh3_plot <- pretty_plot(p_grid, stan_prior, josh3$teacher_posterior, title = "Josh's estimated knowledge" )
summary_stats(josh3)
```

    ##      |0.8      0.8| 
    ## 0.7697698 0.8398398 
    ## [1] 0.8078078
    ## [1] 0.02772858

``` r
## Mikkel
mik3 <- calc_teacher_knowledge(66,132, stan_prior)
mik3_plot <- pretty_plot(p_grid, stan_prior, mik3$teacher_posterior, title = "Mikkel's estimated knowledge" )
summary_stats(mik3)
```

    ##      |0.8      0.8| 
    ## 0.4604605 0.5685686 
    ## [1] 0.5135135
    ## [1] 0.04240518

``` r
## plotting
grid.arrange(ric3_plot, kris3_plot,josh3_plot, mik3_plot, nrow=2)
```

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Question 4

4.  You go back to your teachers and collect more data (multiply the
    previous numbers by 100). Calculate their knowledge with both a
    uniform prior and a normal prior with a mean of 0.8 and a standard
    deviation of 0.2. Do you still see a difference between the results?
    Why?

<!-- end list -->

``` r
## Ricardo
ric4_uni <- calc_teacher_knowledge (300, 600, uni_prior)
ric4_stan <- calc_teacher_knowledge(300, 600, stan_prior)


## Kristian
kris4_uni <- calc_teacher_knowledge (200, 200, uni_prior)
kris4_stan <- calc_teacher_knowledge(200, 200, stan_prior)


## Josh
josh4_uni <- calc_teacher_knowledge (1600, 1980, uni_prior)
josh4_stan <- calc_teacher_knowledge(1600, 1980, stan_prior)


## Mikkel
mik4_uni <- calc_teacher_knowledge (666, 1320, uni_prior)
mik4_stan <- calc_teacher_knowledge(666, 1320, stan_prior)

#summary stats with uniform priors
summary_stats(ric4_uni)
```

    ##      |0.8      0.8| 
    ## 0.4744745 0.5255255 
    ## [1] 0.4994995
    ## [1] 0.02025934

``` r
summary_stats(kris4_uni)
```

    ##     |0.8     0.8| 
    ## 0.991992 1.000000 
    ## [1] 1
    ## [1] 0.005013777

``` r
summary_stats(josh4_uni)
```

    ##      |0.8      0.8| 
    ## 0.7967968 0.8188188 
    ## [1] 0.8078078
    ## [1] 0.008883014

``` r
summary_stats(mik4_uni)
```

    ##      |0.8      0.8| 
    ## 0.4844845 0.5195195 
    ## [1] 0.5045045
    ## [1] 0.01371008

``` r
#summary stats with an informed prior
summary_stats(ric4_stan)
```

    ##      |0.8      0.8| 
    ## 0.4774775 0.5285285 
    ## [1] 0.5035035
    ## [1] 0.02029656

``` r
summary_stats(kris4_stan)
```

    ##     |0.8     0.8| 
    ## 0.991992 1.000000 
    ## [1] 1
    ## [1] 0.00498643

``` r
summary_stats(josh4_stan)
```

    ##      |0.8      0.8| 
    ## 0.7967968 0.8188188 
    ## [1] 0.8078078
    ## [1] 0.008878313

``` r
summary_stats(mik4_stan)
```

    ##      |0.8      0.8| 
    ## 0.4894895 0.5245245 
    ## [1] 0.5055055
    ## [1] 0.01370162

## Question 5

5.  Imagine you’re a skeptic and think your teachers do not know
    anything about CogSci, given the content of their classes. How would
    you operationalize that belief?

### Second part: Focusing on predictions

Questions to be answered (but see guidance below): 1- Write a paragraph
discussing how assessment of prediction performance is different in
Bayesian vs. frequentist models 2- Provide at least one plot and one
written line discussing prediction errors for each of the teachers.

This is the old data: - Riccardo: 3 correct answers out of 6 questions -
Kristian: 2 correct answers out of 2 questions (then he gets bored) -
Josh: 160 correct answers out of 198 questions (Josh never gets bored) -
Mikkel: 66 correct answers out of 132 questions

This is the new data: - Riccardo: 9 correct answers out of 10 questions
(then he freaks out about teaching preparation and leaves) - Kristian: 8
correct answers out of 12 questions - Josh: 148 correct answers out of
172 questions (again, Josh never gets bored) - Mikkel: 34 correct
answers out of 65 questions

``` r
## calculating new knowledge
ric_new <- calc_teacher_knowledge(9, 10, ric3$teacher_posterior)
kris_new <- calc_teacher_knowledge(8, 12, kris3$teacher_posterior)
josh_new <- calc_teacher_knowledge(148, 172, josh3$teacher_posterior)
mik_new <- calc_teacher_knowledge(34, 65, mik3$teacher_posterior)


## function for a plot of posterior comparison
comp_distrib_plot <- function(old_posterior, new_posterior, title = ""){
  # define data
  d <- tibble(old_posterior = old_posterior,
              new_posterior = new_posterior,
              p_grid = p_grid)
  d <- d %>% 
      pivot_longer(cols = c("old_posterior", "new_posterior"), names_to = "name", values_to = "value")
  
  # make a 
  p <- ggplot(d, aes(x = p_grid, y = value, colour = name)) + 
    geom_line()+ 
    labs(y = "Density", x = " Posterior of teacher's knowledge", title = title)
  p
  
  return(p)
}

## running the plot comparison function on data of all 4 teachers
ric_comp <- comp_distrib_plot(ric3$teacher_posterior, ric_new$teacher_posterior, title = "Riccardo's estimated knowledge 2019 & 2020")
kris_comp <- comp_distrib_plot(kris3$teacher_posterior, kris_new$teacher_posterior, title = "Kristian's estimated knowledge 2019 & 2020")
josh_comp <- comp_distrib_plot(josh3$teacher_posterior, josh_new$teacher_posterior, title = "Josh's estimated knowledge 2019 & 2020")
mik_comp <- comp_distrib_plot(mik3$teacher_posterior, mik_new$teacher_posterior, title = "Mikkel's estimated knowledge 2019 & 2020")

all_comparison <- grid.arrange(ric_comp, kris_comp, josh_comp, mik_comp, nrow=2)
```

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
### histogram function

hist_function <- function(posterior, n_questions, bin_width, x = "", title = ""){
  
  #drawing samples
  sampling <- sample(size = 100000,
                     x = p_grid,
                     prob = posterior,
                     replace = T)
  pos_pred <- rbinom(100000,
                     size = n_questions,
                     prob = sampling)
  #doing some nice plots
  d <- as.data.frame(pos_pred)
 
   g <- ggplot(d, aes(pos_pred, fill = "pink"))+ 
     geom_histogram() + 
     labs(x = x, y = "Count", title = title) +
     stat_bin(binwidth = bin_width)
   return(g)
}

ric_his <- hist_function(ric3$teacher_posterior, 10, 1, x = "Correct answers (out of 10)", title = "Predictions of Riccardo")
kris_his <- hist_function(ric3$teacher_posterior, 12, 1, x = "Correct answers (out of 12)", title = "Predictions of Kristian")
josh_his <- hist_function(ric3$teacher_posterior, 172, 1, x = "Correct answers (out of 172)", title = "Predictions of Josh")
mik_his <- hist_function(ric3$teacher_posterior, 65, 1, x = "Correct answers (out of 65)", title = "Predictions of Mikkel")

grid.arrange(ric_his, kris_his, josh_his, mik_his, nrow=2)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# this function makes a histogram of errors
error_hist <- function(posterior, questions, correct, name){
  samples <-  sample(p_grid , prob=posterior , size=100000 , replace=TRUE )
  error_hist <- hist(samples * questions - correct,  main=name, breaks=questions, xlab = "Error Count")
  return(error_hist)
}

# this function makes a histogram of errors (used for large samples to control bin width)
error_hist2 <- function(posterior, questions, correct, name){
  samples <-  sample(p_grid , prob=posterior , size=100000 , replace=TRUE)
  error_hist <- hist(samples * questions - correct,  main=name, xlab = "Error Count")
  return(error_hist)
}
error_hist(ric3$teacher_posterior,10,9, "Riccardo Error")
```

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

    ## $breaks
    ##  [1] -9 -8 -7 -6 -5 -4 -3 -2 -1  0  1
    ## 
    ## $counts
    ##  [1]     1    51   618  3520 11562 23838 30439 21966  7439   566
    ## 
    ## $density
    ##  [1] 0.00001 0.00051 0.00618 0.03520 0.11562 0.23838 0.30439 0.21966 0.07439
    ## [10] 0.00566
    ## 
    ## $mids
    ##  [1] -8.5 -7.5 -6.5 -5.5 -4.5 -3.5 -2.5 -1.5 -0.5  0.5
    ## 
    ## $xname
    ## [1] "samples * questions - correct"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
error_hist(kris3$teacher_posterior,12,8, "Kristian Error")
```

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

    ## $breaks
    ##  [1] -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4
    ## 
    ## $counts
    ##  [1]     4    19   126   526  1738  4405  9257 15385 21371 24188 22981
    ## 
    ## $density
    ##  [1] 0.00004 0.00019 0.00126 0.00526 0.01738 0.04405 0.09257 0.15385 0.21371
    ## [10] 0.24188 0.22981
    ## 
    ## $mids
    ##  [1] -6.5 -5.5 -4.5 -3.5 -2.5 -1.5 -0.5  0.5  1.5  2.5  3.5
    ## 
    ## $xname
    ## [1] "samples * questions - correct"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
error_hist2(josh3$teacher_posterior,172,148, "Josh Error")
```

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

    ## $breaks
    ##  [1] -32 -30 -28 -26 -24 -22 -20 -18 -16 -14 -12 -10  -8  -6  -4  -2   0   2   4
    ## [20]   6   8
    ## 
    ## $counts
    ##  [1]     6    17    46   182   476  1216  2585  4621  8441 11561 15556 17143
    ## [13] 14445 11911  6596  3523  1310   309    50     6
    ## 
    ## $density
    ##  [1] 0.000030 0.000085 0.000230 0.000910 0.002380 0.006080 0.012925 0.023105
    ##  [9] 0.042205 0.057805 0.077780 0.085715 0.072225 0.059555 0.032980 0.017615
    ## [17] 0.006550 0.001545 0.000250 0.000030
    ## 
    ## $mids
    ##  [1] -31 -29 -27 -25 -23 -21 -19 -17 -15 -13 -11  -9  -7  -5  -3  -1   1   3   5
    ## [20]   7
    ## 
    ## $xname
    ## [1] "samples * questions - correct"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
error_hist2(mik$teacher_posterior,65,34, "Mikkel Error")
```

![](assignment-2-last-document_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->

    ## $breaks
    ##  [1] -14 -13 -12 -11 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4
    ## [20]   5   6   7   8   9  10  11  12
    ## 
    ## $counts
    ##  [1]     1     4    14    76   270   615  1384  3050  5094  8416 10984 12877
    ## [13] 14565 12865 10730  8539  5022  3027  1428   671   259    72    27     9
    ## [25]     0     1
    ## 
    ## $density
    ##  [1] 0.00001 0.00004 0.00014 0.00076 0.00270 0.00615 0.01384 0.03050 0.05094
    ## [10] 0.08416 0.10984 0.12877 0.14565 0.12865 0.10730 0.08539 0.05022 0.03027
    ## [19] 0.01428 0.00671 0.00259 0.00072 0.00027 0.00009 0.00000 0.00001
    ## 
    ## $mids
    ##  [1] -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5
    ## [13]  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5
    ## [25]  10.5  11.5
    ## 
    ## $xname
    ## [1] "samples * questions - correct"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"
