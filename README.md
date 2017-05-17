# js-stats
Package provides the implementation of various statistics distribution such as normal distribution, fisher, student-t, and so on

[![Build Status](https://travis-ci.org/chen0040/js-stats.svg?branch=master)](https://travis-ci.org/chen0040/js-stats) [![Coverage Status](https://coveralls.io/repos/github/chen0040/js-stats/badge.svg?branch=master)](https://coveralls.io/github/chen0040/js-stats?branch=master) 

# Features

* Normal Distribution
  
  - cumulativeProbability(Z)
  - invCumulativeProbability(p)

* Student's T Distribution

  - cumulativeProbability(t_df)
  - invCumulativeProbability(p)

* Fisher–Snedecor Distribution

  - cumulativeProbabiliyt(F)

* Chi-Square Distribution

  - cumulativeProbabiliy(ChiSquare)

# Install

Run the following npm command to install

```bash
npm install js-stats
```

# Usage

Sample code is available at [playground](https://runkit.com/chen0040/js-stats-playground)

### Using with nodejs

```javascript
jsstats = require('js-stats');

//====================NORMAL DISTRIBUTION====================//

var mu = 0.0; // mean
var sd = 1.0; // standard deviation
var normal_distribution = new jsstats.NormalDistribution(mu, sd);

var X = 10.0; // point estimate value 
var p = normal_distribution.cumulativeProbability(X); // cumulative probability

var p = 0.7; // cumulative probability
var X = normal_distribution.invCumulativeProbability(p); // point estimate value

//====================T DISTRIBUTION====================//

var df = 10; // degrees of freedom for t-distribution
var t_distribution = new jsstats.TDistribution(df);

var t_df = 10.0; // point estimate or test statistic
var p = t_distribution.cumulativeProbability(t_df); // cumulative probability

var p = 0.7;
var t_df = t_distribution.invCumulativeProbability(p); // point estimate or test statistic


//====================F DISTRIBUTION====================//

var df1 = 10; // degrees of freedom for f-distribution
var df2 = 20; // degrees of freedom for f-distribution
var f_distribution = new jsstats.FDistribution(df1, df2);

var F = 10.0; // point estimate or test statistic
var p = f_distribution.cumulativeProbability(F); // cumulative probability


//====================Chi Square DISTRIBUTION====================//

var df = 10; // degrees of freedom for cs-distribution
var cs_distribution = new jsstats.ChiSquareDistribution(df);

var X = 10.0; // point estimate or test statistic
var p = cs_distribution.cumulativeProbability(X); // cumulative probability



```

### Using with HTML page

