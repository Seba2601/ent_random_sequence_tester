# StringENT test suite — ENT battery revisited

The [Fourmilab Random Sequence Tester](https://www.fourmilab.ch/random/),
**ent**, applies various tests to sequences of bytes stored in files
and reports the results of those tests. The program is useful for
evaluating pseudorandom number generators for encryption and
statistical sampling applications, compression algorithms, and other
applications where the information density of a file is of interest.

The extension of the battery, named **StringENT** [Luengo & al], successfully introduces two additional tests and provides p-values for all of them, which is the most useful way to determine whether the randomness hypothesis holds.

## Description

**StringENT** performs a variety of tests on the stream of bytes in its input
file (or standard input if no input file is specified) and produces
output as follows on the standard output stream:

    StringENT | Results report
    --------------------------------------------------
    Entropy = 7.9999986667 bits per byte.
    Optimum compression would reduce the size
    of this 125000000 byte file by 0 percent.
    --------------------------------------------------
    Chi square distribution for 125000000 samples is 231.02.
    --------------------------------------------------
    Arithmetic mean value of data bytes is 127.5029 (127.5 = random).
    --------------------------------------------------
    Monte Carlo value for Pi is 3.141369650 (error 0.01 percent).
    --------------------------------------------------
    Serial correlation coefficient is -0.000037 (totally uncorrelated = 0.0).
    --------------------------------------------------
    The number of runs test is 62507249 runs.
    --------------------------------------------------
    The local means test's X^2 statistic is 121762.608670 for 122070 blocks.
    --------------------------------------------------

With the **p-value** option enabled, the output is displayed as follows:
    
    StringENT | Results report
    --------------------------------------------------
    Entropy = 7.9999986667 bits per byte.
    Optimum compression would reduce the size
    of this 125000000 byte file by 0 percent.
    --------------------------------------------------
    Chi square distribution for 125000000 samples is 231.02.
    p-value   0.857105
    --------------------------------------------------
    Arithmetic mean value of data bytes is 127.5029 (127.5 = random).
    p-value   0.662362
    --------------------------------------------------
    Monte Carlo value for Pi is 3.141369650 (error 0.01 percent).
    p-value   0.535373
    --------------------------------------------------
    Serial correlation coefficient is -0.000037 (totally uncorrelated = 0.0).
    p-value   0.676014
    --------------------------------------------------
    The number of runs test is 62507249 runs.
    p-value   0.194782
    --------------------------------------------------
    The local means test's X^2 statistic is 121762.608670 for 122070 blocks.
    p-value   0.732796
    --------------------------------------------------

This version of the test suite is designed to be used with byte files.

The values calculated are as follows:

#### Entropy
The information density of the contents of the file, expressed as a
number of bits per character. The results above, which resulted from
processing an image file compressed with JPEG, indicate that the file
is extremely dense in information—essentially random. Hence,
compression of the file is unlikely to reduce its size. By contrast,
the C source code of the program has entropy of about 4.9 bits per
character, indicating that optimal compression of the file would reduce
its size by 38%. \[Hamming, pp. 104–108\]

#### Chi-square Test
The Chi-square test is one of the most established methods for evaluating the randomness of a sequence, and is particularly effective at exposing flaws in pseudorandom generators. It measures how closely the observed frequency of each byte value (from 0 to 255) matches the expected frequency under a perfectly uniform distribution. Specifically, it computes a statistic chi-square by summing, for each byte value, the squared difference between observed and expected counts, divided by the expected count. Under the assumption of randomness, the chi-square follows a Chi-square distribution with 255 degrees of freedom. The p-value derived from this distribution represents the probability of obtaining a deviation at least as large as the one observed. A low p-value indicates that the observed distribution is unlikely to occur by chance, suggesting non-randomness. For the test to be valid, it's generally required that the sequence be long enough to ensure at least five expected occurrences per byte value, ensuring a reliable approximation to the theoretical distribution.
See \[Knuth, pp. 35–40\] for more information on the chi-square test.

#### Arithmetic Mean
This is simply the result of summing the all the bytes (bits if the 
`-b` option is specified) in the file and dividing by the file length. 
If the data are close to random, this should be about 127.5 (0.5 for 
`-b` option output). If the mean departs from this value, the values 
are consistently high or low.
The p-value is computed by measuring how far the observed mean deviates 
from the expected mean under randomness. This deviation is scaled into 
a z-score using the known variance of the uniform distribution. 
The p-value then reflects the probability of observing a result at least
as extreme as this z-score under a standard normal distribution.

#### Monte Carlo Value for Pi
Each successive sequence of six bytes is used as 24 bit X and Y 
co-ordinates within a square.  If the distance of the 
randomly-generated point is less than the radius of a circle inscribed 
within the square, the six-byte sequence is considered a “hit”.  The 
percentage of hits can be used to calculate the value of π.  For very 
large streams (this approximation converges very slowly), the value 
will approach the correct value of π if the sequence is close to 
random.  A 500000 byte file created by radioactive decay yielded:

    Monte Carlo value for Pi is 3.143580574 (error 0.06 percent).

To compute the p-value, we treat each coordinate pair as a random 
point in a square and check whether it falls inside a circle. 
Under the randomness hypothesis, the expected proportion of points 
inside the circle is π/4. We model this as a binomial experiment, 
where each point has a π/4 chance of success. Using the Central Limit 
Theorem, we convert the observed number of points inside the circle 
into a z-score. The p-value is then obtained by measuring how extreme 
this z-score is under a standard normal distribution. A significant 
deviation suggests non-randomness in the spatial distribution of points.

#### Serial Correlation Coefficient
This quantity measures the extent to which each byte in the file 
depends upon the previous byte. For random sequences, this value (which 
can be positive or negative) will, of course, be close to zero.  A 
non-random byte stream such as a C program will yield a serial 
correlation coefficient on the order of 0.5.  Wildly predictable data 
such as uncompressed bitmaps will exhibit serial correlation 
coefficients approaching 1.  See [Knuth, pp. 64–65] for more details.
To compute the p-value for serial correlation, we first measure how 
strongly a sequence is correlated with a version of itself shifted 
by one position.  If the data is random, this correlation should 
be close to zero. The observed correlation is transformed into a 
statistical score that follows a Student t-distribution, assuming 
a large enough sample size. The p-value is then derived by assessing 
how extreme this score is under that distribution. For practical reasons,
a normal approximation is used instead of the exact t-distribution, 
which provides sufficient accuracy for large datasets. A small p-value 
would indicate that the sequence shows more correlation than expected 
by chance, suggesting predictability.

#### Runs Test

The test works by transforming a numerical sequence into a sequence of signs.
Each value is compared to the median: if the value is greater than the median, 
it is assigned a plus sign; if it is smaller, it receives a minus sign. This 
sign sequence is then analyzed for runs—which are uninterrupted sequences of 
identical signs. A new run starts each time the sign changes.
The key idea behind this test is that truly random sequences should exhibit 
a number of runs that falls within a predictable range, based on the number 
of plus and minus signs observed. The expected number of runs, as well as 
its statistical variance and other details, can be found in [NIST: Runs test 
for detecting non-randomness]. Comparing the actual number of runs to this 
expected value allows us to compute a standardized score (Z-value), which 
in turn is used to determine a two-tailed p-value.

#### Local Means Test

This test builds on two core ideas: the use of local arithmetic means within 
the sequence and a classical Chi-square goodness-of-fit analysis.
The procedure begins by dividing the input sequence into equally sized blocks. 
Each block contains a fixed number of bytes (e.g., 1024 by default), and any 
remaining bytes at the end of the sequence that do not fit evenly are discarded. 
For each block, we compute its average value, which is then compared to the 
expected mean for random data.
Because the values within each block should, under the assumption of randomness, 
fluctuate around the expected mean, these block-wise averages can be treated as 
approximately normally distributed. We then assess how well the observed 
distribution of these averages aligns with the theoretical normal distribution 
using a Chi-square goodness-of-fit test.
If the sequence is truly random, the block averages should be relatively uniform 
and centered around the expected value. However, if the sequence contains structure
these local averages can become significantly skewed. The Chi-square statistic will 
then grow accordingly, leading to very small p-values.

## License

This software is licensed under the Creative Commons
Attribution-ShareAlike license.  Please see [LICENSE.md](LICENSE.md) in
this repository for details.

## References

[Hamming]
Hamming, Richard W.  *Coding and Information Theory*. Englewood Cliffs 
NJ: Prentice-Hall, 1980.  ISBN 978-0-13-139139-0.

[Knuth]
Knuth, Donald E.  *The Art of Computer Programming, Volume 2 / 
Seminumerical Algorithms*.  Reading MA: Addison-Wesley, 1969. ISBN 
978-0-201-89684-8.

[Lempel & Ziv]
Ziv J. and A. Lempel.  “A Universal Algorithm for Sequential Data 
Compression”. IEEE Transactions on Information Theory 23, 3, pp. 
337–343.

[Park & Miller]
Park, Stephen K. and Keith W. Miller.  “Random Number Generators: Good 
Ones Are Hard to Find”. Communications of the ACM, October 1988, p. 
1192.

[Luengo & al]
Almaraz Luengo, E. S., Alaña Olivares, B., García Villalba, L. J. et al. «StringENT Test Suite: ENT Battery Revisited for Efficient P Value Computation». Journal of Cryptographic Engineering, vol. 13, n.o 2, junio de 2023, pp. 235-49.

*[Introduction to Probability and
Statistics](https://www.fourmilab.ch/rpkp/experiments/statistics.html)*
at Fourmilab

*[NIST: Runs test for detecting non-randomness](https://www.itl.nist.gov/div898/handbook/eda/section3/eda35d.htm)*
