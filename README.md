Genotype And Genotype-likelihood concordance scripts...

With sequencing duplicate pairs...


Usage:

(GT concordance. Output one line per pair. Can be re-ordered to 3x3 table)

GTpairConcordance-new.pl inVcf Pair.list > GTconcordance.out



(GL concordance, estimated with Bayes factor. Output one line per site)

PairGLconcordance.pl inVcf Pair.list Exclude.list > GLconcordance.out


Pair list format (delimited by tab):

P11  P12

P21  P22

...


Exclude list format:

Exclude-sample1

Exclude-sample2

...
