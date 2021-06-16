Splice Check

This project was inspired by my work in the Genetics department at Boston Children's Hospital. Last year my
lab designed a personalized medicine for a little girl with a rare disease because we found a rare splicing
mutation in her genome that could be targeted with this particular type of gene therapy called antisense
oligonucleotides (ASOs). We approximate that only about 10% of patient's mutations can be fixed with this
type of treatment, but every day we are getting emails from parents with sick children asking for help. It is
my job to check if these missense mutations are causing mis-splicing in a way that could be targeted with ASOs.
This usually requires me using a bunch of different tools and resouces so that's why I created splice check to
pull all of those into one place.

Splice Check starts by taking the user's input and using it to call the Variant Effect Predictor (VEP). VEP is
a super useful tool because it returns some of the scores that I'm interested in evaluating as well as other
nomenclature for the variant. This output is then parsed to store the SIFT and PolyPhen scores (two different
algorithms for predicting the effect a mutation has on protein function), as well as get the genomic coordinates
because that is the form of input required for the next tool.

Next, it takes the genomic coordinates for a variant (the chromosome number and location on that chromosome) and
pulls a small piece of genomic sequence surrounding the mutation from a database at UCSC. This string of nucleotides
is the input needed to run MaxEntScan. MaxEntScan is a tool developed by a lab at MIT to predict the likelihood that
there is splicing occurring in a particular piece of DNA.

Once all the tools are run I have a function that looks at the scores and checks if it can be easily determined
whether or not this mutation is a good candidate for ASO therapy. All of the results are returned to an output
page as well as some links to other genomic databases so you can get more information if needed.