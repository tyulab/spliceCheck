# [Splice Check](http://aso-amenable.ue.r.appspot.com/)
(If test version is broken, please use the above link instead.)
<hr>

12/7 - fixed a minor typo in one of the helper functions that affects the output (had < instead of >)

To use this application you simply enter a gene name, the cDNA coordinates of a mutation in that gene and
the WT and Mutant bases then hit submit.
e.g. CDKL5 c.2152 G>A
Note: gene names are case sensitive

The output page will give you links to information about the gene
and mutation as well as a prediction of whether or not it could be treated with antisense oligonucleotide
gene therapy.
<hr>

# [Current testing version](http://splicecheck.us-east-1.elasticbeanstalk.com/)

#### Recently added:
- Pass transcript ID into list search
- Change version number in config file
- [download the log](https://splicecheck.us-east-1.elasticbeanstalk.com/log.csv) and set max lines in log
#### Features to add:
- Save output to database
    - What type of database? How will it be accessed? Where can it be stored?
- Serverless hosting?
#### To test:
- Transcript ids
- determine_amenability code fixes
- log.csv concurrency
#### Issues:
- Internal Server Error
    - issue with too many requests at once overloading the instances
    - sometimes causes things in static folder to not display
- Gevent monkey patch
    - use when local
- `The system cannot find the path specified.`
    - not sure where anything calls a path, but it has been an error since the beginning
