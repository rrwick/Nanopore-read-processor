# This tool is deprecated!

It was written back in 2016, and much has since changed in the world of Oxford Nanopore basecalling and read processing. I don't use this script anymore, and you probably shouldn't either!

If you're still interested, the original README follows below:


# Nanopore read analysis

This repo contains a script I wrote to handle Oxford Nanopore files on our servers at the Centre for Systems Genomics at the University of Melbourne. They were only ever intended for us, so if anybody else wants to make use of them, they'll need to do a bit of modification!



### Basics

This script assumes that all Nanopore reads will be deposited in this folder structure on Happyfeet:
`/home/UNIMELB/inouye-hpc-sa/nanopore-data/fast5/TYPE/SAMPLE_NAME`

`TYPE` is either `1d` or `2d` and `SAMPLE_NAME` can be any string.

To use the script, go here: `/home/UNIMELB/inouye-hpc-sa/nanopore-data`
and then run this: `./process_reads.py --help`

The script will look for all `*.fast5` files in the appropriate directories, print out some information about what it finds and then carry out the commands issued by the user.

Two arguments are required:
* `--command` tells the script what to do (details on each command below), and you can use `--command all` to execute all commands.
* `--samples` tells the script which samples to process. It will do a partial match on what's given, so if you use `--samples Kleb` it will run on all samples with `Kleb` in their name. You can use `--samples all` to run it on all found samples.



### Command: list

All this does is list information about the samples found - it doesn't actually do anything to the files.

Example: `./process_reads.py --command list --sample all`



### Command: sort

This command moves fast5 files with no basecalling into a `no_basecall` directory. If it finds fast5 files _with_ basecalling that aren't in a directory that makes sense (like `pass`, `fail` or `basecalled`), then it will move them into a `basecalled` directory.

Example: `./process_reads.py --command sort --sample all`

The end result is that the sample directory will not contain any fast5 files directory. Rather, all fast5 files will be in one of the following directories: `pass`, `fail`, `basecalled` or `no_basecall`.



### Command basecall

This command will go through all the reads without basecalling (those in the `no_basecall` directory) and run [Nanonet](https://github.com/nanoporetech/nanonet) on them. It uses the command `nanonetcall` for 1d reads and `nanonet2d` for 2d reads. This can take a while to finish, especially if there are lot of reads lacking basecalling.

It will by default use 80 threads (all of the threads on Happyfeet) but they are run with maximum niceness, so hopefully won't upset other Happyfeet users too much. You can use the `--nanonet_threads` option to reduce that number, if you want.

Examples:
* `./process_reads.py --command basecall --sample all`
* `./process_reads.py --command basecall --sample all --nanonet_threads 20`

When finished, it will move reads that it successfully basecalled into a `nanonet` directory. This is often about half of the reads previously lacking basecalling, but some will be pretty short and not that useful. However, some will be long and high identity, so it's still good to do this step!



### Command fastq

This command extracts the sequences from the fast5 files to make these files:
* `SAMPLE_NAME.fastq.gz`
* `SAMPLE_NAME.fasta.gz`
* `SAMPLE_NAME.tsv`
These will all be deposited in this directory: `/home/UNIMELB/inouye-hpc-sa/nanopore-data/fastq/SAMPLE_NAME`

This step will look for a `fasta` file of Illumina contigs from the same isolate in this directory: `/home/UNIMELB/inouye-hpc-sa/nanopore-data/references`. If it finds a file (the `fasta` file must contain the sample name), then this command will align all of the reads to the contigs. This step can take a while. As with basecalling, it will use all 80 threads on Happyfeet (but they will be nice) and you can adjust this with the `--alignment_threads` option.

If Illumina contigs for the sample are _not_ available, then this command doesn't do anything too interesting: it will extract the reads to `fasta` and `fastq` and make a `tsv` file with this columns: `Filename, Sample name, Library type, Run name, Flowcell ID Channel number, Basecalling, Read name, Read type, Length Mean qscore, zlib compression ratio, Quality group`.

`zlib compression ratio` is a test to catch repetitive nonsense. Some Nanopore reads end up looking like this: `AATAATAATAATAATAATAATAATAATAAT...` for 10+ kb. To catch these, the script compresses the sequence with zlib to see how much smaller it gets. If it can compress _a lot_ (e.g. to less than 15% of its original size), the read is almost certainly junk and is thrown out.

`Quality group` is how the script categorises the read. If there aren't Illumina contigs for this sample, then this is either 'bad' (when the read is too short or failed the zlib test) or 'unknown'.

If Illumina contigs for the sample _are_ available, then this command can do more. It uses [unicycler_align](https://github.com/rrwick/Unicycler#unicycler-align) to align the reads to the contigs and will include a couple more columns in the tsv: `Alignment identity, Alignment reference`. It then scores each read using both its alignment identity and length to categorise them. Now instead of getting a `Quality group` of 'unknown', reads will be classed as 'poor', 'good' or 'very good'.

'Good' reads are decently long and/or decently high identity. 'Very good' reads are both very long and/or very high identity. Then when saving the files, it makes `SAMPLE_NAME_good.fastq.gz` and `SAMPLE_NAME_very_good.fastq.gz` in addition to the other output files. Note that these categories are not mutually exclusive - the `good` file also includes the `very good` reads and the regular `fastq` file contains the 'good' and 'very good' reads.

You can then choose which reads to use for downstream stuff, like assembly. E.g. if the sequencing run did well and there are a bunch of 'very good' reads, you can use only those. Otherwise, you can use the 'good' reads. Reads that didn't make the good category are likely dodgy, so only go for the full set of reads if you're really desperate!



### Command figures

This command will use the information in the `tsv` file to produce a `png` of some plots using ggplot: read length histograms and stuff like that. If there are Illumina contigs, then it makes more plots showing the breakdown of quality groups. The plots will be saved in the same directory as the `fastq` files: `/home/UNIMELB/inouye-hpc-sa/nanopore-data/fastq/SAMPLE_NAME`



### Command tarball

This command tarballs up all of the fast5 files for the sample and puts it here: `/home/UNIMELB/inouye-hpc-sa/nanopore-data/fast5-gz/SAMPLE_NAME.tar.gz`

This step can take a while, especially if there are a lot of reads.



### Command scp

This final step uses `scp` to transfer the tarball to Helix. When it's done, it will move the sample's `fast5` files to this directory: `/home/UNIMELB/inouye-hpc-sa/nanopore-data/fast5_processed/TYPE/SAMPLE_NAME`

Now that the fast5s have been moved out of the normal directory, they will be excluded from further runs of this script. If, however, you want to include these finished samples, you can run the script with `--include_processed`. Then will look in both `fast5` and `fast5_processed` for read files.



### Command: all

This will do all of the above steps, in order!

Example: `./process_reads.py --command all --sample all`



### License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
