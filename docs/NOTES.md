# Notes


## *genomeindex* subcommand timing


```
date && sgains.py genomeindex --genome-pristine hg19_pristine --genome-dir hg19 && date                                                                                                                                                                           

Tue Nov  7 21:16:50 EET 2017
genomeindex subcommand runned with args: Namespace(config=None, dry_run=False, force=False, func=<bound method GenomeIndexCommand.run of <commands.genomeindex_command.GenomeIndexCommand object at 0x10ef8cdd8>>, genome_dir='hg19', genome_index='genomeindex', genome_pristine='hg19_pri
stine', genome_version='hg19', parallel=1, verbose=0)
...
Total time for backward call to driver() for mirror index: 00:54:56
Tue Nov  7 23:06:43 EET 2017
```