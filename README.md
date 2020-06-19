# kma-pipeline

A (somewhat) easy to use microbial pipeline

For now it only supports paired end reads.

It uses trim galore and unicycler to assemble bacterial genomes.

These genomes are indexed in a database.

The database then interfaces with nullarbor in order to compare the assemblies of different microbial isolates.

## Functionality
This microbial genome pipeline executes a number of tasks on each set of reads.

For each set of reads:
 - (QC reads) (FastQC)
 - trim reads (trim galore)
 - species identification on reads (Kraken2)
 - _de novo_ genome assembly (unicycler)
 - (species identification on assembly) (Kraken2
 
For each group of isolates:
 - ([Nullarbor](https://github.com/tseemann/nullarbor)) report


## Development
This project is currently under activate development.
The development is broken down into the following phases which will be completed in listed order:

### Phase 1
 - [x] parse reads
 - [x] merge reads
 - [ ] QC on reads
 
 
### Phase 2
 - [ ] interface with [Nullarbor](https://github.com/tseemann/nullarbor) somehow
 - [ ] use a shared conda environment

### Phase 3
 - [ ] complete database of all isolates
 - [ ] make it searchable
 - [ ] make a nice user-interface (web-gui)
 
### Phase 4
 - [ ] interface with [Bifrost](https://github.com/ssi-dk/bifrost)
 
