<img src="scripts/logo.png" alt="Pipeline2" width="320"/>

This is an automated pipeline for use for many microbial isolates.

Quality control is followed by species detection, assembly and finally comparisons in groups of isolates: Core genome analysis etc.

A simple database is maintained with relevant metadata for the isolates. This database should interface with nullarbor, tormes, bifrost etc. 

## Functionality
This microbial genome pipeline executes a number of tasks on each set of reads.

For each set of reads:
 - QC and trimming reads (trim galore)
 - species identification on reads and assemblies (Kraken2)
 - _de novo_ genome assembly (unicycler)
 - (species identification on assembly) 
 
For each group of isolates:
 - ([Nullarbor](https://github.com/tseemann/nullarbor)) report


## Development
This project is currently under activate development.
The development is broken down into the following phases which will be completed in listed order:

### Phase 1
 - [x] parse reads
 - [x] concatenate reads 
 - [x] QC on reads (fastqc)
 - [x] species detection (kraken2)
 - [x] assembly (unicycler)
 - [x] annotation (prokka)
 - [x] simple database
 
 
### Phase 2
 - [ ] interface with [Nullarbor](https://github.com/tseemann/nullarbor)
 - [ ] interface with [Bifrost](https://github.com/ssi-dk/bifrost)
 - [ ] use a shared conda environment

### Phase 3
 - [ ] complete database of all isolates
 - [ ] make it searchable
 - [ ] make a nice user-interface (web-gui)
 
 
 ## Usage
 
Using this pipeline is easy. 
 
If you have a directory of PE illumina reads, you simply enter the path in the reads_paths.tab file, which should be present in the root directory with the following header:

`#prefix	singular_sample_name	path	tech	lanes_first	comments...`

Then, check if everything is OK by calling `gwf status`

If everything is as you'd like it, the pipeline can be started with `gwf run`

The output will then be written to output/{prefix}_{sample_name}/


 
 
