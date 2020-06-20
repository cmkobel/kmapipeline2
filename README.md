# kmapipeline2

This is an automated pipeline that takes a list of paths with microbial isolates (PE).

Quality control is followed by species detection, assembly and finally comparisons in groups of isolates: Core genome analysis etc.

A simple database is maintained with relevant metadata for the isolates. This database can interface with nullarbor, tormes, bifrost etc. 

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
 
