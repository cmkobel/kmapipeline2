#! /usr/bin/env python3

from gwf import *
import glob
import os 
from carltools import sanify


gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
})





"""

This file reads the list of paths defined in the input_paths_file 

Then it generates a pipeline for each sample found among the reads in those paths


"""


input_paths_file = 'reads_paths.tab'

nl = '\n'
std_setup = {'cores': 8,
             'memory': '1gb', 
             'account': 'clinicalmicrobio'}


reads_paths_parsed = {}
with open(input_paths_file, 'r') as reads_paths:
    for line in reads_paths:

        if line[0] in ['#', '\n']:
            continue # Comments and newlines are OK
        parse = line.strip().split('\t')
        
        prefix = parse[0]
        method = parse[1] # not implemented yet
        path = parse[2]
        
        if len(parse) != 3:
            raise Exception("Error: Make sure that the input file is correctly tab-delimited on all rows")
        reads_paths_parsed[prefix] = path




print("These are the paths where sample-names will be extracted from.")
print(reads_paths_parsed)
print()




# Now we will iterate through the files in each path.
# We assume that all reads have a constant length suffix
# We assume that the reads have the lane number followed by the direction (lanes_first)

# By eating the suffixes from the end until there are (number_of_files)/4/2 unique sample names.
# If this can not be achieved, an error will be thrown.

for prefix, path in reads_paths_parsed.items():

    glob = sorted(glob.glob(f"{path}/*.fastq.gz"))
    glob_basenames = [os.path.basename(i) for i in glob]

    n = len(glob_basenames)
    print('n', n)
    print()

    max_name_len = max([len(i) for i in glob_basenames])

    for i in glob_basenames:
        #print(i)
        pass


    def parallel_end_eating(glob_basenames):
        min_length = 5
        
        oldset = set()
        newset = set()
        has_been_correct = False
        
        for i in range(1, max_name_len-(min_length-1)): # TODO: En alternativ måde at implementere dette på: Tjek fra starten af strengene, og så stop når min_length er oversteget og alle symboler er ens. Jeg ved ikke om det ville være hurtigere på den måde.
            oldset = newset
            newset = set([j[:-i] for j in glob_basenames])
         

            if len(newset) == n/4/2:
                has_been_correct = True

            if len(newset) < n/4/2:
                if has_been_correct:
                    return oldset
                else:  
                    raise Exception('The assumptions on the file names are not met. Or there is a bug.')

            
        if has_been_correct:
            return newset

    sample_names = sorted(parallel_end_eating(glob_basenames))

    print(f"These are the sample_names for prefix '{prefix}':")
    print(sample_names)
    #if not input('Continue? [y]/n ')[0].strip().upper() == 'Y':
    #    print(' user exited...')
    #    exit()
    print()


    # TODO: Make a mechanism that warns if duplicates exist.

    for sample_name in sample_names:
        full_name = prefix + '_' + sample_name
        print('Generating jobs for', full_name, '...')

        reads = sorted([i for i in glob_basenames if i.startswith(sample_name)]) # assuming lanes first
        #print(reads)

        reads_forward = reads[::2]
        reads_reverse = reads[1::2]

        reads_forward_full = [path + i for i in reads_forward]
        reads_reverse_full = [path + i for i in reads_reverse]


        for i, j in zip(reads_forward, reads_reverse):
            print(f"\t{i}\t{j}")


        print()

        # The problem is that it asks also when you write gwf status
        #if not input('Continue? [y]/n ')[0].strip().upper() == 'Y':
        #    print(' user exited...')
        #    exit()
        
        gwf.target(sanify('_0_cat_reads_', full_name),
            inputs = reads_forward_full + reads_reverse_full,
            outputs = [f"output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz",
                       f"output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz"],
            cores = 8,
            memory = '1gb',
            account = 'clinicalmicrobio') << f"""
                mkdir -p output/isolates/{full_name}/cat_reads/

            
                cat_log="output/isolates/{full_name}/cat_reads/cat.log"

                echo "" > $cat_log
                
                echo -e "forward: {' '.join(reads_forward_full)}" >> $cat_log
                cat {' '.join(reads_forward_full)} > output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz

                echo -e "forward: {' '.join(reads_reverse_full)}" >> $cat_log
                cat {' '.join(reads_reverse_full)} > output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz


                echo Completed $(date) >> $cat_log


                """



        gwf.target(sanify('_1_trim_', full_name),
            inputs = [f"output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz",
                      f"output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz"],
            outputs = [f"output/isolates/{full_name}/cat_reads/PE_R1_trimmed.fq.gz",
                      f"output/isolates/{full_name}/cat_reads/PE_R2_trimmed.fq.gz"],
            cores = 1,
            memory = '1gb',
            account = 'clinicalmicrobio') << f"""

                trim_galore --paired --gzip -o output/isolates/{full_name}/trim_reads/ output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz  

                # TODO: check if trim galore supports multiple cores now? 
                # TODO: Find a way to save the fastq results? They should be calculated within trim_galore

                """




        break





    



    # Now we will shorten the suffixes until we get n/4/2 unique sample names




