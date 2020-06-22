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

Then it generates a pipeline for each sample found among the reads in those paths.

TODO: Disable generation of pipeline samples where the full report is completed.
TODO: Also, make the report system.

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
        parse = line.strip().split('\t') # TODO: only strip newlines. An "empty" tab in the end should not be a parsing error, even though it doesn't make sense to not have a path..
        
        prefix = parse[0]
        method = parse[1] # not implemented yet
        path = parse[2]
        # TODO: Check whether a slash has already been set into the path given in reads_paths_tab col 3 (the path = parse[2] variable.)
        
        if len(parse) != 3:
            raise Exception("Error: Make sure that the input file is correctly tab-delimited on all rows")
        reads_paths_parsed[prefix] = path




print("These are the paths where sample-names will be extracted from.")
print(reads_paths_parsed)
print()




# Now we will iterate through the files in each path.
# We assume that all reads have a constant length suffix
# We assume that the reads have the lane number followed by the direction (lanes_first)

# By eating the suffixes from the end until there are (number_of_files)/4/2 unique sample names, we can extract the sample names.
# A number of small sanity checks are performed underway.

for prefix, path in reads_paths_parsed.items():
    check_set = set() # With this set, I'm checking that each each file in the path is used once only.

    glob_ = sorted(glob.glob(f"{path}/*.fastq.gz"))
    glob_basenames = [os.path.basename(i) for i in glob_]

    n = len(glob_basenames)
    print('n', n)
    print()

    max_name_len = max([len(i) for i in glob_basenames])


    def parallel_end_eating(glob_basenames):
        
        min_length = 5
        
        oldset = set()
        newset = set()
        has_been_correct = False
        global suffix_length
        
        for i in range(1, max_name_len-(min_length-1)): # TODO: En alternativ måde at implementere dette på: Tjek fra starten af strengene, og så stop når min_length er oversteget og alle symboler er ens. Jeg ved ikke om det ville være hurtigere på den måde.
            oldset = newset
            newset = set([j[:-i] for j in glob_basenames])
         

            if len(newset) == n/4/2:
                has_been_correct = True

            if len(newset) < n/4/2:
                if has_been_correct:
                    suffix_length = i-1
                    print('oldset used')
                    return oldset
                else:  
                    raise Exception('The assumptions on the file names are not met. Or there is a bug.')

            
        if has_been_correct:
            suffix_length = i
            print('newset used')
            return newset


    sample_names = sorted(parallel_end_eating(glob_basenames))
    print('this is after the end eating, and the suffix length was', suffix_length)

    print(f"These are the sample_names for prefix '{prefix}':")
    print(sample_names)
    #if not input('Continue? [y]/n ')[0].strip().upper() == 'Y':
    #    print(' user exited...')
    #    exit()
    print()


    # TODO: Make a mechanism that warns if duplicates exist.

    
    for sample_name in sample_names:

        #if not sample_name.startswith("Axx_A"):
        #    continue



        full_name = prefix + '_' + sample_name
        print('Generating jobs for', full_name, '...')

        reads = sorted([i for i in glob_basenames if i.startswith(sample_name) and len(i) == len(sample_name) + suffix_length]) # assuming lanes first
        for i in reads:
            check_set.add(i)
        #print(reads)

        reads_forward = reads[::2]
        reads_reverse = reads[1::2]

        reads_forward_full = [path + i for i in reads_forward]
        reads_reverse_full = [path + i for i in reads_reverse]


        for i, j in zip(reads_forward, reads_reverse):
            print(f"  {i}\t{j}")


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
            walltime = '01:00:00',
            account = 'clinicalmicrobio') << f"""
                mkdir -p output/isolates/{full_name}/cat_reads/

            
                cat_log="output/isolates/{full_name}/cat_reads/cat.log"

                
                
                echo -e "forward:{nl}{nl.join(reads_forward_full)}" > $cat_log
                cat {' '.join(reads_forward_full)} > output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz

                echo -e "{nl}reverse:{nl}{nl.join(reads_reverse_full)}" >> $cat_log
                cat {' '.join(reads_reverse_full)} > output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz


                echo "{nl}{nl}Completed $(date)" >> $cat_log


                """



        gwf.target(sanify('_1_trim_', full_name),
            inputs = [f"output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz",
                      f"output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz"],
            outputs = [f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R1_val_1_fastqc.zip",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2_fastqc.zip"],
            cores = 4,
            memory = '4gb',
            walltime = '04:00:00',
            account = 'clinicalmicrobio') << f"""

                trim_galore --paired --cores 4 --gzip --fastqc -o output/isolates/{full_name}/trim_reads/ output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz  

                # TODO: check if trim galore supports multiple cores now? 
                # TODO: Find a way to save the fastq results? They should be calculated within trim_galore

                """
                


        kraken_reads_top_command = """awk -F '\\t' '$4 ~ "(^S$)|(U)" {gsub(/^[ \\t]+/, "", $6); printf("%6.2f%%\\t%s\\n", $1, $6)}'"""
        gwf.target(sanify('_2_kraken_', full_name),
            inputs = [f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz"],
            outputs = [f"output/isolates/{full_name}/kraken2/kraken2_reads_report.txt",
                       f"output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab"],
            cores = 8,
            memory = '8gb',
            account = 'clinicalmicrobio') << f"""

                mkdir -p output/isolates/{full_name}/kraken2
                
                #kraken2 --db /project/ClinicalMicrobio/faststorage/database/minikraken2_v2_8GB_201904_UPDATE/ --report output/isolates/{full_name}/kraken2/kraken2_reads_report_old.txt --paired output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz > /dev/null
                kraken2 --db /project/ClinicalMicrobio/faststorage/database/minikraken_8GB_20200312/ --report output/isolates/{full_name}/kraken2/kraken2_reads_report.txt --paired output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz > /dev/null

                #cat output/isolates/{full_name}/kraken2/kraken2_reads_report_old.txt | {kraken_reads_top_command} | sort -gr | head -n 10 > output/isolates/{full_name}/kraken2/kraken2_reads_top10_old.tab
                cat output/isolates/{full_name}/kraken2/kraken2_reads_report.txt | {kraken_reads_top_command} | sort -gr | head -n 10 > output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab


                """


        gwf.target(sanify('_3_assemble_', full_name),
            inputs = [f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz"],
            outputs = [f"output/isolates/{full_name}/unicycler/assembly.fasta",
                       f"output/isolates/{full_name}/unicycler/{full_name}_assembly.fasta",
                       f"output/isolates/{full_name}/unicycler/assembly-stats.tab"],
            cores = 8,
            memory = '128gb',
            walltime = '2-00:00:00',
            account = 'clinicalmicrobio') << f"""

                mkdir -p output/isolates/{full_name}/unicycler

                unicycler --min_fasta_length 500 -1 output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz -2 output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz -o output/isolates/{full_name}/unicycler

                cp output/isolates/{full_name}/unicycler/assembly.fasta output/isolates/{full_name}/unicycler/{full_name}_assembly.fasta


                assembly-stats -t output/isolates/{full_name}/unicycler/assembly.fasta > output/isolates/{full_name}/unicycler/assembly-stats.tab


                """
                
            

        gwf.target(sanify('_4_prokka', full_name),
            inputs  = [f"output/isolates/{full_name}/unicycler/assembly.fasta"],
            outputs  = [f"output/isolates/{full_name}/prokka/{full_name}.gff"],
            cores = 8,
            memory = '4g',
            walltime = '04:00:00',
            account = 'clinicalmicrobio') << f"""
                mkdir -p output/isolates/{full_name}/prokka

                prokka --cpu 8 --force --prefix {full_name} --outdir output/isolates/{full_name}/prokka output/isolates/{full_name}/unicycler/{full_name}_assembly.fasta
                
                

                """





        gwf.target(sanify('_5_report', full_name),
            inputs = [

                      f"output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab",

                      f"output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz",
                      f"output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz",
                      
                      f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz",

                      f"output/isolates/{full_name}/unicycler/assembly.fasta",
                      f"output/isolates/{full_name}/prokka/{full_name}.gff"
                      ],
            outputs = [f"output/isolates/{full_name}/report/report.txt"],
            cores = 1,
            memory = '2g',
            walltime = '01:00:00',
            account = 'clinicalmicrobio') << f"""
                # Run a series of python, bash and R-scripts in order to create a report that outlines the results for each sample.
                # The main results should be printed to the database.tab file.

                sleep $[ ( $RANDOM % 60 )  + 1 ]s

                # TODO: generate report: Do this first, because if it fails, nothing should be written to the database.tab later in this spec.
                mkdir -p output/isolates/{full_name}/report
                touch output/isolates/{full_name}/report/report.txt


                # Write to database.tab
                kraken2=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{$1 = ""; print $0;}}' | sed -e 's/^[[:space:]]*//')
                kraken2_p=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{print $1}}' | sed -e 's/%//')

                cat_R1="output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz"
                cat_R2="output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz"

                trim_R1="output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz"
                trim_R2="output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz"

                unicycler_assembly="output/isolates/{full_name}/unicycler/assembly.fasta"

                unicycler_sum=$(tail -n 1 output/isolates/{full_name}/unicycler/assembly-stats.tab | awk '{{print $2}}')
                unicycler_ncontigs=$(tail -n 1 output/isolates/{full_name}/unicycler/assembly-stats.tab | awk '{{print $3}}')
                unicycler_longest=$(tail -n 1 output/isolates/{full_name}/unicycler/assembly-stats.tab | awk '{{print $5}}')


                prokka_gff="output/isolates/{full_name}/prokka/{full_name}.gff"
                prokka_CDS=$(cat output/isolates/{full_name}/prokka/{full_name}.txt | awk '$1 == "CDS:" {{print $0}}' | awk '{{print $2}}')


                # full_name sample_name tech    kraken2_p   kraken2 cat_reads   trim_reads  unicycler_assembly   unicycler_ncontigs unicycler_sum   unicycler_longest        prokka_gff  prokka_CDS date

                echo -e "{full_name}\t{sample_name}\tPE4\t$kraken2_p\t$kraken2\t[\\"$cat_R1\\", \\"$cat_R2\\"]\t[\\"$trim_R1\\", \\"$trim_R2\\"]\t$unicycler_assembly\t$unicycler_ncontigs\t$unicycler_sum\t$unicycler_longest\t$prokka_gff\t$prokka_CDS\t$(date +%F_%H-%M-%S)\t{prefix}\t{path}" >> database.tab


                """

    #    break
    #break


    # Sanity check on the number of reads available and used.
    if len(check_set) != n:
        raise Exception(f"(Fatal) The number of reads available ({n}) in the path ({path}) is not equal to the number of reads used ({len(check_set)}).")        




