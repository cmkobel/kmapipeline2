#! /usr/bin/env python3

from gwf import *
import glob
import os 
from carltools import *



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


# TODO: Make a mechanism that warns if duplicates exist.
# This will be implemented using the paths_done.tab file.
# could be moved to a script? To keep the pipeline lean.
paths_done = []
with open('other/paths_done.tab') as paths_done_open:
    first_line = True
    for line in paths_done_open:
        if first_line == True:
            first_line = False
            continue
        path = line.split('\t')[0]
        paths_done.append(path)


# Find the paths to parse:
# Underway, we will check whether the path has already been completely processed.
reads_paths_parsed = {}
with open(input_paths_file, 'r') as reads_paths:
    for line in reads_paths:

        if line[0] in ['#', '\n']:
            continue # Comments and newlines are OK
        parse = line.strip('\n').split('\t') # TODO: only strip newlines. An "empty" tab in the end should not be a parsing error, even though it doesn't make sense to not have a path..
        

        dprint('the len', len(parse))
        try:
            prefix = parse[0]
            method = parse[1] # not implemented yet
            path = parse[2]
            lanes_first = parse[3]
        except Exception as e:
            print(e, '\nMake sure there are four tab-delimited columns in the', input_paths_file, 'file.')
            exit() # necessary?
        # TODO: Check whether a slash has already been set into the path given in reads_paths_tab col 3 (the path = parse[2] variable.)
        
        # lanes_first switching is not implemented yet.
        if lanes_first == 'no': 
            lanes_first_bool = False
        else:
            lanes_first_bool = True # Defaults to Truef

        
        #reads_paths_parsed[prefix] = path 
        # deprecated

        reads_paths_parsed[prefix] = {'method': method, 'path': path, 'lanes_first': lanes_first_bool}



print("These are the paths where sample-names will be extracted from.")
for key, dict_ in reads_paths_parsed.items():
    if dict_['path'] in paths_done:
        status_msg = 'complete:  will be skipped  :'
    else:
        status_msg = 'incomplete: will be queued  :'
    print(f"  {status_msg}  \"{key}\" @ {dict_['path']}")
print()




# Now we will iterate through the files in each path.
# We assume that all reads have a constant length suffix
# We assume that the reads have the lane number followed by the direction (lanes_first)

# By eating the suffixes from the end until there are (number_of_files)/4/2 unique sample names, we can extract the sample names.
# A number of small sanity checks are performed underway.





print('Now iterating over the remaining paths:')



for prefix, dict_ in reads_paths_parsed.items():
    dprint('items', prefix, dict_)

    if dict_['path'] in paths_done:
        continue

    print()
    print()
    print()
    print(f"#                   {dict_['path']}                   #") 
    print(f"{(40 + len(dict_['path']))*'~'}")
    print(f"prefix: \"{prefix}\"")

    check_set = set() # With this set, I'm checking that each each file in the dict_['path'] is used once only.

    glob_ = sorted(glob.glob(f"{dict_['path']}/*.fastq.gz"))
    glob_basenames = [os.path.basename(i) for i in glob_]

    n = len(glob_basenames)
    print(f"number of files: {n}")
    print()

    #max_name_len = max([len(i) for i in glob_basenames])



    def parallel_end_eating(glob_basenames):
        # deprecated
        
        min_length = 5
        
        oldset = set()
        newset = set()
        has_been_correct = False
        global suffix_length # We need this much later in the program.
        
        for i in range(1, max_name_len-(min_length-1)): # TODO: En alternativ måde at implementere dette på: Tjek fra starten af strengene, og så stop når min_length er oversteget og alle symboler er ens. Jeg ved ikke om det ville være hurtigere på den måde.
            oldset = newset
            newset = set([j[:-i] for j in glob_basenames])
         

            if len(newset) == n/4/2:
                has_been_correct = True

            if len(newset) < n/4/2:
                if has_been_correct:
                    suffix_length = i-1 # The suffix length is used later, to check that the correct sample names and files are linked.
                    dprint('oldset used')
                    return oldset
                else:  
                    raise Exception('The assumptions on the file names are not met. Or there is a bug.')

            
        if has_been_correct:
            suffix_length = i
            dprint('newset used')
            return newset

    min_name_len = min([len(i) for i in glob_basenames])

    


    # This is the new one
    def parallel_right_align_eating(glob_basenames):
        """ This in an improvement of the old 'parallel_end_eating()'.
            'parallel_right_align_eating()' eats from the start, relative to a common length. 
            When n/4/2 unique sample names have been found, the parallel eating continues until a common character occurs in _all_ samples.
            Then, the sample names are found.
            This, I think is the most robust way of inferring the sample names. The old function had problems recognizing the complete sample names if they had wildly different endings.
        """

        global suffix_length # The suffix length is used later, to check that the correct sample names and files are linked.

        n = len(glob_basenames)
        min_name_len = min([len(i) for i in glob_basenames])
        max_name_len = max([len(i) for i in glob_basenames])
        diff = max_name_len - min_name_len

        # Først tjekker has_been_correct at antallet passer med n/4/2
        # Derefter kigges der efter 2 ens bogstaver i alle.

        # This routine checks that the last 2 letters are the same in all strings in the set.
        def check2(set_):
            set_l = list(set_)

            for j in set_l:
                if j[-2:] != set_l[0][-2:]:
                    return False

            return True

            
        set_0 = set()
        set_1 = set()
        set_2 = set()

        has_been_correct = False # This is a one way switch.

        for i in range(diff+1, max_name_len+1):

            # shift the saved sets:
            set_2 = set_1
            set_1 = set_0
            set_0 = set()
            
            set_0 = set(j[:len(j) - max_name_len  + i] for j in glob_basenames)

            # Check that the number of samples is correct (n/4/2):
            if len(set_0) == n/4/2:
                has_been_correct = True
                
                # Debug prints:
                if not True:
                    dprint(len(set_0))
                    for k in set_0:
                        dprint('', k)
                    print()

                #continue # No need to do more in this iteration
                
            # Check that the last 2 letters in all suffixes is matching
            if has_been_correct and check2(set_0):
                dprint()
                dprint('this should be it ')
                suffix_length = max_name_len - i + 2
                dprint('the suffix length is given at i', i, 'as', suffix_length)
                return set_2



    #sample_names = sorted(parallel_end_eating(glob_basenames))
    sample_names = sorted(parallel_right_align_eating(glob_basenames))
    # I think they have to be sorted, because of the way the grouping into lanes and directions is done

    dprint('this is after the end eating, and the suffix length was', suffix_length)

    print(f"These are the {len(sample_names)} sample_names for prefix '{prefix}': ({n}/{len(sample_names)} = {n/len(sample_names)})")
    for sn in sample_names:
        print(f"  {sn}")
    #if not input('Continue? [y]/n ')[0].strip().upper() == 'Y':
    #    print(' user exited...')
    #    exit()
    print()
    print()



    


    
    for sample_name in sample_names:

        #if not sample_name.startswith("Axx_A"):
        #    continue


        full_name = prefix + '_' + sample_name
        print(' Generating jobs for', full_name, '...')

        reads = sorted([i for i in glob_basenames if i.startswith(sample_name) and len(i) == len(sample_name) + suffix_length]) # assuming lanes first
        for i in reads:
            check_set.add(i)
        #print(reads)


        reads_forward = reads[::2]
        reads_reverse = reads[1::2]

        reads_forward_full = [dict_['path'] + i for i in reads_forward]
        reads_reverse_full = [dict_['path'] + i for i in reads_reverse]

        # TODO: An idea for a sanity check: Add the number of forward+reverse reads for each sample to a set. Check tha the size of the set is 1.


        for i, j in zip(reads_forward, reads_reverse):
            print(f"   {i}  {j}")


        print()

        # The problem is that it asks also when you write gwf status
        #if not input('Continue? [y]/n ')[0].strip().upper() == 'Y':
        #    print(' user exited...')
        #    exit()

        # Skip blacklisted samples. Some samples just don't comply.
        blacklist = ['HA_101']
        if sample_name in blacklist:
            continue
        
        gwf.target(sanify('_0_catrds_', full_name),
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



        gwf.target(sanify('_1_trim___', full_name),
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



        gwf.target(sanify('_3_asembl_', full_name),
            inputs = [f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz"],
            outputs = [f"output/isolates/{full_name}/unicycler/assembly.fasta",
                       f"output/isolates/{full_name}/unicycler/{full_name}_assembly.fasta",
                       f"output/isolates/{full_name}/unicycler/assembly-stats.tab"],
            cores = 4,
            memory = '64gb', 
            walltime = '2-00:00:00',
            account = 'clinicalmicrobio') << f"""

                mkdir -p output/isolates/{full_name}/unicycler

                unicycler --min_fasta_length 500 -1 output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz -2 output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz -o output/isolates/{full_name}/unicycler

                cp output/isolates/{full_name}/unicycler/assembly.fasta output/isolates/{full_name}/unicycler/{full_name}_assembly.fasta


                assembly-stats -t output/isolates/{full_name}/unicycler/assembly.fasta > output/isolates/{full_name}/unicycler/assembly-stats.tab


                """
                

            
        gwf.target(sanify('_4_prokka_', full_name),
            inputs  = [f"output/isolates/{full_name}/unicycler/assembly.fasta"],
            outputs  = [f"output/isolates/{full_name}/prokka/{full_name}.gff"],
            cores = 8,
            memory = '4g',
            walltime = '04:00:00',
            account = 'clinicalmicrobio') << f"""
                mkdir -p output/isolates/{full_name}/prokka

                prokka --cpu 8 --force --prefix {full_name} --outdir output/isolates/{full_name}/prokka output/isolates/{full_name}/unicycler/{full_name}_assembly.fasta
                
                

                """



        gwf.target(sanify('_5_report_', full_name),
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

                mkdir -p database/backup/
                cp database.tab database/backup/database_backup_$(date +%F_%H-%M-%S).tab

                # full_name sample_name tech    kraken2_p   kraken2 cat_reads   trim_reads  unicycler_assembly   unicycler_ncontigs unicycler_sum   unicycler_longest        prokka_gff  prokka_CDS date
                echo -e "{full_name}\t{sample_name}\tPE4\t$kraken2_p\t$kraken2\t[\\"$cat_R1\\", \\"$cat_R2\\"]\t[\\"$trim_R1\\", \\"$trim_R2\\"]\t$unicycler_assembly\t$unicycler_ncontigs\t$unicycler_sum\t$unicycler_longest\t$prokka_gff\t$prokka_CDS\t$(date +%F_%H-%M-%S)\t{prefix}\t{dict_['path']}" >> database.tab


                """


    # Kill dict_['path']
    # When all isolates from a dict_['path'] has been correctly processed, the dict_['path'] is added to a file (paths_done.tab). This disables the dict_['path'] from being processed in the pipeline.
    if True: # Can be easily disabled for debugging.
        gwf.target(sanify('_99_kill__', full_name),
            inputs = [f"output/isolates/{prefix + '_' + sample_name}/report/report.txt" for sample_name in sample_names if sample_name not in blacklist],
            outputs = [],
            cores = 1,
                memory = '1g',
                walltime = '00:10:00',
                account = 'clinicalmicrobio') << f"""

                mkdir -p other
                echo -e "{dict_['path']}\t{prefix}\t$(date +%F_%H-%M-%S)" >> other/paths_done.tab

                """

    #    break
    #break


    # Sanity check on the number of reads available and used.
    if len(check_set) != n:
        raise Exception(f"(Fatal) The number of reads available ({n}) in the path ({dict_['path']}) is not equal to the number of reads used ({len(check_set)}).")        




#TODO: count to 8 when adding reads for catting. 
       #Well, that is actually a bad idea: It will make the pipeline less flexible. 
       #Let's say you want to add samples which are not PE4 but PE1 ?
