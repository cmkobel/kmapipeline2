#! /usr/bin/env python3
from gwf import *
import glob
import os 
from cutils import *
import json



gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
})


print("""
                                      _____ _         _ _         ___ 
                                     |  _  |_|___ ___| |_|___ ___|_  |
                                     |   __| | . | -_| | |   | -_|  _|
                                     |__|  |_|  _|___|_|_|_|_|___|___|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|_|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""")

"""

This file reads the list of paths defined in the input_paths_file 

Then it generates a pipeline for each sample found among the reads in those paths.

TODO: Disable generation of pipeline samples where the full report is completed.
TODO: Also, make the report system.

"""


input_paths_file = 'reads_paths.tab'
input_paths_done_file = 'other/paths_done.tab'

nl = '\n'
std_setup = {'cores': 8,
             'memory': '1gb', 
             'account': 'clinicalmicrobio'}


# TODO: Make a mechanism that warns if duplicates exist.
# This will be implemented using the paths_done.tab file.
# could be moved to a script? To keep the pipeline lean.

def parse_paths_done(input_paths_done_file):
    paths_done = []
    with open(input_paths_done_file) as input_paths_done_file_open:
        #first_line = True # skip header line
        for line in input_paths_done_file_open:
            #if first_line == True:
            #    first_line = False
            #    continue
            if line[0] in ['#', '\n']: # comments and newlines are welcome.
                continue
            path = line.split('\t')[0]
            paths_done.append(path)

    return paths_done


paths_done = parse_paths_done(input_paths_done_file)

# TODO: implement reading the blacklist
#       Do that next time something fails.


# Find the paths to parse:
# Underway, we will check whether the path has already been completely processed.

# TODO: functionize this, so the variables can't bleed over into other parts of the program.
def parse_input(input_paths_file):
    reads_paths_parsed = {}
    with open(input_paths_file, 'r') as reads_paths:
        for _i, line in enumerate(reads_paths):

            if line[0] in ['#', '\n']: # If the first symbol is a hash or newline.
                continue # Comments and newlines are OK
            parse = line.strip('\n').split('\t') 
            
            
            try:
                prefix = parse[0]
                singular_sample_name = parse[1]
                path = parse[2]
                method = parse[3] # better name: technology or tech
                
                if len(parse) > 4:
                    lanes_first = parse[4]
                else:
                    lanes_first = 'yes'
                    dprint(f"No information given about lanes first for {prefix}. Assuming lanes first.")



            except Exception as e:
                print(f"error at line {_i+1}:\n {line}")
                print(e, '\nMake sure there are four tab-delimited columns in the', input_paths_file, 'file.')
                exit() # necessary?
            
            # lanes_first switching is not implemented yet. But the parsing is done:
            if lanes_first == 'no': 
                lanes_first_bool = False
            else:
                lanes_first_bool = True # Defaults to Truef


            # TODO: Check whether a slash has already been set into the path given in reads_paths_tab col 3 (the path = parse[2] variable.)
            if path[-1] != "/":
                path += '/'


            if " " in singular_sample_name:
                raise exception('Please remove any spaces in singular_sample_name.') # TODO: test this


            reads_paths_parsed[prefix] = {'method': method, 'path': path, 'singular_sample_name': singular_sample_name, 'lanes_first': lanes_first_bool}

    return reads_paths_parsed

reads_paths_parsed = parse_input(input_paths_file)



def overview_presentation():
    print("These are the paths where sample-names will be extracted from.")
    for key, dict_ in reads_paths_parsed.items():

        # shorten the path if it is too long to show on screen
        if len(dict_['path']) > 90:
            shortened_path = dict_['path'][:11] + '…' + dict_['path'][-69:]
        else:
            shortened_path = dict_['path']


        if dict_['path'] in paths_done:
            status_msg = ' complete:   will be skipped: '
        else:
            status_msg = '> incomplete: will be enqueued:'
        print(f"{status_msg} \"{key}\" {dict_['singular_sample_name']} @ {shortened_path} ({dict_['method']})")
    print()

overview_presentation()



# Now we will iterate through the files in each path.
# We assume that all reads have a constant length suffix
# We assume that the reads have the lane number followed by the direction (lanes_first)

# By eating the suffixes from the end until there are (number_of_files)/4/2 unique sample names, we can extract the sample names.
# A number of small sanity checks are performed underway.




print('Now iterating over the remaining paths, if any:')


for prefix, dict_ in reads_paths_parsed.items():
    
    # Skip paths which are done.
    if dict_['path'] in paths_done:
        continue

    if dict_['method'] == "PE4":
        PE_method = 4
    elif dict_['method'] == "PE2": # Warning; has not been tested yet.
        PE_method = 2
    elif dict_['method'] == "PE1": # Attention, only tested on 200626
        PE_method = 1
    else:
        raise Exception(f"Fatal: method {dict_['method']} is not yet implemented.")


    if dict_['singular_sample_name'] != "":
        singular_sample_name_bool = True
    else:
        singular_sample_name_bool = False

    print()
    print()
    print()
    print(f"#                   {dict_['path']}                   #") 
    print(f"{(40 + len(dict_['path']))*'~'}")
    print(f"prefix: \"{prefix}\"")
    print('attributes:', dict_)

    check_set = set() # With this set, I'm checking that each each file in the dict_['path'] is used once only.

    # Man skal passe på med at tro at der nødvendigvis er en bug. Nogengange er der faktisk ikke. Måske er det noget helt andet der gør at ens kode ikke virker.
    glob_ = sorted(glob.glob(f"{dict_['path']}/*.fastq.gz"))
    #print('glob', glob_) # Bliver der globbet noget?
    glob_basenames = [os.path.basename(i) for i in glob_]

    
    n = len(glob_basenames)
    print(f"number of reads files: {n}")
    print()



    if n == 0:
        raise Exception('Fatal: glob_basenames is empty. Have you specified the correct path?')
    #print('these are the glob-basenames', glob_basenames) # I wonder why it is empty when there is only one sample.
    min_name_len = min([len(i) for i in glob_basenames])

    global suffix_length 

    # This is the new one
    def parallel_right_align_eating(glob_basenames, PE_method):
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

            for set_item in set_0:
                dprint('', set_item)
            dprint()
            # Check that the number of samples is correct (n/4/2):
            if len(set_0) == n/PE_method/2:
                has_been_correct = True
                
                # Debug prints:
                if not True:
                    dprint(len(set_0))
                    for k in set_0:
                        dprint('', k)
                    print()

                #continue # No need to do more in this iteration
                
            # Check that the last 2 letters in all suffixes are matching
            if has_been_correct and check2(set_0):
                dprint()
                dprint('this should be it ')
                suffix_length = max_name_len - i + 2
                dprint('the suffix length is given at i', i, 'as', suffix_length)
                return set_2
        raise Exception('for loop exhausted')



    #sample_names = sorted(parallel_end_eating(glob_basenames))

    # If there is only one sample, the pattern can not be automatically parsed.
    if singular_sample_name_bool:
        sample_names = [dict_['singular_sample_name']]
        suffix_length = len(glob_basenames[0]) - len(dict_['singular_sample_name'])

        print('sn-debug: the singular_sample_name has been finally set')


        # Check that the pe-number is correct

        if len(glob_basenames)/PE_method/2 != 1:
            raise Exception(f"The PE number ({dict_['method']}) seems to be incorrect as the number of files ({len(glob_basenames)}) is not equal to the expected ({PE_method}*2 = {PE_method*2}).")


    # if (len(glob_basenames)/PE_method/2 == 1.0):
    #     with open(dict_['path'] + '/pipeline2.input.txt') as manual_text_file:
    #         manual_sampleid = [line.strip() for line in manual_text_file][0]
    #     print('loaded manual sampleid:', manual_sampleid)
        
    #     suffix_length = len(manual_sampleid) - len(glob_basenames[0]) -1

    #     sample_names = [manual_sampleid]
    else:
        sample_names = sorted(parallel_right_align_eating(glob_basenames, PE_method))
        print('sn-debug: the sample names have been set with parallel eating')
    #print('sn', sample_names)
    # I think they have to be sorted, because of the way the grouping into lanes and directions is done

    # This cannot be printed when there is only one sample.
    #dprint('this is after the end eating, and the suffix length was', suffix_length)
    

    print(f"These are the {len(sample_names)} sample_name(s) for prefix '{prefix}': ({n}/{len(sample_names)} = {n/len(sample_names)} files per isolate <fits PE{int(n/len(sample_names)/2)}>)")
    for sn in sample_names:
        print(f"  {sn}")
    #if not input('Continue? [y]/n ')[0].strip().upper() == 'Y':
    #    print(' user exited...')
    #    exit()
    print()
    print()



    


    
    for i_, sample_name in enumerate(sample_names):

        #if not sample_name.startswith("Axx_A"):
        #    continue


        full_name = prefix + '_' + sample_name

        
        reads = sorted([i for i in glob_basenames if i.startswith(sample_name) and len(i) == len(sample_name) + suffix_length]) # assuming lanes first
        for i in reads:
            check_set.add(i)
        #print(reads)


        reads_forward = reads[::2]
        reads_reverse = reads[1::2]
        

        reads_forward_full = [dict_['path'] + i for i in reads_forward]
        reads_reverse_full = [dict_['path'] + i for i in reads_reverse]
        
        #print(len(reads))
        if len(reads) == 0:
            raise Exception(f"Fatal: No reads have matched the sample_name. Below is a comparison of (1) the given sample name and (2) the first globbed file for this path:\
                \n\t(1) {sample_name}\
                \n\t(2) {glob_basenames[0]}\
                \nPlease correct the sample name such that it matches the globbed files.")

        spacer = (len(reads_forward[0]) - len(sample_name)-5) * ' ' # -5 is for the length of "(R1)"
        print(f" {i_+1}) Generating jobs for {full_name} ...\n  >{sample_name} (R1){spacer} >{sample_name} (R2)")

        # TODO: An idea for a sanity check: Add the number of forward+reverse reads for each sample to a set. Check tha the size of the set is 1.


        for i, j in zip(reads_forward, reads_reverse):
            print(f"   {i}  {j}")


        print()

        # The problem is that it asks also when you write gwf status
        #if not input('Continue? [y]/n ')[0].strip().upper() == 'Y':
        #    print(' user exited...')
        #    exit()

        # Skip blacklisted samples. Some samples just don't comply.
        #blacklist = ['HA_101', '16-138518_S10', '16-149463_S7', 'D10_00037', 'D_02357', 'HA_11']

        # New, better blacklist with full_name's instead of sample_name's
        blacklist = ['200622_HA_11',
                     '180511_D10_00037',
                     '180828_D10_00037',
                     #'180914_D_02357',
                     '200622_HA_38', '200622_HA_94', '200622_HA_101']
        #blacklist = []
        if prefix + '_' +sample_name in blacklist:
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
            walltime = '16:00:00',
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
            memory = '32gb',
            walltime = '8:00:00',
            account = 'clinicalmicrobio') << f"""

                mkdir -p output/isolates/{full_name}/kraken2
                

                kraken2 --threads 8 --db /project/ClinicalMicrobio/faststorage/database/minikraken_8GB_20200312/ --report output/isolates/{full_name}/kraken2/kraken2_reads_report.txt --paired output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz > /dev/null

                
                cat output/isolates/{full_name}/kraken2/kraken2_reads_report.txt | {kraken_reads_top_command} | sort -gr | head -n 10 > output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab
                
                # Write to the kraken-database
                # TODO: Manually update the old samples
                kraken2=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{$1 = ""; print $0;}}' | sed -e 's/^[[:space:]]*//')
                kraken2_p=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{print $1}}' | sed -e 's/%//')
                echo -e "{full_name}\t$kraken2_p\t$kraken2" >> kraken2.tab




                """



        gwf.target(sanify('_2.1_coverg_', full_name),
            inputs = [f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz",
                      f"output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab"],
            outputs = [f"output/isolates/{full_name}/coverage/coverage_unfiltered.tab"],
            cores = 8,
            memory = '16g', #'16g',
            walltime = '2:00:00', #'4:00:00',
            account = 'clinicalmicrobio') << f"""
                echo $(pwd)
                
                # I failed to get sambamba to work in the p222 environment, so I use a different one.
                source /home/cmkobel/miniconda3/etc/profile.d/conda.sh 
                conda activate antihum # I ought to change its name to "coverage"


                # Pick the reference automatically
                kraken2=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{$1 = ""; print $0;}}' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed -e "s/ /_/g")
                echo "The species name is ${{kraken2}}"

                if [[ "$kraken2" == "unclassified" ]]; then
                    echo "Because the top hit is unclassified, the coverage will not be measured."
                    touch output/isolates/{full_name}/coverage/coverage_unfiltered.tab # disable from workflow
                    exit 1
                fi

                if [[ $(compgen -G "references/${{kraken2}}/*.fa") ]]; then 
                    echo reference found
                    reference=$(compgen -G "references/${{kraken2}}/*.fa" | head -n 1)
                    echo $reference

                    reference_basename=$(basename $reference)
                    reference_basename_stem=${{reference_basename%.*}}

                    echo "$reference_basename $reference_basename_stem"


                else 
                    echo reference not found
                    echo $kraken2 >> other/missing_references.tab
                    exit 1
                fi 



                # Most of the rest is taken from anti_hum 


                
                mkdir -p output/isolates/{full_name}/coverage
                cd output/isolates/{full_name}/coverage
                
                # Clear old reference
                touch reference.fa
                rm reference.*

                # Copy newest reference and index
                cp ../../../../$reference ${{reference_basename_stem}}.fa
                echo "$reference" >> reference_source.txt
                echo "indexing..."
                bwa index -a bwtsw ${{reference_basename_stem}}.fa
                



                # Map
                echo mapping...
                bwa mem -M -t 8 ${{reference_basename_stem}}.fa ../trim_reads/PE_R1_val_1.fq.gz ../trim_reads/PE_R2_val_2.fq.gz  \
                | sambamba view -f bam -F "proper_pair" -S -t 8 /dev/stdin > unsorted.bam
                #| sambamba view -f bam -T ${{reference_basename_stem}}.fa /dev/stdin > unsorted.bam
                
                  
                # Sort
                echo sorting...
                sambamba sort -t 8 -m 16GB --out=sorted.bam --tmpdir=/scratch/$GWF_JOBID/ unsorted.bam
                sambamba index sorted.bam 

                rm unsorted.bam 


                # Get coverage of the non-filtered mapping.
                sambamba depth window -w 1000 sorted.bam > coverage_unfiltered.tab

                echo -e "{full_name}\t$kraken2\t$reference_basename_stem\treads\t$(awk '{{ total += $5; count++ }} END {{ print total/count }}' coverage_unfiltered.tab)" >> coverage_report.tab

                echo done


                # is filtering relevant?
                

                """



        gwf.target(sanify('_3_unicyc_', full_name),
            inputs = [f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz",
                      f"output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab"], # Not strictly necessary, just makes the database so much more integrated.
            outputs = [f"output/isolates/{full_name}/unicycler/assembly.fasta",
                       f"output/isolates/{full_name}/unicycler/final_assembly/{full_name}.fasta",
                       f"output/isolates/{full_name}/unicycler/assembly-stats.tab"],
            cores = 16,
            memory = '256g', 
            walltime = '12:00:00',
            account = 'clinicalmicrobio') << f"""

                # Somehow, samtools (dependency of unicycler) has a problem together with some other package in p222. I don't know which. What a mess.
                source /home/cmkobel/miniconda3/etc/profile.d/conda.sh 
                conda activate unicycler 

                # Create a report before doing anything
                # If unicycler fails, then there will at least be a report which indicates that the assembly has failed.

                mkdir -p output/isolates/{full_name}/report

                # Collect data for report
                kraken2=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{$1 = ""; print $0;}}' | sed -e 's/^[[:space:]]*//')
                kraken2_p=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{print $1}}' | sed -e 's/%//')

                cat_R1="output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz"
                cat_R2="output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz"

                trim_R1="output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz"
                trim_R2="output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz"

                unicycler_assembly="" 

                unicycler_sum=""
                unicycler_ncontigs=""
                unicycler_longest=""

                prokka_gff=""
                prokka_CDS=""

                method='{dict_['method']}'


                # full_name sample_name tech kraken2_p kraken2 cat_reads trim_reads unicycler unicycler_ncontigs unicycler_sum unicycler_longest prokka_gff prokka_CDS pipeline_date prefix path
                echo -e "{full_name}\t{sample_name}\t$method\t$kraken2_p\t$kraken2\t[\\"$cat_R1\\", \\"$cat_R2\\"]\t[\\"$trim_R1\\", \\"$trim_R2\\"]\t$unicycler_assembly\t$unicycler_ncontigs\t$unicycler_sum\t$unicycler_longest\t$prokka_gff\t$prokka_CDS\t$(date +%F_%H-%M-%S)\t{prefix}\t{dict_['path']}" > output/isolates/{full_name}/report/meta_report.txt






                # temp
                # cp output/isolates/{full_name}/report/meta_report.txt output/isolates/{full_name}/report/check_meta_report.txt
                # exit 0






                mkdir -p output/isolates/{full_name}/unicycler
                mkdir -p other # for errors


                unicycler --threads 16 --min_fasta_length 500 -1 output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz -2 output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz -o output/isolates/{full_name}/unicycler || echo -e "{full_name}\t{sample_name}\t{dict_['path']}\t$(date +%F_%H-%M-%S)" >> other/blacklist.tab

                # a file with a name might be easier to work with.     
                # This is also the file that p2assemblyextractor works with.           
                mkdir -p output/isolates/{full_name}/unicycler/final_assembly
                cp output/isolates/{full_name}/unicycler/assembly.fasta output/isolates/{full_name}/unicycler/final_assembly/{full_name}.fasta


                assembly-stats -t output/isolates/{full_name}/unicycler/assembly.fasta > output/isolates/{full_name}/unicycler/assembly-stats.tab



                # It would make some sense to update the meta report here and put in the assembly stats.
                # But I think it is not worth the debugging effort, so the user will have to wait for the summary job to run (after prokka).




                """
                

		# Alternative to unicycler: skesa
        gwf.target(sanify('_3_skesa__', full_name),
            inputs = [f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz",
                      f"output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab"], # Not strictly necessary, just makes the database so much more integrated.
            outputs = [f"output/isolates/{full_name}/skesa/{full_name}.fa",
                       f"output/isolates/{full_name}/skesa/assembly-stats.tab"],
            cores = 16,
            memory = '128g', 
            walltime = '12:00:00',
            account = 'clinicalmicrobio') << f"""

           		mkdir -p output/isolates/{full_name}/skesa
           		cd output/isolates/{full_name}/skesa

           		skesa --reads ../trim_reads/PE_R1_val_1.fq.gz,../trim_reads/PE_R2_val_2.fq.gz --cores 16 --memory 128 --min_contig 500 > {full_name}.fa

           		assembly-stats -t {full_name}.fa > assembly-stats.tab

           		"""


            
        gwf.target(sanify('_4_prokka_', full_name),
            inputs  = [f"output/isolates/{full_name}/unicycler/assembly.fasta"],
            outputs  = [f"output/isolates/{full_name}/prokka/{full_name}.gff"],
            cores = 16,
            memory = '8g',
            walltime = '04:00:00',
            account = 'clinicalmicrobio') << f"""
                


                mkdir -p output/isolates/{full_name}/prokka

                # prokka --cpu 8 --force --prefix {full_name} --outdir output/isolates/{full_name}/prokka output/isolates/{full_name}/unicycler/final_assembly/{full_name}.fasta









                # Copied from assemblycomparator

                


                cd output/isolates/{full_name}/prokka #
            
                # In reality, it is never possible to create and identical assembly, so the if-statement doesn't really need to check for an existing hash-dir.
                # But in order to make the code more portable, i wish to keep it as it is..

                # hash tables
                # Generate hash key from assembly
                hash=$(cat ../unicycler/final_assembly/{full_name}.fasta | sha256sum | awk '{{print $1}}')
                # Set up directories
                hash_base="/faststorage/project/ClinicalMicrobio/database/prokka_hash/keys"
                hash_key_dir="${{hash_base}}/${{hash}}"
                # Check if the key exists
                if [[ -d "${{hash_key_dir}}" ]]; then
                    echo "Hash key exists"
                    echo $hash
                    #mkdir -p prokka
                    
                    # log usage
                    echo -e "copying from ${{hash}}" > prokka_hash.txt
                    echo -e "copy\t$(pwd)/\t${{hash}}\t{full_name}\t$(date +%F_%H-%M-%S)" > "${{hash_key_dir}}"/usage_log.tab
                    
                    cp "${{hash_key_dir}}/prokka."* .
                    # Touch it all to update the modified dates
                    touch prokka.*
                    
                    echo -e "${{hash}}" > hash.txt
                else
                    
                    prokka --cpu 16 --force --outdir . --prefix prokka ../unicycler/final_assembly/{full_name}.fasta 
                    if [[ -d "${{hash_base}}" ]]; then
                        mkdir -p "${{hash_key_dir}}"
                        cp prokka.* "${{hash_key_dir}}"
                        echo -e "{full_name}\t${{hash}}\tsha256sum\tpipeline2\t$(date +%F_%H-%M-%S)" >> ${{hash_key_dir}}/index.tab
                        touch ${{hash_key_dir}}/{full_name}.sample_name
                        cat ${{hash_key_dir}}/index.tab >> ${{hash_base}}/index.tab
                    fi
                fi
                cp prokka.gff {full_name}.gff
                
                

                """



        gwf.target(sanify('_5_report_', full_name),
            inputs = [

                      f"output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab",

                      f"output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz",
                      f"output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz",
                      
                      f"output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz",
                      f"output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz",

                      f"output/isolates/{full_name}/unicycler/final_assembly/{full_name}.fasta",
                      f"output/isolates/{full_name}/prokka/{full_name}.gff"
                      ],
            outputs = [f"output/isolates/{full_name}/report/meta_report.txt"],
            cores = 1,
            memory = '1g',
            walltime = '00:10:00',
            account = 'clinicalmicrobio') << f"""
                # Run a series of python, bash and R-scripts in order to create a report that outlines the results for each sample.
                # The main results should be printed to the database.tab file.

                #sleep $[ ( $RANDOM % 2060 )  + 1 ]s




                # TODO: generate report: Do this first, because if it fails, nothing should be written to the database.tab later in this spec.
                mkdir -p output/isolates/{full_name}/report

                # Collect data for report
                kraken2=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{$1 = ""; print $0;}}' | sed -e 's/^[[:space:]]*//')
                kraken2_p=$(head -n 1 output/isolates/{full_name}/kraken2/kraken2_reads_top10.tab | awk '{{print $1}}' | sed -e 's/%//')

                cat_R1="output/isolates/{full_name}/cat_reads/PE_R1.fastq.gz"
                cat_R2="output/isolates/{full_name}/cat_reads/PE_R2.fastq.gz"

                trim_R1="output/isolates/{full_name}/trim_reads/PE_R1_val_1.fq.gz"
                trim_R2="output/isolates/{full_name}/trim_reads/PE_R2_val_2.fq.gz"

                unicycler_assembly="output/isolates/{full_name}/unicycler/final_assembly/{full_name}.fasta"

                unicycler_sum=$(tail -n 1 output/isolates/{full_name}/unicycler/assembly-stats.tab | awk '{{print $2}}')
                unicycler_ncontigs=$(tail -n 1 output/isolates/{full_name}/unicycler/assembly-stats.tab | awk '{{print $3}}')
                unicycler_longest=$(tail -n 1 output/isolates/{full_name}/unicycler/assembly-stats.tab | awk '{{print $5}}')


                prokka_gff="output/isolates/{full_name}/prokka/{full_name}.gff"
                prokka_CDS=$(cat output/isolates/{full_name}/prokka/prokka.txt | awk '$1 == "CDS:" {{print $0}}' | awk '{{print $2}}')

                method='{dict_['method']}'

                # full_name sample_name tech    kraken2_p   kraken2 cat_reads   trim_reads  unicycler_assembly   unicycler_ncontigs unicycler_sum   unicycler_longest        prokka_gff  prokka_CDS date
                echo -e "{full_name}\t{sample_name}\t$method\t$kraken2_p\t$kraken2\t[\\"$cat_R1\\", \\"$cat_R2\\"]\t[\\"$trim_R1\\", \\"$trim_R2\\"]\t$unicycler_assembly\t$unicycler_ncontigs\t$unicycler_sum\t$unicycler_longest\t$prokka_gff\t$prokka_CDS\t$(date +%F_%H-%M-%S)\t{prefix}\t{dict_['path']}" > output/isolates/{full_name}/report/meta_report.txt



                """


    # Kill dict_['path']
    # When all isolates from a dict_['path'] has been correctly processed, the dict_['path'] is added to a file (paths_done.tab). This disables the dict_['path'] from being processed in the pipeline.
    if True: # Can be easily disabled for debugging.
        input_list = [f"output/isolates/{prefix + '_' + sample_name}/report/meta_report.txt" for sample_name in sample_names if prefix + '_' + sample_name not in blacklist]
        gwf.target(sanify('_99_kill__', prefix),
            inputs = input_list,
            outputs = [],
            cores = 1,
                memory = '1g',
                walltime = '01:00:00',
                account = 'clinicalmicrobio') << f"""



                ./update_db.sh



                # Exclude the path from the pipeline when all are complete.
                mkdir -p other
                echo -e "{dict_['path']}\t{prefix}\t$(date +%F_%H-%M-%S)" >> other/paths_done.tab

                """

    #    break
    #break


    # Sanity check on the number of reads available and used.
    if len(check_set) != n:
        raise Exception(f"(Fatal) The number of reads available ({n}) in the path ({dict_['path']}) is not equal to the number of reads used ({len(check_set)}).")        


print() # a nice newline between the reads-files and the gwf jobs.

#TODO: count to 8 when adding reads for catting. 
       #Well, that is actually a bad idea: It will make the pipeline less flexible. 
       #Let's say you want to add samples which are not PE4 but PE1 ?
