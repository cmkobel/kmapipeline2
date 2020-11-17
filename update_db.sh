#!/bin/bash


# Backup
mkdir -p database/backup/
echo "Backing up:"

# TODO: Warn the user if a backupable file does not exist.

backup() {
	var=$(basename $1)
	touch "$1"
	if [[ ! -f "$1" ]]; then
		touch "$1"
		echo "  Warning: ${1} didn't exist. Empty file created."
	fi
	cp "$1" database/backup/${var}_$(date +%F_%H-%M-%S).tab && echo "  $var OK"
}

#touch collected_database.tab # if it doesn't exist
#cp collected_database.tab database/backup/collected_database_backup_$(date +%F_%H-%M-%S).tab && echo "  collected_database OK"

backup collected_database.tab

backup reads_paths.tab

backup other/paths_done.tab


echo ""
echo "Collecting meta reports:"

# Reset content and write header
a=$(cat collected_database.tab | wc -l)
echo -e "full_name\tsample_name\ttech\tkraken2_p\tkraken2\tcat_reads\ttrim_reads\tunicycler\tunicycler_ncontigs\tunicycler_sum\tunicycler_longest\tprokka_gff\tprokka_CDS\tpipeline_date\tprefix\tpath\tcoverage\tmedian_inssize\tmode_inssize\tmad_inssize" > collected_database.tab

# Collect

j=0
N=$(ls output/isolates/*/report/meta_report.txt | wc -l)
printf "\r  $j    / $N"
for i in output/isolates/*/report/meta_report.txt; do
	if [[ $(($j %3)) == 0 ]]; then
		printf "\r  $j"
	fi
	cat $i >> collected_database.tab
	j=$((j+1))
done
printf "\r  $j"

b=$(cat collected_database.tab | wc -l)

# Finally
#printf " OK\n"
printf "\n"
#echo "  All meta reports have been collected into collected_database.tab"
echo "$((b-a)) new lines added to collected_database.tab"
