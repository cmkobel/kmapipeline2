#!/bin/bash


# Backup
touch collected_database.tab # if it doesn't exist
mkdir -p database/backup/
cp collected_database.tab database/backup/collected_database_backup_$(date +%F_%H-%M-%S).tab && echo "Backup: collected_database OK"
cp reads_paths.tab database/backup/reads_paths_backup_$(date +%F_%H-%M-%S).tab && echo "Backup: reads_paths.tab OK"


# Reset
touch collected_database.tab
rm collected_database.tab

# Collect
j=1
N=$(ls output/isolates/*/report/meta_report.txt | wc -l)
for i in output/isolates/*/report/meta_report.txt; do
	printf "\r$j / $N"
	cat $i >> collected_database.tab
	j=$((j+1))
done

# Finally
printf " OK\n"
echo "All meta reports have been collected into collected_database.tab"