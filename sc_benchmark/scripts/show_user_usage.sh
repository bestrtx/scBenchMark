#!/bin/bash

while getopts n:t: opt
do 
	case "${opt}" in
		n) name=${OPTARG};;
		t) time=${OPTARG};;
	esac
done
echo "NameVal: $name";

# Replace USERNAME with the username of the user you want to monitor
USERNAME="ljr"


totle_time=0
while true
do
    # Get the memory usage and CPU usage of the user
    MEMORY_USAGE=$(ps -u $USERNAME -o rss | awk '{sum+=$1} END {print sum/1024}')
    CPU_USAGE=$(ps -u $USERNAME -o %cpu | awk '{sum+=$1} END {print sum}')

    # Get the current date and time
    DATE=$(date +"%Y-%m-%d %H:%M:%S")

    # Write the memory usage and CPU usage to a file
    echo "$DATE Memory Usage: $MEMORY_USAGE MB CPU Usage: $CPU_USAGE%" >> $name.log

    # Wait for 1 minute before checking again
    sleep $time
    totle_time=`expr $totle_time + $time`
    if [ $totle_time -gt 43200 ]
      then
        break
      fi
done

