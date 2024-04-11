#!/bin/bash

echo $(ps -ef | grep sh | grep show_user_usage.sh)
tmp=$(ps -ef | grep sh | grep show_user_usage.sh |awk NR==1'{print $2}')
echo $tmp
$(kill -9 $tmp)

