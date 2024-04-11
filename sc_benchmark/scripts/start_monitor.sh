#!/bin/bash
while getopts n:t: opt
do 
	case "${opt}" in
		n) name=${OPTARG};;
		t) time=${OPTARG};;
	esac
done

for (( i=0; i<1; i++ ))
do
    {
	    bash ./sc_benchmark/scripts/show_user_usage.sh -n $name -t $time
    } &  #将上述程序块放到后台执行
done
#wait    #等待上述程序结束
echo "开始监控！"
