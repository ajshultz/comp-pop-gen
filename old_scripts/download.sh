#!/bin/bash

list=`cat accessions.lst`
njobs=24
path_to_ascp="/n/home13/sayshwarya/.aspera/connect/bin/ascp"
path_to_ascp_openssh="/n/home13/sayshwarya/.aspera/connect/etc/asperaweb_id_dsa.openssh"

for i in $list ; do
 echo -n "--> STARTING ASSESSION $i DOWNLOAD @ " ; date
 /n/holylfs/EXTERNAL_REPOS/SOFTWARE/sratoolkit.2.5.5-centos_linux64/bin/prefetch --force yes --max-size 500000000 -a "${path_to_ascp}|${path_to_ascp_openssh}" --ascp-options "-QT -l 10G" $i ./ > $i.log & 
 running=`jobs -r | grep -c Running`
 while [ "$running" -gt "$njobs" ] ; do
        sleep 1
        running=`jobs -r | grep -c Running`
        done
done
