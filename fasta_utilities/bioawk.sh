#!/bin/bash
/nfs/esnitkin/bin_group/bioawk-master/bioawk -c fastx '{ print $name, length($seq) }' < $1
