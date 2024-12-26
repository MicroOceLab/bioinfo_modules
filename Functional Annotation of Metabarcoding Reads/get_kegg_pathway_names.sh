#!/bin/bash

kos_list=${1}
kos_w_desc=${2}

kegg_pathways=$(curl -# https://rest.kegg.jp/list/pathway)

while read line
do
	hit=$(echo -e "${kegg_pathways}" | grep "${line/ko/map}")
	echo -e "${line}\t${hit}"
	echo -e "${line}\t${hit}" | cut -f1,3 >> ${kos_w_desc}
done < ${kos_list}
