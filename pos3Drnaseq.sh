#!/bin/bash

#  Copyright (c) 2022 Miquéias Fernandes

#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:

#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.

#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.

N_ARGS=$#
if [ $N_ARGS -lt 6 ]
then
 echo "Usage:  $> bash pos3Drnaseq.sh      seu@email    pre_results.zip   3Drnaseq_out.zip  http://...ptnas.faa.gz   http://...genome.gff3.gz temp_dir"
 exit 1
fi

export TZ=America/Sao_Paulo
EMAIL=$(echo $1 | sed s/@/%40/)
OUT_PRE=$2
OUT_3D=$3
PTNAS=$4
GFF=$5
TMP=$6
TIMEOUT=60

preparar () {
    echo "preparando ..."
    apt install curl wget 1>/dev/null 2>/dev/null
    pip install biopython deeptools 1>/dev/null 2>/dev/null
    if ! grep 'pysam.index(bamFile)' /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py 1>/dev/null 2>/dev/null
    then
        cp /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py .
        cp bamHandler.py /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py.old
        grep -B100000 'bam = pysam.Samfile(bamFile' bamHandler.py | grep -v 'bam = pysam.Samfile(bamFile' > xtemp
        echo '        pysam.index(bamFile)' >> xtemp
        tail bamHandler.py -n+`grep -n 'bam = pysam.Samfile(bamFile' bamHandler.py | cut -d: -f1` >> xtemp
        cp xtemp /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py
        rm xtemp bamHandler.py
    fi
    if [ ! -d $TMP ] ; then mkdir $TMP && echo "diretorio criado: $TMP" ; fi
}

importar() {
    echo "importando ..."
    rm tmp tmp_dir -rf && mkdir tmp_dir
    cp $OUT_3D tmp_dir/out3d.zip
    unzip tmp_dir/out3d.zip >/dev/null && cd tmp && mv file* out3d && mv out3d ../tmp_dir && cd .. && rm tmp -r && cd tmp_dir
    cut -d, -f1 out3d/result/Significant\ DAS\ genes\ list\ and\ statistics.csv | tail +2 | uniq > das_genes
    echo "DAS genes: $(wc -l das_genes)"
    mkdir outpre && cd outpre && unzip ../../$OUT_PRE 1>/dev/null && cd ../
    grep -f <(cut -d, -f2 das_genes | tr -d \" | awk '{print ","$0","}') <(sed s/$/,/ outpre/transcript_gene_mapping.csv) > das_transcripts
    echo "DAS transcripts: $(wc -l das_transcripts)"
    paste -d',' das_transcripts <(sed 's/.*cds_//' das_transcripts | cut -d_ -f-2)| tr -s , , > das_ptnas
    wget -qO ptnas.faa.gz $PTNAS && gunzip ptnas.faa.gz

    echo "das_ptnas = 'das_ptnas'" > script.py
    echo "ptnas = 'ptnas.faa'" >> script.py
    cat >> script.py << EOF
from Bio import SeqIO, Seq, SeqRecord
pts = set([x.strip().split(',')[2] for x in open(das_ptnas).readlines()])
ss= [s for s in SeqIO.parse(ptnas, 'fasta') if s.id in pts]
ptna_name = [f"{x.id} {x.description.replace(x.id, '').strip()}\n" for x in ss]
ptna_seq = [f'{x.id},{str(x.seq)}\n' for x in ss]
open('ptna_name', 'w').writelines(ptna_name)
open('ptna_seq', 'w').writelines(ptna_seq)
EOF
    python3 script.py 
    rm script.py
    wget -qO genes.gff3.gz $GFF && gunzip genes.gff3.gz
}

anotar () {
    echo "anotando ..."
    if [ ! -d ../$TMP/anotacoes ] ; then mkdir ../$TMP/anotacoes ; fi
    ## https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=iprscan5
    API='https://www.ebi.ac.uk/Tools/services/rest/iprscan5/'
    Q='goterms=true&pathways=true&appl=PfamA'
    k=1
    TT=$(grep -c , ptna_seq)
    while read l
    do 
        ID=$(echo $l | cut -d, -f1)
        SEQ=$(echo $l | cut -d, -f2)
        if [ ! -f ../$TMP/anotacoes/$ID.tsv ]
        then
            JOB=$(curl -sSX POST --header 'Content-Type: application/x-www-form-urlencoded' --header 'Accept: text/plain' -d "email=$EMAIL&$Q&sequence=$SEQ" $API/run)
            echo "[$k de $TT em $( date +%D.%H:%M:%S)] rodando $ID pelo job $JOB ..."
            sleep 30s

            for i in $(seq $TIMEOUT)
            do
                if grep FINISHED <(curl -sSX GET --header 'Accept: text/plain' "$API/status/$JOB") >/dev/null
                then 
                    curl -sSX GET --header 'Accept: text/tab-separated-values' "$API/result/$JOB/tsv" > ../$TMP/anotacoes/$ID.tsv
                    break
                else sleep 1m
                fi
            done
        else echo "[$k de $TT] $ID restaurado de $TMP/anotacoes ..."
        fi
        (( k=k+1 ))
    done < ptna_seq 
    cp -r ../$TMP/anotacoes .
}

cobertura () {
    if [ ! -d ../$TMP/beds ] ; then mkdir ../$TMP/beds ; fi
    k=1
    TT="$(( $(tail +2 outpre/experimental_design.csv | grep -c , )  * $(grep -vc '^$' das_genes) ))"
    while read g
    do 
        for SAMPLE in `cut -d, -f2 outpre/experimental_design.csv | tail +2`
        do
            if [ ! -f ../$TMP/beds/$g.$SAMPLE.bed ]
            then
                bamCoverage -b outpre/out_$SAMPLE.sorted.bam -o ../$TMP/beds/$g.$SAMPLE.bed --outFileFormat bedgraph --binSize 3 -p 2 -r $g 
            fi 
            echo "[$k de $TT] $g / $SAMPLE  OK ..."
            (( k=k+1 ))
        done
    done < das_genes 
}

main () {
    preparar
    importar
    anotar & cobertura &
    wait
}

main $@