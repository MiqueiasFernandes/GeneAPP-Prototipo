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
if [ $N_ARGS -lt 4 ]
then
 echo "Usage:  $> bash pos3Drnaseq.sh      seu@email    pre_results.zip   3Drnaseq_out.zip  http://...ptnas.faa.gz   http://...genome.gff3.gz temp_dir"
 exit 1
fi

EMAIL=$(echo $1 | sed s/@/%40/)
OUT_PRE=$2
OUT_3D=$3
PTNAS=$4
GFF=$4
TMP=$5

preparar () {
    echo "preparando ..."
    apt install curl wget 1>/dev/null 2>/dev/null
    echo "email: $EMAIL"
    echo "tmp: $TMP"
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
    wget -qO ptnas.faa.gz $PTNAS && gunzip ptnas.faa.gz
    wget -qO genes.gff3.gz $GFF && gunzip genes.gff3.gz
}

anotar () {
    echo "anotando ..."
    ## https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=iprscan5
    API='https://www.ebi.ac.uk/Tools/services/rest/iprscan5/'
    SEQ=$(echo $1 | tr "[:lower:]" "[:upper:]" | tr -cd "[:alpha:]")
    JOB=$(curl -X POST --header 'Content-Type: application/x-www-form-urlencoded' --header 'Accept: text/plain' -d "email=$EMAIL&sequence=$SEQ" $API/run 1> log 2> err)
    echo rodando o job $JOB ...
    ## status: curl -X GET --header 'Accept: text/plain' 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/iprscan5-R20220430-195319-0285-44991048-p1m'
    ## result: !curl -X GET --header 'Accept: text/tab-separated-values' 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/iprscan5-R20220430-195319-0285-44991048-p1m/tsv'
}

main () {
    preparar
    importar
}

main $@