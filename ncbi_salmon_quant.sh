#!/bin/bash

echo "[1] $(date) verificando pacotes."

if ! command -v salmon &> /dev/null
then
    echo "instalar pacotes!"
    echo "apt install sra-toolkit trimmomatic fastqc salmon"
    exit
fi

tid=t$(date +%s)

echo '[2.1] baixando o GTF ...'
wget -O gene.$tid.gtf.gz $1 $> _2.1_transcripts.download.log
echo '[2.2] descompactando o GTF ...'
gunzip gene.$tid.gtf.gz $> _2.2_transcripts.download.log

echo '[3.1] baixando os transcritos ...'
wget -O cds.$tid.fa.gz $2 $> _3.1_transcripts.download.log
echo '[3.2] descompactando os transcritos ...'
gunzip cds.$tid.fa.gz $> _2.2_transcripts.unzip.log
echo '[3.3] indexando os transcritos ...'
salmon index -t cds.$tid.fa --index idx$tid $> _2.3_transcripts.index.log

for x in $@
    do 
        if [[ $x doesNotContain *","* ]]
        then
         continue
        else
            RUN=`echo $x | cut -d, -f1`
            SAMPLE=`echo $x | cut -d, -f2`
            echo "$RUN: $SAMPLE" >> processed
            i=$( wc -l processed | tr -cs 0-9 , | cut -d, -f2 )

            echo "[4.$i.1] obtendo a amostra $SAMPLE pelo acesso $RUN no sra ..."
            fastq-dump --split-3 $RUN $> _4.1_download.$RUN.$SAMPLE.log
            
            echo "[4.$i.2] fazendo controle de qualidade da amostra $SAMPLE com o TrimmomaticPE ..."
            !TrimmomaticPE \
                $RUN_1.fastq $RUN_2.fastq \
                $SAMPLE.1.fq $SAMPLE.1.unp.fq \
                $SAMPLE.2.fq $SAMPLE.2.unp.fq \
                ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                $> _4.2_qc.$SAMPLE.log
            
            echo "[4.$i.3] reportando controle de qualidade da amostra $SAMPLE com fastqc ..."
            rm qc_$SAMPLE -rf && mkdir qc_$SAMPLE
            fastqc $SAMPLE.1.fq $SAMPLE.2.fq -o qc_$SAMPLE $> _4.3_stats.$SAMPLE.log
            
            echo "[4.$i.4] quantificando a amostra $SAMPLE com salmon ..."
            salmon quant -1 $SAMPLE.1.fq -2 $SAMPLE.2.fq \
            -o quant_$SAMPLE --libType IU --index idx$tid $> _4.4_quant.$SAMPLE.log

            echo "[4.$i.4] limpando dados de $SAMPLE ..."
            mkdir out_$SAMPLE
            mv qc_$SAMPLE out_$SAMPLE -r
            mv quant_$SAMPLE out_$SAMPLE -r
            rm *.fastq *.fq       
        fi
done 

echo '[5] compactando para RESULTS.zip ...'
zip RESULTS.zip out_*/* *.log

echo $(date) terminado.
