#!/bin/bash

echo "[1] $(date +%D.%H-%M-%S) prepando o ambiente..."

p=1
for prog in sra-toolkit trimmomatic fastqc salmon
    do
    if ! command -v $prog 1> /dev/null 2> /dev/null
    then
        echo "[1.$p] instalando o $prog ..."
        apt install $prog 1> _1.$p\_install.$prog.log 2> _1.$p\_install.$prog.err
        (( p=p+1 ))
    fi
done

echo "usando o salmon ! Versão: $( salmon --version )" 1> _1.0_pacotes.log 2> _1.0_pacotes.log
echo "usando o sra-toolkit ! Versão: $( fastq-dump --version )" >>  _1.0_pacotes.log
echo "usando o trimmomatic ! Versão: $( TrimmomaticPE -version )" >> _1.0_pacotes.log
echo "usando o fastqc ! Versão: $( fastqc --version )" >> _1.0_pacotes.log

for pkg in multiqc
    do
    if ! command --version $pkg 1> /dev/null 2> /dev/null
    then
        echo "[1.$p] instalando o $pkg ..."
        pip install $pkg 1> _1.$p\_install.$pkg.log 2> _1.$p\_install.$pkg.err
        (( p=p+1 ))
    fi
done

echo "usando o multiqc ! $( multiqc --version )" >> _1.0_pacotes.log

echo '[2] importando a anotação ...'
tid=t$(date +%s)
echo '[2.1] baixando o GTF ...'
wget -O gene.$tid.gtf.gz $1 1> _2.1_gtf.download.log 2> _2.1_gtf.download.err
echo '[2.2] descompactando o GTF ...'
gunzip gene.$tid.gtf.gz 1> _2.2_gtf.unzip.log 2> _2.2_gtf.unzip.err

echo '[3] importando transcritos ...'
echo '[3.1] baixando os transcritos ...'
wget -O cds.$tid.fa.gz $2 1> _3.1_transcripts.download.log 2> _3.1_transcripts.download.err
echo '[3.2] descompactando os transcritos ...'
gunzip cds.$tid.fa.gz 1> _3.2_transcripts.unzip.log 2> _3.2_transcripts.unzip.err
echo '[3.3] indexando os transcritos ...'
salmon index -t cds.$tid.fa --index idx$tid 1> _3.3_transcripts.index.log 2> _3.3_transcripts.index.err

echo '[4] quantificando amostras ...'
i=1
for x in $@
    do 
        if [[ `echo $x | grep ,` ]]
        then
            RUN=`echo $x | cut -d, -f1`
            SAMPLE=`echo $x | cut -d, -f2`

            echo "[4.$i.1] obtendo a amostra $SAMPLE pelo acesso $RUN no sra ..."
            fastq-dump --split-3 $RUN 1> _4.$i.1_download.$RUN.$SAMPLE.log 2> _4.$i.1_download.$RUN.$SAMPLE.err
            
            echo "[4.$i.2] fazendo controle de qualidade da amostra $SAMPLE com o TrimmomaticPE ..."
            TrimmomaticPE \
                $RUN\_1.fastq $RUN\_2.fastq \
                $SAMPLE.1.fq $SAMPLE.1.unp.fq \
                $SAMPLE.2.fq $SAMPLE.2.unp.fq \
                ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                1> _4.$i.2_qc.$SAMPLE.log 2> _4.$i.2_qc.$SAMPLE.err
            
            echo "[4.$i.3] reportando controle de qualidade da amostra $SAMPLE com fastqc ..."
            rm qc_$SAMPLE -rf && mkdir qc_$SAMPLE
            fastqc $SAMPLE.1.fq $SAMPLE.2.fq -o qc_$SAMPLE 1> _4.$i.3_stats.$SAMPLE.log 2> _4.$i.3_stats.$SAMPLE.err
            
            echo "[4.$i.4] quantificando com a amostra $SAMPLE com salmon ..."
            salmon quant -1 $SAMPLE.1.fq -2 $SAMPLE.2.fq \
            -o quant_$SAMPLE --libType IU --index idx$tid 1> _4.$i.4_quant.$SAMPLE.log 2> _4.$i.4_quant.$SAMPLE.err

            echo "[4.$i.5] limpando dados de $SAMPLE ..."
            mkdir out_$SAMPLE
            cp quant_$SAMPLE/quant.sf out_$SAMPLE/$SAMPLE.quant.sf
            mv qc_$SAMPLE out_$SAMPLE
            mv quant_$SAMPLE out_$SAMPLE
            rm *.fastq *.fq -f
            (( i=i+1 ))       
        fi
done 

echo '[5] executando o multiqc ...'
multiqc out_*/qc_*

echo '[6] compactando para RESULTS.zip ...'
zip RESULTS.zip out_*/*  multiqc_* *.log *.err

echo $(date +%D.%H-%M-%S) terminado.
