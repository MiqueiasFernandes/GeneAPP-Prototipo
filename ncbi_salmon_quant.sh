#!/bin/bash

MIN_READ_LEN=80

tid=t$(date +%s)
mkdir results$tid && cd results$tid
echo "[1    ] $(date +%D.%H-%M-%S) prepando o ambiente results$tid/ ..."

p=1

## sra-toolkit : https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
## trimmomatic : http://www.usadellab.org/cms/?page=trimmomatic   http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
## fastqc      : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## salmon      : https://salmon.readthedocs.io/en/latest/salmon.html
## samtools    : https://www.htslib.org/doc/samtools.html
## hisat2      : http://daehwankimlab.github.io/hisat2/manual/

for prog in sra-toolkit trimmomatic fastqc salmon samtools bamtools hisat2
    do
    if ! command -v $prog 1> /dev/null 2> /dev/null
    then
        echo "[1.$p  ] instalando o $prog ..."
        apt install $prog -y 1> _1.$p\_install.$prog.log 2> _1.$p\_install.$prog.err
        (( p=p+1 ))
    fi
done

echo "usando o salmon ! Versão: " > _1.0_pacotes.log && salmon -v 2>> _1.0_pacotes.log
echo "usando o hisat2 ! Versão: $( hisat2 --version )" >>  _1.0_pacotes.log
echo "usando o bamtools ! Versão: $( bamtools --version )" >>  _1.0_pacotes.log
echo "usando o sra-toolkit ! Versão: $( fastq-dump --version )" >>  _1.0_pacotes.log
echo "usando o trimmomatic ! Versão: $( TrimmomaticPE -version )" >> _1.0_pacotes.log
echo "usando o fastqc ! Versão: $( fastqc --version )" >> _1.0_pacotes.log

## multiqc   : https://multiqc.info/docs/
## biopython : https://biopython.org/
## deeptools : https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

for pkg in multiqc biopython deeptools
    do
    if [[ ! `pip list | grep $pkg` ]]
        then 
        echo "[1.$p  ] instalando o $pkg ..."
        pip install $pkg 1> _1.$p\_install.$pkg.log 2> _1.$p\_install.$pkg.err
        (( p=p+1 ))
    fi
    if [[ `pip list | grep $pkg` ]]
        then 
        echo "usando o $pkg ! $( pip list | grep $pkg | tr -s \  \   )" >> _1.0_pacotes.log
    else
        echo ERRO: ao instalar pacote $pkg
    fi
done

if ! grep 'pysam.index(bamFile)' /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py
    then
    cp /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py .
    cp bamHandler.py /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py.old
    grep -B100000 'bam = pysam.Samfile(bamFile' bamHandler.py | grep -v 'bam = pysam.Samfile(bamFile' > xtemp
    echo '        pysam.index(bamFile)' >> xtemp
    tail bamHandler.py -n+`grep -n 'bam = pysam.Samfile(bamFile' bamHandler.py | cut -d: -f1` >> xtemp
    cp xtemp /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py
    rm xtemp bamHandler.py
fi

echo "[2    ] $(date +%D.%H-%M-%S) importando genoma e a anotação ..."
echo '[2.1  ] baixando o genoma ...'
wget -O genoma.$tid.fa.gz $1 1> _2.1_genoma.download.log 2> _2.1_genoma.download.err
echo '[2.2  ] descompactando o genoma ...'
gunzip genoma.$tid.fa.gz 1> _2.2_genoma.unzip.log 2> _2.2_genoma.unzip.err
echo '[2.3  ] baixando o GTF ...'
wget -O gene.$tid.gtf.gz $2 1> _2.3_gtf.download.log 2> _2.3_gtf.download.err
echo '[2.4  ] descompactando o GTF ...'
gunzip gene.$tid.gtf.gz 1> _2.4_gtf.unzip.log 2> _2.4_gtf.unzip.err

echo "[3    ] $(date +%D.%H-%M-%S) importando transcritos ..."
echo '[3.1  ] baixando os transcritos ...'
wget -O cds.$tid.fa.gz $3 1> _3.1_transcripts.download.log 2> _3.1_transcripts.download.err
echo '[3.2  ] descompactando os transcritos ...'
gunzip cds.$tid.fa.gz 1> _3.2_transcripts.unzip.log 2> _3.2_transcripts.unzip.err

echo '[3.3  ] filtrando os transcritos ...'
echo "cds = 'cds.$tid.fa'" > script.py
cat >> script.py << EOF 
seqs = [(l.strip(), l[1:-1].split()) for l in open(cds).readlines() if l.startswith('>')]
print(len(seqs), 'sequencias de CDS')
pars = [[a, b[0], c] for a, b, c in [[x[0], [z for z in x if 'gene=' in z], y] for y, x in seqs] if len(b) == 1]
conts = {g: [0, []] for g in set([x[1] for x in pars])}
print(len(conts), 'genes')
for a, b, c in pars:
  conts[b][0]+=1
  conts[b][1].append(c)
print(len(set([k for k, v in conts.items() if v[0] < 2])), 'genes sem AS')
print(len(set([k for k, v in conts.items() if v[0] > 1])), 'genes com AS')
ok = []
for k, v in conts.items():
  if v[0] > 1:
    for seq in v[1]:
      ok.append(seq)
print(len(ok), 'CDS de genes com AS')
k=False
as_cds = []
for l in open(cds).readlines():
  if l.startswith('>'):
    k = l.strip() in ok
  if k:
    as_cds.append(l)
open(cds, 'w').writelines(as_cds)
EOF
python3 script.py 1> _3.3_transcripts.filter.log 2> _3.3_transcripts.filter.err
rm script.py

echo '[3.4  ] extraindo sequencia de genes ...'
echo "cds = 'cds.$tid.fa'" > script.py
echo "genoma = 'genoma.$tid.fa'" >> script.py
echo "gtf = 'gene.$tid.gtf'" >> script.py
cat >> script.py << EOF
from Bio import SeqIO, Seq, SeqRecord
gen_acecc = set([l.split('gene=')[1].split()[0].replace(']', '') for l in open(cds).readlines() if l.startswith('>')])
gns = [l.strip().split('\t') for l in open(gtf).readlines() if '\tgene\t' in l]
cords = [[x[0], int(x[3]), int(x[4]), x[6] == '+', x[-1].split('"')[1]] for x in gns]
print(len(cords), 'genes no GTF')
cords = [x for x in cords if x[-1] in gen_acecc]
print(len(cords), 'genes filtrados')
seqs = SeqIO.to_dict(SeqIO.parse(genoma, 'fasta'))
gseqsF = [SeqRecord.SeqRecord(seqs[s[0]].seq[s[1]-1:s[2]], id=s[-1], description='') for s in cords if s[3]]
gseqsR = [SeqRecord.SeqRecord(seqs[s[0]].seq[s[1]:s[2]+1].reverse_complement(), id=s[-1], description='') for s in cords if not s[3]]
SeqIO.write(gseqsF+gseqsR, 'gene_seqs.fa', 'fasta')
print('finalizado.')
EOF
python3 script.py 1> _3.4_genes.extract.log 2> _3.4_genes.extract.err
rm script.py
echo '[3.5  ] indexando sequencia de genes ...'
hisat2-build gene_seqs.fa idxgenes 1> _3.5_genes.index.log 2> _3.5_genes.index.err

echo '[3.6  ] indexando os transcritos ...'
salmon index -t cds.$tid.fa --index idx$tid 1> _3.6_transcripts.index.log 2> _3.6_transcripts.index.err

echo "[4    ] $(date +%D.%H-%M-%S) quantificando amostras ..."
i=1
for x in $@
    do 
        if [[ `echo $x | grep ,` ]]
        then
            RUN=`echo $x | cut -d, -f1`
            SAMPLE=`echo $x | cut -d, -f2`

            echo "[4.$i.1] $(date +%D.%H-%M-%S) obtendo a amostra $SAMPLE pelo acesso $RUN no sra ..."
            fastq-dump --split-3 --minReadLen $MIN_READ_LEN $RUN 1> _4.$i.1_download.$RUN.$SAMPLE.log 2> _4.$i.1_download.$RUN.$SAMPLE.err
            
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                echo "[4.$i.2] fazendo controle de qualidade da amostra $SAMPLE com o TrimmomaticPE ..."
                TrimmomaticPE \
                    $RUN\_1.fastq $RUN\_2.fastq \
                    $SAMPLE.F.fq $SAMPLE.1.unp.fq \
                    $SAMPLE.R.fq $SAMPLE.2.unp.fq \
                    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                    1> _4.$i.2_qc.$SAMPLE.log 2> _4.$i.2_qc.$SAMPLE.err
                else
                echo "[4.$i.2] fazendo controle de qualidade da amostra $SAMPLE com o TrimmomaticSE ..."
                TrimmomaticSE \
                    $RUN.fastq $SAMPLE.fq \
                    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                    1> _4.$i.2_qc.$SAMPLE.log 2> _4.$i.2_qc.$SAMPLE.err
            fi
            
            echo "[4.$i.3] reportando controle de qualidade da amostra $SAMPLE com fastqc ..."
            rm qc_$SAMPLE -rf && mkdir qc_$SAMPLE
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    fastqc $SAMPLE.F.fq $SAMPLE.R.fq -o qc_$SAMPLE 1> _4.$i.3_stats.$SAMPLE.log 2> _4.$i.3_stats.$SAMPLE.err
                else
                    fastqc $SAMPLE.fq -o qc_$SAMPLE 1> _4.$i.3_stats.$SAMPLE.log 2> _4.$i.3_stats.$SAMPLE.err
            fi

            echo "[4.$i.4] quantificando com a amostra $SAMPLE com salmon ..."
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    salmon quant -1 $SAMPLE.F.fq -2 $SAMPLE.R.fq \
            -o quant_$SAMPLE --libType IU --index idx$tid 1> _4.$i.4_quant.$SAMPLE.log 2> _4.$i.4_quant.$SAMPLE.err
                else
                    salmon quant -r $SAMPLE.fq \
            -o quant_$SAMPLE --libType IU --index idx$tid 1> _4.$i.4_quant.$SAMPLE.log 2> _4.$i.4_quant.$SAMPLE.err
            fi

            echo "[4.$i.5] mapeando com a amostra $SAMPLE com hisat2 ..."
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    hisat2 -x idxgenes -1 $SAMPLE.F.fq -2 $SAMPLE.R.fq --no-unal -S $SAMPLE.maped.sam  1> _4.$i.5_map.$SAMPLE.log 2> _4.$i.5_map.$SAMPLE.err
                else
                    hisat2 -x idxgenes -U $SAMPLE.fq --no-unal -S $SAMPLE.maped.sam  1> _4.$i.5_map.$SAMPLE.log 2> _4.$i.5_map.$SAMPLE.err
            fi
            
            echo "[4.$i.6] transformando sam em bam ordenado do $SAMPLE com samtools ..."
            samtools view -S -b $SAMPLE.maped.sam > $SAMPLE.maped.bam  2> _4.$i.6_bam.$SAMPLE.err
            bamtools sort -in $SAMPLE.maped.bam -out $SAMPLE.sorted.bam  1> _4.$i.6_bam.$SAMPLE.log 2>> _4.$i.6_bam.$SAMPLE.err

            echo "[4.$i.7] gerando arquivo de cobertura para a amostra $SAMPLE com deeptools ..."
            GENE=$(grep \> gene_seqs.fa | head -1000 | tail -1 | tr -d \> | cut -d\   -f1)
            bamCoverage -b $SAMPLE.sorted.bam -o $SAMPLE.bed --outFileFormat bedgraph --binSize 3 -p 2 -r $GENE 1> _4.$i.7_cov.$SAMPLE.log 2> _4.$i.7_cov.$SAMPLE.err

            echo "[4.$i.8] limpando dados de $SAMPLE ..."
            mkdir out_$SAMPLE
            cp quant_$SAMPLE/quant.sf out_$SAMPLE/$SAMPLE.quant.sf
            mv qc_$SAMPLE out_$SAMPLE
            mv quant_$SAMPLE out_$SAMPLE
            mv $SAMPLE.bed out_$SAMPLE
            mv $SAMPLE.sorted.bam out_$SAMPLE
            rm *.fastq *.fq *.bam* -f
            (( i=i+1 ))       
        fi
done 

echo "[5    ] $(date +%D.%H-%M-%S) executando o multiqc ..."
multiqc out_*/qc_* 1> _5_multiqc.log 2> _5_multiqc.err

echo "[6    ] $(date +%D.%H-%M-%S) compactando para RESULTS.zip ..."
zip RESULTS.zip out_*/*  multiqc_* *.log *.err 1> _6_zip.log 2> _6_zip.err
cp RESULTS.zip ../ && cd ..

echo $(date +%D.%H-%M-%S) terminado.
