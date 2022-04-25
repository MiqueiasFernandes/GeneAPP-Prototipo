#!/bin/bash

echo "[  1  ] $(date +%D.%H-%M-%S) prepando o ambiente..."
tid=t$(date +%s)
mkdir results$tid && cd results$tid

p=1

## sra-toolkit : https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
## trimmomatic : http://www.usadellab.org/cms/?page=trimmomatic   http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
## fastqc      : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## salmon      : https://salmon.readthedocs.io/en/latest/salmon.html
## bamtools    : https://www.htslib.org/doc/samtools.html
## rna-star    : https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

for prog in sra-toolkit trimmomatic fastqc salmon bamtools rna-star
    do
    if ! command -v $prog 1> /dev/null 2> /dev/null
    then
        echo "[1.$p  ] instalando o $prog ..."
        apt install $prog -y 1> _1.$p\_install.$prog.log 2> _1.$p\_install.$prog.err
        (( p=p+1 ))
    fi
done

echo "usando o salmon ! Versão: " > _1.0_pacotes.log && salmon -v 2>> _1.0_pacotes.log
echo "usando o star ! Versão: $( STAR --version )" >>  _1.0_pacotes.log
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
        echo "usando o $pkg ! $( pip list | grep $pkg )" >> _1.0_pacotes.log
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

echo '[  2  ] $(date +%D.%H-%M-%S) importando genoma e a anotação ...'
echo '[2.1  ] baixando o genoma ...'
wget -O genoma.$tid.fa.gz $1 1> _2.1_genoma.download.log 2> _2.1_genoma.download.err
echo '[2.2  ] descompactando o genoma ...'
gunzip genoma.$tid.fa.gz 1> _2.2_genoma.unzip.log 2> _2.2_genoma.unzip.err
echo '[2.3  ] baixando o GTF ...'
wget -O gene.$tid.gtf.gz $2 1> _2.3_gtf.download.log 2> _2.3_gtf.download.err
echo '[2.4  ] descompactando o GTF ...'
gunzip gene.$tid.gtf.gz 1> _2.4_gtf.unzip.log 2> _2.4_gtf.unzip.err

echo '[  3  ] $(date +%D.%H-%M-%S) importando transcritos ...'
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
lns = [l.strip().split('\t') for l in open(gtf).readlines()]
genes = [[x[0], x[6] == '+', int(x[3]), int(x[4]), x[-1].split('gene_id ')[1].split(';')[0].replace('"', '').strip()] 
 for x in lns if len(x) > 2 and x[2] == 'gene']
gen_acecc = set([l.split('gene=')[1].split()[0].replace(']', '') for l in open(cds).readlines() if l.startswith('>')])
gns = [x for x in genes if x[-1] in gen_acecc]
print(len(gen_acecc) - len(gns), 'genes nao encontrados na busca por ID')
seqs = SeqIO.to_dict(SeqIO.parse(genoma, "fasta"))
gene_seqs = []
for seq, strand, ini, end, name in gns:
  s = Seq.Seq(str(seqs[seq].seq[ini-1:end]))
  if not strand:
    s = s.reverse_complement()
  s = SeqRecord.SeqRecord(s, id=name, description="")
  gene_seqs.append(s)
SeqIO.write(gene_seqs, 'gene_seqs.fa', 'fasta')
print('finalizado.')
EOF
python3 script.py 1> _3.4_genes.extract.log 2> _3.4_genes.extract.err
rm script.py
mkdir idxgenes
echo '[3.5  ] indexando sequencia de genes ...'
STAR --runMode genomeGenerate --genomeDir idxgenes  --genomeFastaFiles gene_seqs.fa 1> _3.5_genes.index.log 2> _3.5_genes.index.err

echo '[3.6  ] indexando os transcritos ...'
salmon index -t cds.$tid.fa --index idx$tid 1> _3.6_transcripts.index.log 2> _3.6_transcripts.index.err

echo '[  4  ] quantificando amostras ...'
i=1
for x in $@
    do 
        if [[ `echo $x | grep ,` ]]
        then
            RUN=`echo $x | cut -d, -f1`
            SAMPLE=`echo $x | cut -d, -f2`

            echo "[4.$i.1] $(date +%D.%H-%M-%S) obtendo a amostra $SAMPLE pelo acesso $RUN no sra ..."
            fastq-dump --split-3 --minReadLen 80 $RUN 1> _4.$i.1_download.$RUN.$SAMPLE.log 2> _4.$i.1_download.$RUN.$SAMPLE.err
            
            echo "[4.$i.2] fazendo controle de qualidade da amostra $SAMPLE com o trimmomatic ..."
            TrimmomaticPE \
                $RUN\_1.fastq $RUN\_2.fastq \
                $SAMPLE.F.fq $SAMPLE.1.unp.fq \
                $SAMPLE.R.fq $SAMPLE.2.unp.fq \
                ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                1> _4.$i.2_qc.$SAMPLE.log 2> _4.$i.2_qc.$SAMPLE.err
            
            echo "[4.$i.3] reportando controle de qualidade da amostra $SAMPLE com fastqc ..."
            rm qc_$SAMPLE -rf && mkdir qc_$SAMPLE
            fastqc $SAMPLE.F.fq $SAMPLE.R.fq -o qc_$SAMPLE 1> _4.$i.3_stats.$SAMPLE.log 2> _4.$i.3_stats.$SAMPLE.err
            
            echo "[4.$i.4] quantificando com a amostra $SAMPLE com salmon ..."
            salmon quant -1 $SAMPLE.F.fq -2 $SAMPLE.R.fq \
            -o quant_$SAMPLE --libType IU --index idx$tid 1> _4.$i.4_quant.$SAMPLE.log 2> _4.$i.4_quant.$SAMPLE.err

            echo "[4.$i.5] mapeando com a amostra $SAMPLE com star ..."
            STAR --genomeDir idxgenes \
            --readFilesIn $SAMPLE.F.fq $SAMPLE.R.fq \
            --outSAMtype BAM SortedByCoordinate 1> _4.$i.5_map.$SAMPLE.log 2> _4.$i.5_map.$SAMPLE.err

            echo "[4.$i.6] gerando arquivo de cobertura para a amostra $SAMPLE com deeptools ..."
            bamCoverage --bam Aligned.sortedByCoord.out.bam -o $SAMPLE.bw --binSize 3 1> _4.$i.6_cov.$SAMPLE.log 2> _4.$i.6_cov.$SAMPLE.err

            echo "[4.$i.7] limpando dados de $SAMPLE ..."
            mkdir out_$SAMPLE
            cp quant_$SAMPLE/quant.sf out_$SAMPLE/$SAMPLE.quant.sf
            mv qc_$SAMPLE out_$SAMPLE
            mv quant_$SAMPLE out_$SAMPLE
            mv Aligned.sortedByCoord.out.bam out_$SAMPLE/$SAMPLE.sortedByCoord.out.bam
            mv $SAMPLE.bw out_$SAMPLE
            rm *.fastq *.fq *.bam* -f
            (( i=i+1 ))       
        fi
done 

echo '[  5  ] $(date +%D.%H-%M-%S) executando o multiqc ...'
multiqc out_*/qc_* 1> _5_multiqc.log 2> _5_multiqc.err

echo '[  6  ] $(date +%D.%H-%M-%S) compactando para RESULTS.zip ...'
zip RESULTS.zip out_*/*  multiqc_* *.log *.err 1> _6_zip.log 2> _6_zip.err
cp RESULTS.zip ../ && cd ..

echo $(date +%D.%H-%M-%S) terminado.
