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

echo "usando o salmon ! Versão: " > _1.0_pacotes.log && salmon -v 2>> _1.0_pacotes.log
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

echo '[3.4] filtrando os transcritos ...'
cat > script.py << EOF 
seqs = [(l.strip(), l[1:-1].split()) for l in open('cds.fa').readlines() if l.startswith('>')]
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
for l in open('cds.fa').readlines():
  if l.startswith('>'):
    k = l.strip() in ok
  if k:
    as_cds.append(l)
open('cds_filtrada.fna', 'w').writelines(as_cds)
EOF
python3 script.py 1> _3.4_transcripts.filter.log 2> _3.4_transcripts.filter.err
mv cds_filtrada.fna cds.fa
rm script.py

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
multiqc out_*/qc_* 1> _5_multiqc.log 2> _5_multiqc.err

echo '[6] compactando para RESULTS.zip ...'
zip RESULTS.zip out_*/*  multiqc_* *.log *.err 1> _6_zip.log 2> _6_zip.err

echo $(date +%D.%H-%M-%S) terminado.
