import NCBI from "../api/NCBI"
import Fasta from "./Fasta";


const MARGIN = 60

export default class AS {
    constructor(gene) {
        this.gene = gene;
        this.se = null
        this.getSE()
    }


    getSE() {
        if (this.se != null)
            return this.se
        const isoformas = this.gene.getIsoformas();
        const ses = []
        isoformas.forEach(iso1 => {
            const exons = iso1.getExons();
            isoformas.forEach(iso2 => {
                if (iso1 === iso2)
                    return
                const introns = iso2.getIntrons();
                exons
                    .filter(e => e.inicio !== iso1.inicio && e.fim !== iso1.fim) /// se for o 1 exon é inciador alt: LOC101245278
                    .map(e => introns.map(i => [e.inicio >= i.inicio && e.fim <= i.fim, e, i, iso1, iso2]))
                    .forEach(se => se.forEach(ei => ei[0] ? ses.push(ei.slice(1)) : 0))
            })
        });
        return this.se = ses.map(s => new SE(s[0], s[1], s[2], s[3]));
    }
}

// um exon dentro de intron de outra isoforma
class SE {
    constructor(exon, intron, isoform1, isoform2) {
        this.exon = exon
        this.intron = intron
        this.isoform1 = isoform1
        this.isoform2 = isoform2
        this.exon_canonico = isoform1
            .getExons()
            .filter(e => 
                exon.fita ? 
                (e.fim === isoform1.getIntrons().filter(i => i.fim + 1 === exon.inicio)[0].inicio - 1) :
                (e.inicio === isoform1.getIntrons().filter(i => i.inicio -1 === exon.fim)[0].fim +1)
            )[0]
        this.exons = isoform2.getExons()
            .filter(e => e.fim + 1 === intron.inicio || e.inicio - 1 === intron.fim)
            .sort((a, b) => a.fita ? a.inicio - b.inicio : b.inicio - a.inicio)

        this.alt = [isoform1.gene.sequencia, isoform1.fita, [
            ///juncao q inclui o exon alternativo
            [this.exon_canonico.inicio, this.exon_canonico.fim, this.exon_canonico.size],
            [this.exon.inicio, this.exon.fim, this.exon.size]
        ], [
            /// juncao q exclui o exon alternativo
            [this.exons[0].inicio, this.exons[0].fim, this.exons[0].size],
            [this.exons[1].inicio, this.exons[1].fim, this.exons[1].size],
        ]];

        ///verificar se tem 100 up/down stream

        this.canonic_seq = ""
        this.alt_seq = ""
        this.fastas = []

    }

    loadSeqs() {
        /// os exons com < de 100 bp ficam ruin na analise de cov junc
        NCBI.api()
            .get_locus(this.isoform1.gene.sequencia, this.exon_canonico.fim - MARGIN, this.exon_canonico.fim, this.isoform1.fita)
            .then(r => {
                this.canonic_seq = r.seq;
                NCBI.api()
                    .get_locus(this.isoform1.gene.sequencia, this.exon.inicio, this.exon.inicio + MARGIN, this.isoform1.fita)
                    .then(r => this.canonic_seq += r.seq)
                    .then(() => this.fastas.push(new Fasta('>canonic\n' + this.canonic_seq)))
            })
        NCBI.api()
            .get_locus(this.isoform1.gene.sequencia, this.exons[0].fim - MARGIN, this.exons[0].fim, this.isoform1.fita)
            .then(r => {
                this.alt_seq = r.seq;
                NCBI.api()
                    .get_locus(this.isoform1.gene.sequencia, this.exons[1].inicio, this.exons[1].inicio + MARGIN, this.isoform1.fita)
                    .then(r => this.alt_seq += r.seq)
                    .then(() => this.fastas.push(new Fasta('>alternative\n' + this.alt_seq)))
            })

    }
}