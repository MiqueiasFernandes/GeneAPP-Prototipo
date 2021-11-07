import Exon from "./Exon";
import Intron from "./Intron";
import Locus from "./Locus";

export default class Gene extends Locus {

    __isoformas = new Array();
    __loci = new Array();

    constructor(
        nome,
        inicio,
        fim,
        fita,
        sequencia,
        nota
    ) {
        super(nome, inicio, fim, fita, 'Gene', sequencia, nota);
    }

    addIsoforma(isoforma) {
        isoforma.gene = this
        this.__isoformas.push(isoforma)
    }

    addLocus(locus) {
        this.__loci.push(locus)
    }

    getLoci(tipos) {
        return this.__loci.filter(l => !tipos || tipos.includes(l.tipo))
    }

    getIsoformas() {
        return this.__isoformas;
    }

    getPrimaryTranscript() {
        // 1. pegar todos pontos

        const exons = this.__isoformas.map(i => i.getExons().map(e => e.inicio + ',' + e.fim).join(',')).join(',');
        const introns = this.__isoformas.map(i => i.getIntrons().map(e => e.inicio + ',' + e.fim).join(',')).join(',');
        const all_pos = [...new Set((exons + ',' + introns).split(','))].map(p => parseInt(p)).sort();

        const const_exons = all_pos.filter(p => this.__isoformas.every(i => i.getExonByPos(p)));
        const const_introns = all_pos.filter(p => this.__isoformas.every(i => i.getIntronByPos(p)));

        const tipo = all_pos.map(p => [p, const_exons.includes(p) ? 'E' : const_introns.includes(p) ? 'I' : 'A']);

        const final = [];
        tipo.forEach(pos => {
            if (final.length < 1) {
                final.push([pos[0], 0, pos[1]])
            } else if (final[final.length - 1][2] === pos[1]) {
                final[final.length - 1][1] = pos[0]
            } else {
                final[final.length - 1][1] = pos[0] - 1
                final.push([pos[0], 0, pos[1]])
            }
        })

        return final.map(p => {
            const inicio = p[0], fim = p[1], tipo = p[2]
            const fita = this.fita;
            if (tipo === 'E') return new Exon(null, inicio, fim, fita)
            if (tipo === 'I') return new Intron(null, inicio, fim, fita)
            return new Locus(null, inicio, fim, fita)
        });
    }
}
