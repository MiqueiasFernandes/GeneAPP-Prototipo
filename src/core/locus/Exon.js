import Locus from "./Locus";

export default class Exon extends Locus {

    cds;
    intron;

    constructor(
        nome,
        inicio,
        fim,
        fita,
        sequencia
    ) {
        super(nome, inicio, fim, fita, 'Exon', sequencia);
    }

    extrair(gene, seq) {
        const c = this.relativo(gene)
        this.sequence = seq.substr(c[0], c[1]-c[0])
    }
}
