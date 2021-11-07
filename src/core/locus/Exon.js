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
}
