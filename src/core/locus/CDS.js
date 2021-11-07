import Locus from "./Locus";

export default class CDS extends Locus {

    exon;

    constructor(
        nome,
        inicio,
        fim,
        fita,
        sequencia
    ) {
        super(nome, inicio, fim, fita, 'CDS', sequencia);
    }
}
