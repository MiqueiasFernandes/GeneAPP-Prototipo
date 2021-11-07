import  Locus  from "./Locus";

export default class Intron extends Locus {

    exon;

    constructor(
        nome,
        inicio,
        fim,
        fita,
        sequencia
    ) {
        super(nome, inicio, fim, fita, 'Intron', sequencia);
    }
}
