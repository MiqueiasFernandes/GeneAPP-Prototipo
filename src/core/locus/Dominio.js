export default class Dominio {

    constructor(
        nome,
        inicio,
        fim,
        color='#ffa'
    ) {

        this.nome = nome
        this.inicio = inicio
        this.fim = fim
        this.isoforma = undefined
        this.loci = []
        this.color = color
    }
}
