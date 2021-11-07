export default class Locus {

    constructor(
        nome = 'unknown',
        inicio,
        fim,
        fita,
        tipo = 'Locus',
        sequencia, nota) {
        this.nome = nome;
        this.inicio = inicio;
        this.fim = fim;
        this.size = 1 + fim - inicio;
        this.fita = fita;
        this.tipo = tipo;
        this.sequencia = sequencia;
        this.nota = nota
        this.expressao = {};
    }

    contem(posicao) {
        return posicao >= this.inicio && posicao <= this.fim;
    }

    naoContem(posicao) {
        return !this.contem(posicao);
    }

    toString() {
        return `${this.tipo}: ${this.nome} ${this.fita ? '+' : '-'}[${this.inicio}:${this.fim}]`
    }

    getSequence() {
        return this.sequencia;
    }

    setExpressao(tpm, cond) {
        this.expressao[cond] = tpm
        return this
    }
    getExpressao(cond) {
        return this.expressao[cond]
    }
    getCondicoes() {
        return Object.keys(this.expressao)
    }
    set_draw_strategy(draw_strategy) {
        this.draw_strategy = draw_strategy;
        return this;
    }

}