import Locus from "./Locus";
import Intron from './Intron';

export default class Isoforma extends Locus {

    exons = [];
    introns = [];
    cds = [];
    cds_sorted = [];
    coords = [];
    dominios = [];

    constructor(
        nome,
        inicio,
        fim,
        fita,
        sequencia
    ) {
        super(nome, inicio, fim, fita, 'Isoforma', sequencia);
    }

    addExon(exon) {
        this.exons.push(exon);
        return this;
    }

    addCDS(cds) {
        this.cds.push(cds);
        this.exons.forEach(e => e.contem(cds.inicio) && (e.cds = cds));
        return this;
    }

    updateIntrons() {
        if (this.exons.length < 1 || this.introns.length > 0) {
            return;
        }
        this.exons = this.exons.sort((e1, e2) => e1.inicio - e2.inicio);
        this.introns = this.exons
            .map(e => `${e.inicio};${e.fim}`)
            .join(',')
            .split(';')
            .slice(1, -1)
            .map(c => c.split(','))
            .map(i => new Intron('intron', parseInt(i[0]) + 1, parseInt(i[1]) - 1, this.fita));
        this.exons.forEach((e, i) => (e.intron = this.introns[i]) && (this.introns[i].exon = this.exons[i + 1]));
    }

    getExons() {
        this.updateIntrons();
        return this.exons;
    }

    getIntrons() {
        this.updateIntrons();
        return this.introns;
    }

    getCDS() {
        if (this.cds.length !== this.cds_sorted.length) {
            this.cds_sorted = this.cds.sort((a, b) => this.fita ? (a.inicio - b.inicio) : (b.fim - a.fim));
        }
        return this.cds_sorted;
    }

    getProteinName() {
        if (this.cds.length < 1)
            return null
        return this.cds[0].nome;
    }

    getExonByPos(pos) {
        for (let i = 0; i < this.exons.length; i++) {
            if (this.exons[i].contem(pos))
                return this.exons[i]
        }
    }

    getIntronByPos(pos) {
        for (let i = 0; i < this.introns.length; i++) {
            if (this.introns[i].contem(pos))
                return this.introns[i]
        }
    }


    addDominio(dominio) {

        // calcular coordernadas dos nucleotideos pela proteina nas duas fitas
        if (this.cds.length > 0 && this.coords < 1) {
            const coords = []
            this.getCDS().forEach((cds, i) => {
                if (coords.length < 1) {
                    coords.push({
                        cds,
                        nt_ini: 1,
                        nt_fim: cds.size,
                        nt_size: cds.size,
                        aa_ini: 1,
                        aa_fim: Math.ceil(cds.size / 3),
                        fase_ini: 0,
                        fase_fim: cds.size % 3
                    })
                } else {
                    const ult = coords[i - 1]
                    const rest = (ult.fase_fim === 0) ? 0 : ((ult.fase_fim === 1) ? 2 : 1);
                    const util = cds.size - rest
                    const ini = ult.aa_fim + (rest > 0 ? 0 : 1)
                    coords.push({
                        cds,
                        nt_ini: ult.nt_fim + 1,
                        nt_fim: ult.nt_fim + cds.size,
                        nt_size: cds.size,
                        aa_ini: ini,
                        aa_fim: ini + Math.ceil(util / 3),
                        fase_ini: rest,
                        fase_fim: util % 3
                    })
                }
            })
            this.coords = coords;
        }

        const dom_ini = dominio.inicio
        const dom_fim = dominio.fim

        dominio.loci = this.coords.map(c => {
            const c_ini = c.aa_ini
            const c_fim = c.aa_fim
            const inicio = c.cds.inicio
            const fim = c.cds.fim
            if (dom_ini <= c_ini && dom_fim >= c_fim) { /// cds dentro do dominio
                return { tipo: 'CDS', ini: inicio, fim, cds: c }
            }
            const dom_ini_dentro_cds = dom_ini >= c_ini && dom_ini <= c_fim
            const dom_fim_dentro_cds = dom_fim >= c_ini && dom_fim <= c_fim
            // distancia entre o inicio e dominio
            const offset_ini = (c.fase_ini > 0 ? c.fase_ini : 3) + (3 * (dom_ini - c_ini - 1));
            let nt_ini = this.fita ? (inicio + offset_ini) : (fim - offset_ini)
            // distancia entre o fim e dominio
            const offset_fim = (c.fase_fim > 0 ? c.fase_fim : 3) + (3 * (c_fim - dom_fim - 1));
            let nt_fim = this.fita ? (fim - offset_fim) : (inicio + offset_fim)

            if (dom_ini_dentro_cds && dom_fim_dentro_cds) { /// dominio dentro de cds
                return { tipo: 'DOMINIO', ini: Math.min(nt_ini, nt_fim), fim: Math.max(nt_ini, nt_fim), cds: c }
            }
            if (dom_ini_dentro_cds) {
                return { tipo: 'INICIO', ini: this.fita ? nt_ini : inicio, fim: this.fita ? fim : nt_ini, cds: c }
            }
            if (dom_fim_dentro_cds) {
                return { tipo: 'FIM', ini: this.fita ? inicio : nt_fim, fim: this.fita ? nt_fim : fim, cds: c }
            }

        }).filter(l => l)
        dominio.isoforma = this;
        this.dominios.push(dominio)
    }

    getDominios = () => this.dominios;
}
