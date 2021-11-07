import Drawable from './Drawable'
import DrawableExon from './DrawableExon'
import DrawableIntron from './DrawableIntron'
import DrawableIsoform from './DrawableIsoform'
import DrawableLocus from './DrawableLocus'
import Histogram from './Histogram'

export default class DrawableGene extends Drawable {
    constructor(
        drawable, gene, bounds
    ) {
        super(drawable, null, bounds, null, gene.nome)
        this.gene = gene
        this.canvas = this;

        if (!gene.fita) {
            this.canvas = new Drawable({ svg: this.flipX() }, null, this.bounds)
        }

        const primary = gene.getPrimaryTranscript();
        const primary_drawable = new Drawable(this.canvas, null, this.bounds, null, 'primary_transcript');

        this.alternative_locus = primary.filter(l => l.tipo === 'Locus').map(i => {
            return new DrawableLocus(primary_drawable, i,
                bounds.scaleX(this.gene.inicio, this.gene.fim, i.inicio, i.fim).withHeight(10))
        })

        this.constitutive_exons = primary.filter(l => l.tipo === 'Exon').map(e => {
            return new DrawableExon(primary_drawable, e,
                bounds.scaleX(this.gene.inicio, this.gene.fim, e.inicio, e.fim).withHeight(10), this.gradient_par('red'))
        })

        this.constitutive_introns = primary.filter(l => l.tipo === 'Intron').map(i => {
            return new DrawableIntron(primary_drawable, i,
                bounds.scaleX(this.gene.inicio, this.gene.fim, i.inicio, i.fim).withHeight(10))
        })

        let y = 0;

        const isoformBound = bounds.withHeight(30);
        this.drawableIsoforms = this.gene.getIsoformas().map(
            (i, idx) => this.drawableIsoform(i, isoformBound.down(y = (25 + 30 * idx))))

        const hadle_loci = (v, d) => {
            return v.map(i => {
                return new DrawableLocus(primary_drawable, i,
                    bounds.scaleX(this.gene.inicio, this.gene.fim, i.inicio, i.fim).withHeight(10).down(d),
                    null, i.draw_strategy)
            })
        }

        this.drawableLociGenomic = hadle_loci(gene.getLoci(['SNP', 'SNV', 'TRANSPOSON']), 80 + y)
        this.drawableLociTranscriptomic = hadle_loci(gene.getLoci(['LNCRNA', 'MIRNA', 'SPLICING']), 150 + y)
        this.drawableLociProteomic = hadle_loci(gene.getLoci(['']))
    }

    drawableIsoform(isoforma, _bounds) {
        const bounds = _bounds.scaleX(this.gene.inicio, this.gene.fim, isoforma.inicio, isoforma.fim)
        return new DrawableIsoform(this.canvas, isoforma, bounds, null)
    }

    draw() {

        this.text(-15, this.bounds.y + 10, this.gene.fita ? '5´' : '3´', { b: true, serif: true })
        this.text(0, this.bounds.y - 10, this.gene.nome, { serif: true })

        this.alternative_locus.forEach(l => l.draw())
        this.constitutive_exons.forEach(l => l.draw())
        this.constitutive_introns.forEach(l => l.draw())

        const amostra = []
        this.gene.getCondicoes().forEach(c => amostra.push(this.gene.getExpressao(c)))
        this.gene.getIsoformas().forEach(i => i.getCondicoes().forEach(c => amostra.push(i.getExpressao(c))))

        let spectro = this.spectro(null, amostra[0], amostra)
        const spct = spectro[0]

        const cd1 = 'Controle'
        const cd2 = 'Tratamento'

        this.text(this.bounds.r + 15, this.bounds.y - 5, cd1, { serif: true, fs: '.8em', r: -45 })
        this.circ(this.bounds.r + 15, this.bounds.y + 5, 8, this.spectro(spct, this.gene.getExpressao(cd1))[1]).style('stroke', 'gray')

        this.text(this.bounds.r + 35, this.bounds.y - 5, cd2, { serif: true, fs: '.8em', r: -45 })
        this.circ(this.bounds.r + 35, this.bounds.y + 5, 8, this.spectro(spct, this.gene.getExpressao(cd2))[1]).style('stroke', 'gray')

        this.drawableIsoforms.forEach(i => {
            i.draw();
            this.circ(this.bounds.r + 15, i.bounds.mv() - 5, 8, this.spectro(spct, i.isoforma.getExpressao(cd1))[1]).style('stroke', 'gray')
            this.circ(this.bounds.r + 35, i.bounds.mv() - 5, 8, this.spectro(spct, i.isoforma.getExpressao(cd2))[1]).style('stroke', 'gray')
            this.text(this.bounds.r + 48, i.bounds.mv(), i.isoforma.nome, { serif: true })
        })

        //Contexto Genomico (Genes, Microsatelite, Transposons, SNV/SNP)
        if (this.drawableLociGenomic.length > 0) {
            this.text(0, this.drawableLociGenomic[0].bounds.y - 20, 'Contexto Genomico', 
            { b: true, s: '#ffa000', c: 'yellow' })
            this.drawableLociGenomic.forEach(l => l.draw())
        }

        //Contexto transcriptomico (ncRNA mi, lnc ...) expressao
        if (this.drawableLociTranscriptomic.length > 0) {
            this.text(0, this.drawableLociTranscriptomic[0].bounds.y - 20, 'Contexto Transcriptomico', 
            { b: true, s: '#ffa000', c: 'yellow' })
            this.drawableLociTranscriptomic.forEach(l => l.draw())
            const hist = new Histogram(this.drawable,
                this.drawableLociTranscriptomic[0].bounds
                    .down(30)
                    .withX(this.bounds.x)
                    .withWidth(this.bounds.width)
                    .withHeight(100)
            )
            hist.draw()
        }

        //Context proteomico (dominios de proteinas)
        const all_doms = []
        this.drawableIsoforms.forEach(i => i.drawableDominios.forEach(d => d.forEach(x => all_doms.push(x))))
        const yd = 15
        if (all_doms.length > 0) {
            this.text(0, 485, 'Contexto Proteomico', { b: true, s: '#ffa000', c: 'yellow' })
            all_doms.sort((a, b) => a.bounds.x - b.bounds.x).forEach((dom, i) => {
                const x = dom.bounds.x, y = 490 + i * yd
                this.rect(x, y, dom.bounds.width, 20, dom.dominio.color, 8).style('stroke', 'gray')
                this.text(x + dom.bounds.width / 2, y + 14, dom.dominio.nome, { hc: true })
            });
        }

    }

}