import Drawable from './Drawable'
import DrawableIntron from './DrawableIntron'
import DrawableExon from './DrawableExon'
import DrawableCDS from './DrawableExon'
import DrawableDominio from './DrawableDominio'

export default class DrawableIsoform extends Drawable {
    constructor(
        drawable, isoforma, bounds
    ) {
        super(drawable, null, bounds, null, isoforma.nome)
        this.isoforma = isoforma;
        const _bounds = bounds.down(5).withHeight(bounds.height - 20);
        this.drawableIntrons = isoforma.getIntrons().map((i) => this.drawableIntron(i, _bounds))
        this.drawableExons = isoforma.getExons().map((e) => this.drawableExon(e, _bounds))
        this.drawableCDSs = isoforma.getCDS().map((c) => this.drawableCDS(c, _bounds.withHeight(bounds.height - 16).up(2)))
        this.drawableDominios = isoforma.getDominios().map((d) => this.drawableDominio(d, _bounds.withHeight(bounds.height - 16).up(2)))
    }

    drawableIntron(intron, _bounds) {
        const bounds = _bounds.scaleX(this.isoforma.inicio, this.isoforma.fim, intron.inicio, intron.fim)
        return new DrawableIntron(this, intron, bounds, null)
    }

    drawableExon(exon, _bounds) {
        const bounds = _bounds.scaleX(this.isoforma.inicio, this.isoforma.fim, exon.inicio, exon.fim)
        return new DrawableExon(this, exon, bounds, this.gradient_par('green'))
    }

    drawableCDS(cds, _bounds) {
        const bounds = _bounds.scaleX(this.isoforma.inicio, this.isoforma.fim, cds.inicio, cds.fim)
        return new DrawableCDS(this, cds, bounds, this.gradient_par('blue'))
    }

    drawableDominio(dominio, _bounds) {
        return dominio.loci.map(l => new DrawableDominio(
            this, dominio,
            _bounds.scaleX(this.isoforma.inicio, this.isoforma.fim, l.ini, l.fim)
        ))
    }

    draw() {
        this.drawableIntrons.forEach(i => i.draw())
        this.drawableExons.forEach(e => e.draw())
        this.drawableCDSs.forEach(c => c.draw())
        this.drawableDominios.forEach(d => d.forEach(l => l.draw()))
    }

}