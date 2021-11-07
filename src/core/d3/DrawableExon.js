import Drawable from './Drawable'

export default class DrawableExon extends Drawable {
    constructor(
        drawable, exon, bounds, color='blue'
    ) {
        super(drawable, null, bounds)
        this.exon = exon;
        this.color=color;
    }

    draw() {
        this.rect(this.bounds.x, this.bounds.y, this.bounds.width, this.bounds.height, this.getBasicColor(this.color), 4)
        this.rect(this.bounds.x, this.bounds.y, this.bounds.width, this.bounds.height, this.color, 4)
    }

}