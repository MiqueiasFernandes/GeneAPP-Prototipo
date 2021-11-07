import Drawable from './Drawable'

export default class DrawableLocus extends Drawable {
    constructor(
        drawable, locus, bounds, fill, draw_strategy
    ) {
        super(drawable, null, bounds)
        this.locus = locus
        if (fill)
            this.rect(bounds.x, bounds.y, this.bounds.width, this.bounds.height, fill)
        this.draw_strategy = draw_strategy
    }

    draw() {
        if (this.draw_strategy) {
            return this.draw_strategy(this)
        }
        const cy = this.bounds.y + this.bounds.height / 2;
        this.line({ x1: this.bounds.x, h: cy, x2: this.bounds.r, sw: 3 })
    }

}