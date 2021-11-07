import Drawable from './Drawable'

export default class DrawableIntron extends Drawable {
    constructor(
        drawable, intron, bounds
    ) {
        super(drawable, null, bounds)
    }

    draw() {
        this.wave(this.bounds.x, this.bounds.y, this.bounds.width, 'black', 3);
        this.circ(this.bounds.x, this.bounds.y + this.bounds.height / 2, 3)
        this.circ(this.bounds.r, this.bounds.y + this.bounds.height / 2, 3)
    }

}