import Drawable from './Drawable'

export default class DrawableCDS extends Drawable {
    constructor(
        drawable, cds, bounds, color='green'
    ) {
        super(drawable, null, bounds)
        this.cds = cds;
        this.color = color;
    }

    draw() {
        this.rect(this.bounds.x, this.bounds.y, this.bounds.width, this.bounds.height, this.getBasicColor(this.color), 4)
        this.rect(this.bounds.x, this.bounds.y, this.bounds.width, this.bounds.height, this.color, 4)
    }

}