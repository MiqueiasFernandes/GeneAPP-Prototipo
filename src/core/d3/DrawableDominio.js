import Drawable from './Drawable'

export default class DrawableDominio extends Drawable {
    constructor(
        drawable, dominio, bounds, c='#b3b3b3'
    ) {
        super(drawable, null, bounds)
        this.dominio = dominio;
        this.color = this.pattern(c);
    }

    draw() {
        //this.rect(this.bounds.x, this.bounds.y, this.bounds.width, this.bounds.height, this.getBasicColor(this.color), 4)
        this.rect(this.bounds.x, this.bounds.y, this.bounds.width, this.bounds.height, this.color, 4)
    }

}