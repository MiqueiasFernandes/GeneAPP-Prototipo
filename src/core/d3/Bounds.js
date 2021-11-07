export default class Bounds {
    constructor(w, h, x=0, y=0, margin = {top: 0, right: 0, bottom: 0, left: 0}) {
        this.x = x;
        this.y = y;
        this.margin = margin;
        this.width = w - margin.left - margin.right;
        this.height = h - margin.top - margin.bottom;
        this.total_with = w;
        this.total_height = h;
        this.r = x + this.width;
    }

    down(pts) {
        return new Bounds(this.width, this.height, this.x, this.y + pts)
    }

    up(pts) {
        return new Bounds(this.width, this.height, this.x, this.y - pts)
    }

    left(pts) {
        return new Bounds(this.width, this.height, this.x + pts, this.y)
    }

    withHeight(h) {
        return new Bounds(this.width, h, this.x, this.y)
    }
    
    withWith(w) {
        return new Bounds(w, this.height, this.x, this.y)
    }

    scaleX(_a1, _a2, _b1, _b2) {
        const a1 = Math.min(_a1, _a2)
        const a2 = Math.max(_a1, _a2)
        const b1 = Math.min(_b1, _b2)
        const b2 = Math.max(_b1, _b2)
        const my_w = 1+a2 - a1
        const o_w = 1+b2 - b1
        const rel = this.width / my_w
        return new Bounds(o_w*rel, this.height, this.x+(b1-a1)*rel, this.y)
    }

    mv() {
        return this.y + this.height / 2;
    }

    mh() {
        return this.x + this.width / 2;
    }

    toString() {
        const margin = this.margin;
        const m = `[^${margin.top} >${margin.right} v${margin.bottom} <${margin.left}]`
        return `x: ${this.x} y: ${this.y} w: ${this.width} h: ${this.height} ${m}`
    }
}