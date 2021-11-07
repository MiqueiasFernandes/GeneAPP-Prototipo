import * as d3 from "d3";
import Drawable from './Drawable'

export default class Histogram extends Drawable {

    constructor(
        drawable, bounds, data = {
            samp1: { rep1: { 1: 0, 2: 5, 3: 2, 6: 1, 7:0}, rep2: { 2: 0,3:2, 4:4,  5: 1.5,6:2, 7:0 }, rep3: { 6:0,7:5,8:0 } },
            samp2: { rep1: { 2: 0, 3: 2, 6: 5, 7:0}, rep2: { 1: 0,3:2, 4:4,  5: 1.5,6:2, 7:0 }, rep3: { 5:0,6:5,7:0 } }
        }
    ) {
        super(drawable, null, bounds, null, 'histograma')
        this.data = data;
    }

    getAllData(dt, fn = (x) => Object.values(x), conv = x => x) {
        let all_data = [];
        Object.values(dt).forEach( // samples
            s => Object.values(s).forEach( // replicates
                r => fn(r).forEach( // positions || valores
                    x => all_data = all_data.concat(x)
                )
            )
        );
        return all_data.map(x => conv(x))
    }

    eixo(amostra, modo, y, range, config = x => x, label) {
        const eixo = d3.scaleLinear()
            .domain([Math.min(...amostra), Math.max(...amostra)])
            .range(range);
        this.svg.append("g")
            .attr("transform", `translate(0,${y})`)
            .call(config(modo(eixo)));
        if (label) {
            const h = Math.abs(range[0]-range[1])
            if (label.x) this.text(this.bounds.mh(), y+30, label.x, {hc: true, b: true})
            if (label.y) this.text(this.bounds.x - 20, y + h/2, label.y, {r: -90, b: true, hc: true})
        }
        return eixo;
    }

    grafico(X, Y, y, data, c = "#69b3a2") { // padrao de dados [x, y -> altura]
        this.svg.append("path")
            .attr("transform", `translate(${this.bounds.x},${y})`)
            .datum(data)
            .attr("fill", c)
            .attr("opacity", ".6")
            .attr("stroke", "#000")
            .attr("stroke-width", 1)
            .attr("stroke-linejoin", "round")
            .attr("d", d3.line()
                .curve(d3.curveBasis)
                .x(d => X(d[0]))
                .y(d => Y(d[1]))
            );
    }

    draw() {
        const x_dados = [... new Set(this.getAllData(this.data, (x) => Object.keys(x), parseInt))]
        const s1_dados = [... new Set(this.getAllData({ sp: this.data.samp1 }))]
        const s2_dados = [... new Set(this.getAllData({ sp: this.data.samp2 }))]
        const h = this.bounds.height / 2
        const X = this.eixo(x_dados, d3.axisBottom, this.bounds.d, [0, this.bounds.width], x => x)
        const Y1 = this.eixo(s1_dados, d3.axisLeft, this.bounds.y, [h, 0], (x) => x.ticks(3), {y:'CTRL'})
        const Y2 = this.eixo(s2_dados, d3.axisLeft, this.bounds.mv(), [0, h], (x) => x.ticks(3), {y:'TRT'})

        Object.values(this.data.samp1).forEach((r) =>{
            const data = Object.entries(r).map(e => [parseInt(e[0]),e[1]])
            this.grafico(X, Y1, this.bounds.y, data, "#69b3a2")
        })

        Object.values(this.data.samp2).forEach((r) =>{
            const data = Object.entries(r).map(e => [parseInt(e[0]),e[1]])
            this.grafico(X, Y2, this.bounds.mv(), data, "#404080")
        })

    }
}
