import Drawable from "./Drawable";
import * as d3 from "d3";

export default class Violin extends Drawable {

    constructor(
        drawable, bounds, data = {
            tipo1: [1, 2, 2, 3, 3, 3, 4, 4, 5],
            tipo2: [1, 2, 2, 2, 3, 3, 4, 5, 5],
            tipo3: [1, 1, 2, 3, 3, 4, 4, 4, 5]
        }
    ) {
        super(drawable, null, bounds, null, 'violin')
        this.data = data;
    }

    plot() {

        const all = Object.values(this.data)
            .map(x => x.join(';')).join(';').split(';').map(x => parseFloat(x))

        const tipos = Object.keys(this.data);

        var y = d3.scaleLinear()
            .domain([d3.min(all), d3.max(all)])
            .range([this.bounds.height, 0])
        this.svg.append("g").call(d3.axisLeft(y))

        // Build and Show the X scale. It is a band scale like for a boxplot: each group has an dedicated RANGE on the axis. This range has a length of x.bandwidth
        var x = d3.scaleBand()
            .range([0, this.bounds.width])
            .domain(tipos)
            .padding(0.05)     // This is important: it is the space between 2 groups. 0 means no padding. 1 is the maximum.
        this.svg.append("g")
            .attr("transform", "translate(0," + this.bounds.height + ")")
            .call(d3.axisBottom(x))

        // // Features of the histogram
        var histogram = d3.bin()
            .domain(y.domain())
            //.thresholds(y.ticks(3))  
            .value(d => d.value)

        var sumstat = tipos.map(t => ({ key: t, value: histogram(this.data[t].map(x => ({ key: t, value: x }))) }))

        console.log(sumstat)

        var maxNum = d3.max(Object.values(sumstat.map(x => x.value)).map(x => d3.max(x.map(z => z.length))));
        console.log(maxNum)

        // // The maximum width of a violin must be x.bandwidth = the width dedicated to a group
        var xNum = d3.scaleLinear()
            .range([0, x.bandwidth()])
            .domain([-maxNum, maxNum])

        // // Add the shape to this svg!
        this.svg
            .selectAll("myViolin")
            .data(sumstat)
            .enter()
            .append("g")
            .attr("transform", (d) => `translate(${x(d.key)},0)`)
            .append("path")
            .datum(d => d.value)
            .style("stroke", "none")
            .style("fill", "#1f77b4")
            .attr("d", d3.area()
                .x0((d) => xNum(d.length))
                .x1((d) => xNum(-d.length))
                .y((d) => y(d.x0))
                .curve(d3.curveCatmullRom)
            )

    }
}