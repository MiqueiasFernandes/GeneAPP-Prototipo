export default class Fasta {
    constructor(data) {
        this.lines = data.split('\n')
        this.id = this.lines[0].split(" ")[0].substr(1);
        this.seq = this.lines.filter(l => !l.startsWith(">")).join("");
    }
}