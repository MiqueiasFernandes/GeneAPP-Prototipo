import File from "../../util/File";

export default class Fasta {
    constructor(data, id) {
        this.lines = data.split('\n')
        this.id = id ? id : this.lines[0].split(" ")[0].substr(1);
        this.seq = this.lines.filter(l => !l.startsWith(">")).join("").replaceAll(/[^\w]/g, '');
        this.size = this.seq.length;
    }

    download() {
        const filename = this.id
        const data = this.lines.join('\n')
        File.download(filename, data)
    }

    static fromMultiFasta(file) {
        return file.split('>').filter(x => x.length > 3).map(f => new Fasta('>' + f))
    }
}