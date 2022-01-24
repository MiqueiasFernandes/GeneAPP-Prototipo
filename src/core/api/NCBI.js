import axios from "axios";
import Isoforma from "../locus/Isoforma";
import Exon from "../locus/Exon";
import CDS from "../locus/CDS";
import Gene from "../locus/Gene";
import Fasta from "../locus/Fasta";
import LocalStorage from "../../util/LocalStorage";


export default class NCBI {

    static HOST = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
    static API_KEY = "50d6c4f1c9e2808a7f22b6023cdeccfa5809";
    static INSTANCE = null;
    fila = []

    config() {
        /// configurar concorrencia
        // https://github.com/axios/axios/issues/994
        // http://bluebirdjs.com/docs/api/promise.map.html
        // axios.interceptors.request.use((config) => {
        //     if (config.url === NCBI.HOST) {
        //         const prom = new Promise((rs, rj) => {
        //             rs(config)
        //         });
        //         this.fila.push([config.params, prom])
        //         return prom;
        //     }
        //     return config;
        // });

        // const dequeue = (config) => {
        //     if (config.url === NCBI.HOST) {
        //         this.fila.pop();
        //     }
        //     return config;
        // }

        // axios.interceptors.response.use(dequeue, dequeue);
    }

    static api() {
        if (!NCBI.INSTANCE) {
            NCBI.INSTANCE = new NCBI();
            NCBI.INSTANCE.config();
        }
        NCBI.HOST = LocalStorage.instance().get("ncbi", NCBI.HOST);
        NCBI.API_KEY = LocalStorage.instance().get("API_KEY", NCBI.API_KEY);
        return NCBI.INSTANCE;
    }



    get_sequence(id) {
        return axios
            .get(NCBI.HOST, {
                params: {
                    api_key: NCBI.API_KEY,
                    db: 'protein',
                    rettype: 'fasta',
                    id
                }
            })
    }

    get_locus(id, ini, fim, strand = true, label) {
        return axios
            .get(NCBI.HOST, {
                params: {
                    api_key: NCBI.API_KEY,
                    db: 'nuccore',
                    rettype: 'fasta',
                    retmode: 'text',
                    strand: strand ? 1 : 2,
                    seq_start: ini,
                    seq_stop: fim,
                    id
                }
            })
            .then(r => new Fasta(r.data, label))

    }

    async_load(genes) {
        return genes.map((g) => new GeneTable(g));
    }

    getGenes(genes) {
        return genes.map((g) => this.parseGene(g));
    }

    parseGene(gene_table) {
        const { name, parsed } = gene_table;
        const { organism, seq, gid, from, to, strand, mrnas, gene } = parsed;
        const _gene = new Gene(
            name,
            from,
            to,
            strand,
            seq,
            `${organism} ${gene} ${gid}`
        );

        _gene.organism = organism;

        const isos = mrnas.map((mrna) => {
            const exons = mrna.exons.map(
                (e, idx) => new Exon("e" + idx, e[0], e[1], strand)
            );

            const isoforma = new Isoforma(
                mrna.mrna_id,
                Math.min(...exons.map((e) => e.inicio)),
                Math.max(...exons.map((e) => e.fim)),
                strand
            );
            exons.forEach((e) => isoforma.addExon(e));
            mrna.cds.forEach((c) =>
                isoforma.addCDS(new CDS(mrna.protein, c[0], c[1], strand))
            );
            return isoforma;
        });
        isos.forEach((iso) => _gene.addIsoforma(iso));
        return _gene;
    }

}

class GeneTable {
    data = null;
    parsed = null;
    status = "unload";
    constructor(name) {
        this.name = name;
        axios
            .get(NCBI.HOST, {
                params: {
                    api_key: NCBI.API_KEY,
                    db: 'gene',
                    rettype: 'gene_table',
                    retmode: 'text',
                    id: name
                }
            })
            .then((resp) => {
                this.status = "loaded";
                this.data = resp.data;
                this.getParsed();
            })
            .catch((e) => {
                console.error(e)
                this.status = "error"
            });
    }

    parse() {
        const parsed = this.data.split("\n");
        const gene = parsed[0].split("[")[0];
        const organism = parsed[0].split("[")[1].split("]")[0];
        const gid = parsed[1].split("Gene ID:")[1].split(",")[0].trim();
        const head = parsed.filter(l => l.startsWith("Reference ") && l.includes(" Assembly ") && l.includes("from: ") && l.includes(" to:"))[0]
        const ref = head.split("Primary Assembly")[1];
        const seq = ref.trimStart().split(" ")[0].trim();
        const strand = !ref.includes("(minus strand)");
        const a = parseInt(ref.split("from:")[1].split("to")[0].trim());
        const b = parseInt(ref.split("to: ")[1].split(" ")[0]);
        const from = a < b ? a : b;
        const to = a > b ? a : b;
        let mrnas = [];
        parsed
            .slice(parsed.indexOf(head) + 1)
            .join("\n")
            .split("\n\n")
            .forEach(
                (e, i) => (i % 2 < 1 && mrnas.push([e])) || mrnas[(i - 1) / 2].push(e)
            );

        mrnas = mrnas
            .filter((m) => m.length > 1)
            .map((m) => m[0].split("\n").concat(m[1].split("\n")))
            .filter(
                (m) =>
                    m[1].startsWith("protein") &&
                    m[2].startsWith("Exon table for ") &&
                    m[3].startsWith("Genomic Interval Exon") &&
                    m[4].startsWith("---------------------------")
            );

        mrnas = mrnas.map((m) => {
            const mrna_id = m[2]
                .split("Exon table for ")[1]
                .split("and protein")[0]
                .replace("mRNA ", "")
                .trim();
            const protein = m[2].split("and protein")[1].trim();
            let exons = [];
            let cds = [];
            if (m[0].includes(mrna_id)) {
                const n_exons = parseInt(
                    m[0].split(`${mrna_id}, `)[1].split(" exons, ")[0]
                );
                exons = m
                    .slice(5, 5 + n_exons)
                    .map((e) => e.trim().replace(/\s+/g, ",").split(","));
                exons = exons.map((e) =>
                    e.filter((x) => x.includes("-")).length < 4
                        ? e[0].split("-")
                        : (e[0] + "-" + e[1]).split("-")
                );
                cds = exons
                    .filter((e) => e.length > 2)
                    .map((e) => [parseInt(e[2]), parseInt(e[3])]);
                cds = cds.map(c => c[0] < c[1] ? c : [c[1], c[0]]);
                exons = exons.map((e) => [parseInt(e[0]), parseInt(e[1])]);
                exons = exons.map(c => c[0] < c[1] ? c : [c[1], c[0]])
            }
            return { mrna_id, protein, exons, cds };
        });

        this.parsed = { gene, organism, gid, seq, strand, from, to, mrnas };
    }

    getParsed() {
        if (!this.parsed) this.parse();
        return this.parsed;
    }
}