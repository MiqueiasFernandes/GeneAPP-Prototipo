<template>
  importt view <Button @click="getGeneTable()">carregar</Button>
  <ProgressBar
    v-if="genes.length > 0 && !carregado"
    ref="progress"
    label
    :min="0"
    :max="genes.length"
  ></ProgressBar>

  <div v-if="carregado">
    <ol v-for="gene of genes" :key="gene.name">
      <li>
        <b>{{ gene.name }} {{ gene.status }}</b
        ><br />
        {{ gene.getParsed() }}
      </li>
    </ol>
  </div>
</template>

<script>
import axios from "axios";
const efetch_host = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
const efetch_query = "db=gene&rettype=gene_table&retmode=text&id=";
import Isoforma from "../core/locus/Isoforma";
import Exon from "../core/locus/Exon";
import CDS from "../core/locus/CDS";
import Gene from "../core/locus/Gene";

class GeneTable {
  data = null;
  parsed = null;
  status = "unload";
  constructor(name) {
    this.name = name;
    axios
      .get(`${efetch_host}?${efetch_query}${name}`)
      .then((resp) => {
        this.status = "loaded";
        this.data = resp.data;
      })
      .catch(() => (this.status = "error"));
  }

  parse() {
    const parsed = this.data.split("\n");
    const gene = parsed[0].split("[")[0];
    const organism = parsed[0].split("[")[1].split("]")[0];
    const gid = parsed[1].split("Gene ID:")[1].split(",")[0].trim();
    const ref = parsed[4].split("Primary Assembly")[1];
    const seq = ref.trimStart().split(" ")[0].trim();
    const strand = !ref.includes("(minus strand)");
    const from = parseInt(ref.split("from:")[1].split("to")[0].trim());
    const to = parseInt(ref.split("to: ")[1].split(" ")[0]);
    let mrnas = [];
    parsed
      .slice(5)
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
        exons = m.slice(5, 5+n_exons).map(e => e.trim().replace(/\s+/g, ',').split(','));
        exons = exons.map(e => e.filter(x => x.includes('-')).length < 4 ? e[0].split('-') : (e[0]+'-'+e[1]).split('-'));
        cds = exons.filter(e => e.length > 2).map(e => [parseInt(e[2]), parseInt(e[3])])
        exons = exons.map(e => [parseInt(e[0]), parseInt(e[1])])
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

export default {
  data() {
    return {
      carregado: false,
      genes: [],
    };
  },
  methods: {
    getGeneTable() {
      this.genes = "LOC101258920,LOC101245278,LOC104644837,57502"
        .split(",")
        .map((g) => new GeneTable(g));

      const ti = setInterval(() => {
        const qtd = this.genes.filter((g) => g.status === "loaded").length;
        this.$refs.progress.setValue(qtd);
        if (qtd === this.genes.length || !this.$refs.progress) {
          clearInterval(ti);
          this.loaded();
        }
      }, 1000);
    },

    loaded() {
      this.carregado = true;
    },

    getGenes() {
      return this.genes.map(g => this.parseGene(g));
    },

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
          isoforma.addCDS(new CDS(null, c[0], c[1], strand))
        );
        return isoforma;
      });
      isos.forEach((iso) => _gene.addIsoforma(iso));
      return _gene;
    },
  },
};
</script>

<style>
</style>