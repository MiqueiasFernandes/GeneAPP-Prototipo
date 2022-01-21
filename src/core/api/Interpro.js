import axios from "axios";

const INTERPRO_HOST = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
const WAIT = 5000

export default class Interpro {
    static runs = {};
    static INSTANCE = null;
    static api() {
        if (!Interpro.INSTANCE)
            Interpro.INSTANCE = new Interpro();
        return Interpro.INSTANCE;
    }

    post(id, sequence, email) {
        Interpro.runs[id] = { id, job: null, status: 0, cbk: null, job_status: null }
        const str = `email=${email}&goterms=true&pathways=true&appl=PfamA&sequence=`;
        axios
            .post(`${INTERPRO_HOST}/run`, str + sequence)
            .then(r => {
                Interpro.runs[id].job = r.data;
                Interpro.runs[id].status = 1;
                const int = setInterval(() => {
                    if (Interpro.runs[id].status > 1) {
                        clearInterval(int);
                        Interpro.api().load(Interpro.runs[id].job, id)
                    } else {
                        Interpro.api().get(Interpro.runs[id].job, id)
                    }
                }, WAIT)
            });
        return Interpro.runs[id];
    }

    get(job, id) {
        axios
            .get(`${INTERPRO_HOST}/status/${job}`)
            .then((r) => {
                Interpro.runs[id].job_status = r.data;
                if (r.data === "FINISHED") {
                    Interpro.runs[id].status++;
                }
            })
    }

    load(job, id) {
        axios
            .get(`${INTERPRO_HOST}/result/${job}/json`)
            .then((r) => {
                Interpro.runs[id].data = r.data;
                Interpro.runs[id].status++;
                if (Interpro.runs[id].cbk) {
                    Interpro.runs[id].cbk(Interpro.runs[id], Interpro.runs)
                }
            })
    }
}