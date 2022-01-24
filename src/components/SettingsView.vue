<template>
  <Display subtitle ico="gear">Settings</Display>
  <div>
    <div class="mb-3">
      <label for="exampleInputEmail1" class="form-label">Email address</label>
      <input
        type="email"
        class="form-control"
        id="exampleInputEmail1"
        aria-describedby="emailHelp"
        v-model="email"
      />
      <div id="emailHelp" class="form-text">
        Email sera usado para informar o usuario da API.
      </div>
    </div>

    <div class="mb-3">
      <label for="exampleInputPassword1" class="form-label"
        >NCBI API Host</label
      >
      <input
        type="text"
        class="form-control"
        id="exampleInputPassword1"
        v-model="ncbi"
      />
    </div>

    <div class="mb-3">
      <label for="exampleInputPassword1" class="form-label"
        >NCBI API Key</label
      >
      <input
        type="text"
        class="form-control"
        id="exampleInputPassword1"
        v-model="api_key"
      />
      <div id="emailHelp" class="form-text">
        Visite this <a target="_blank" href="https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/">link for more details</a>.
      </div>
    </div>

    <div class="mb-3">
      <label for="exampleInputPassword1" class="form-label"
        >Interpro API Host</label
      >
      <input
        type="text"
        class="form-control"
        id="exampleInputPassword1"
        v-model="interpro"
      />
    </div>

    <Button ico="save" @click="save()" :disabled="saved">Save</Button>
  </div>
</template>

<script>
import NCBI from "../core/api/NCBI";
import Interpro from "../core/api/Interpro";
export default {
  data() {
    return {
      email: "",
      ncbi: "",
      interpro: "",
      api_key: "",
      saved: false,
    };
  },
  mounted() {
    this.email = this.storage.get("email");
    ///http://192.168.64.3:8088/efetch
    this.ncbi = this.storage.get("ncbi", NCBI.HOST);
    this.api_key = this.storage.get("API_KEY", NCBI.API_KEY);
    this.interpro = this.storage.get("interpro", Interpro.INTERPRO_HOST);
  },
  methods: {
    save() {
      this.storage.set("email", this.email);
      this.storage.set("ncbi", (NCBI.HOST = this.ncbi));
      this.storage.set("API_KEY", (NCBI.API_KEY = this.api_key));
      this.storage.set("interpro", (Interpro.INTERPRO_HOST = this.interpro));
      this.saved = true;
    },
  },
};
</script>

<style>
</style>