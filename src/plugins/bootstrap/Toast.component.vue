<template>
  <div class="toasts">
    <div
      v-for="toast in toasts"
      :key="toast.id"
      :ref="toast.id"
      :id="toast.id"
      class="toast d-flex align-items-center border-0"
      :class="[
        toast.color === 'light' ? '' : 'text-white',
        `bg-${toast.color}`,
      ]"
      role="alert"
      aria-live="assertive"
      aria-atomic="true"
    >
      <div class="toast-body">
        <strong class="mr-2" v-if="toast.title">{{ toast.title }}</strong
        >{{ toast.text }}
      </div>
      <button
        type="button"
        :class="toast.color === 'light' ? '' : 'btn-close-white'"
        class="btn-close ml-auto mr-2 fix-multiline"
        data-dismiss="toast"
        aria-label="Close"
      ></button>
    </div>
  </div>
</template>
<script>
export default {
  props: {
    global: Boolean,
  },
  data: () => ({ toasts: [], index: 1 }),
  mounted() {
    if (this.global) {
      this.$toast_handler(this);
    }
  },
  updated() {
    this.toasts.forEach((toast) => {
      const el = this.$refs[toast.id];
      if (!toast.instance && el) {
        toast.instance = new this.$bootstrap.Toast(this.$refs[toast.id], {
          autohide: toast.delay > 0,
          delay: toast.delay * 1000,
        });
        toast.instance.show();
        el.addEventListener("hidden.bs.toast", (el) => {
          const index = this.toasts.findIndex((t) => t.id === el.target.id);
          const toast_ = this.toasts[index];
          toast_.instance.hide();
          this.toasts.splice(index, 1);
        });
      }
    });
  },
  methods: {
    notify(text, title, color = "secondary", delay = 5) {
      this.toasts.push({
        text,
        title,
        color,
        delay,
        id: `toast-${this.index++}`,
      });
    },
  },
};
</script>
<style scoped>
.toasts {
  width: 25rem;
  z-index: 100;
  position: fixed;
  top: 6rem;
  right: 8rem;
}
.fix-multiline {
  padding: .4rem;
}
</style>