<template>
  <div class="progress" :style="`${'height: ' + height + ';'}`">
    <div
      class="progress-bar"
      :class="classObj"
      role="progressbar"
      :style="`${'width: ' + width + '%;'}`"
      :aria-valuenow="width"
      :aria-valuemin="min"
      :aria-valuemax="max"
    >
      <span v-if="label">{{ width }}%</span>
    </div>
  </div>
</template>
<script>
export default {
  props: {
    strip: {
      type: Boolean,
      default: true,
    },
    label: {
      type: Boolean,
      default: false,
    },
    anime: {
      type: Boolean,
      default: true,
    },
    min: {
      type: Number,
      default: 0,
    },
    max: {
      type: Number,
      default: 100,
    },
    height: {
      type: String,
      default: "unset",
    },
    color: {
      type: String,
      default: "primary",
    },
  },
  computed: {
    classObj: (t) => {
      const pre = {
        "progress-bar-striped": t.strip,
        "progress-bar-animated": t.anime,
      };
      pre["bg-" + t.color] = true;
      return pre;
    },
  },
  data: () => ({ width: 70 }),
  methods: {
    setValue(value) {
      this.setPercent(parseInt((100 * value) / this.max));
    },
    setPercent(percent) {
      this.width = percent;
    },
  },
};
</script>