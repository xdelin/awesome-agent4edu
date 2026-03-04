# Chart.js with Vue 3 (vue-chartjs) Examples

Complete examples for integrating Chart.js with Vue 3 using vue-chartjs.

## Installation

```bash
npm install vue-chartjs chart.js
```

## Basic Setup with Composition API

```vue
<!-- components/BarChart.vue -->
<template>
  <Bar :data="chartData" :options="chartOptions" />
</template>

<script setup lang="ts">
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend,
  type ChartData,
  type ChartOptions
} from 'chart.js';
import { Bar } from 'vue-chartjs';

// Register required components
ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
);

interface Props {
  chartData: ChartData<'bar'>;
  chartOptions?: ChartOptions<'bar'>;
}

const props = withDefaults(defineProps<Props>(), {
  chartOptions: () => ({
    responsive: true,
    maintainAspectRatio: true,
    plugins: {
      legend: { position: 'top' }
    }
  })
});
</script>
```

## Usage Example

```vue
<!-- views/Dashboard.vue -->
<template>
  <div class="p-6">
    <h1 class="text-2xl font-bold mb-4">Dashboard</h1>
    <div class="grid grid-cols-2 gap-6">
      <div class="bg-white p-4 rounded-lg shadow">
        <h2 class="text-lg font-semibold mb-2">Monthly Sales</h2>
        <BarChart :chart-data="salesData" />
      </div>
      <div class="bg-white p-4 rounded-lg shadow">
        <h2 class="text-lg font-semibold mb-2">Revenue Trend</h2>
        <LineChart :chart-data="revenueData" />
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref } from 'vue';
import BarChart from '../components/BarChart.vue';
import LineChart from '../components/LineChart.vue';

const salesData = ref({
  labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'],
  datasets: [
    {
      label: 'Sales 2024',
      data: [12, 19, 3, 5, 2, 3],
      backgroundColor: 'rgba(54, 162, 235, 0.5)',
      borderColor: 'rgb(54, 162, 235)',
      borderWidth: 1
    },
    {
      label: 'Sales 2023',
      data: [8, 15, 7, 12, 6, 9],
      backgroundColor: 'rgba(255, 99, 132, 0.5)',
      borderColor: 'rgb(255, 99, 132)',
      borderWidth: 1
    }
  ]
});

const revenueData = ref({
  labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'],
  datasets: [
    {
      label: 'Revenue',
      data: [1200, 1900, 300, 500, 200, 300],
      borderColor: 'rgb(75, 192, 192)',
      tension: 0.1
    }
  ]
});
</script>
```

## Reactive Updates with ref

```vue
<!-- components/ReactiveChart.vue -->
<template>
  <div>
    <Line :data="chartData" :options="chartOptions" />
    <button @click="addDataPoint" class="mt-4 px-4 py-2 bg-blue-500 text-white rounded">
      Add Data Point
    </button>
  </div>
</template>

<script setup lang="ts">
import { ref, computed } from 'vue';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip
} from 'chart.js';
import { Line } from 'vue-chartjs';

ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip
);

const labels = ref<string[]>(['Jan', 'Feb', 'Mar']);
const dataPoints = ref<number[]>([10, 20, 30]);

// Computed property ensures reactivity
const chartData = computed(() => ({
  labels: labels.value,
  datasets: [
    {
      label: 'Data',
      data: dataPoints.value,
      borderColor: 'rgb(75, 192, 192)',
      backgroundColor: 'rgba(75, 192, 192, 0.5)'
    }
  ]
}));

const chartOptions = {
  responsive: true,
  plugins: {
    legend: { position: 'top' as const }
  }
};

function addDataPoint() {
  const months = ['Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  const nextMonth = months[labels.value.length - 3] || `Point ${labels.value.length + 1}`;
  labels.value.push(nextMonth);
  dataPoints.value.push(Math.floor(Math.random() * 100));
}
</script>
```

## Live Data Updates

```vue
<!-- components/LiveChart.vue -->
<template>
  <div>
    <Line :data="chartData" :options="chartOptions" />
    <div class="mt-4 space-x-2">
      <button @click="startStream" :disabled="isStreaming" class="px-4 py-2 bg-green-500 text-white rounded disabled:opacity-50">
        Start
      </button>
      <button @click="stopStream" :disabled="!isStreaming" class="px-4 py-2 bg-red-500 text-white rounded disabled:opacity-50">
        Stop
      </button>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, onUnmounted, computed } from 'vue';
import { Line } from 'vue-chartjs';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement
} from 'chart.js';

ChartJS.register(CategoryScale, LinearScale, PointElement, LineElement);

const MAX_POINTS = 20;

const labels = ref<string[]>([]);
const dataPoints = ref<number[]>([]);
const isStreaming = ref(false);
let intervalId: number | null = null;

const chartData = computed(() => ({
  labels: labels.value,
  datasets: [
    {
      label: 'Live Data',
      data: dataPoints.value,
      borderColor: 'rgb(75, 192, 192)',
      fill: false,
      tension: 0.1
    }
  ]
}));

const chartOptions = {
  responsive: true,
  animation: { duration: 0 },
  scales: {
    y: { beginAtZero: true, max: 100 }
  }
};

function addPoint() {
  const now = new Date().toLocaleTimeString();
  const value = Math.floor(Math.random() * 100);

  labels.value.push(now);
  dataPoints.value.push(value);

  if (labels.value.length > MAX_POINTS) {
    labels.value.shift();
    dataPoints.value.shift();
  }
}

function startStream() {
  isStreaming.value = true;
  intervalId = window.setInterval(addPoint, 1000);
}

function stopStream() {
  isStreaming.value = false;
  if (intervalId) {
    clearInterval(intervalId);
    intervalId = null;
  }
}

onUnmounted(() => {
  stopStream();
});
</script>
```

## Data Fetching

### Vanilla Vue 3

```vue
<!-- components/FetchedChart.vue -->
<template>
  <div>
    <div v-if="loading" class="text-gray-500">Loading chart...</div>
    <div v-else-if="error" class="text-red-500">Error: {{ error }}</div>
    <Bar v-else :data="chartData" :options="chartOptions" />
  </div>
</template>

<script setup lang="ts">
import { ref, computed, onMounted } from 'vue';
import { Bar } from 'vue-chartjs';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip
} from 'chart.js';

ChartJS.register(CategoryScale, LinearScale, BarElement, Title, Tooltip);

interface SalesData {
  month: string;
  revenue: number;
}

const data = ref<SalesData[]>([]);
const loading = ref(true);
const error = ref<string | null>(null);

onMounted(async () => {
  try {
    const response = await fetch('/api/sales');
    if (!response.ok) throw new Error('Failed to fetch');
    data.value = await response.json();
  } catch (err) {
    error.value = err instanceof Error ? err.message : 'Unknown error';
  } finally {
    loading.value = false;
  }
});

const chartData = computed(() => ({
  labels: data.value.map((d) => d.month),
  datasets: [
    {
      label: 'Revenue',
      data: data.value.map((d) => d.revenue),
      backgroundColor: 'rgba(54, 162, 235, 0.5)',
      borderColor: 'rgb(54, 162, 235)',
      borderWidth: 1
    }
  ]
}));

const chartOptions = {
  responsive: true,
  plugins: {
    title: { display: true, text: 'Monthly Revenue' }
  }
};
</script>
```

### Nuxt 3 (useFetch)

```vue
<!-- components/FetchedChart.vue (Nuxt 3) -->
<template>
  <div>
    <div v-if="pending" class="text-gray-500">Loading chart...</div>
    <div v-else-if="error" class="text-red-500">Error: {{ error.message }}</div>
    <Bar v-else :data="chartData" :options="chartOptions" />
  </div>
</template>

<script setup lang="ts">
import { computed } from 'vue';
import { Bar } from 'vue-chartjs';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip
} from 'chart.js';

ChartJS.register(CategoryScale, LinearScale, BarElement, Title, Tooltip);

interface SalesData {
  month: string;
  revenue: number;
}

// Nuxt 3 auto-imported composable
const { data, pending, error } = await useFetch<SalesData[]>('/api/sales');

const chartData = computed(() => ({
  labels: data.value?.map((d) => d.month) ?? [],
  datasets: [
    {
      label: 'Revenue',
      data: data.value?.map((d) => d.revenue) ?? [],
      backgroundColor: 'rgba(54, 162, 235, 0.5)',
      borderColor: 'rgb(54, 162, 235)',
      borderWidth: 1
    }
  ]
}));

const chartOptions = {
  responsive: true,
  plugins: {
    title: { display: true, text: 'Monthly Revenue' }
  }
};
</script>
```

## Reusable Chart Composable

```typescript
// composables/useChart.ts
import { ref, computed, onUnmounted, type Ref } from 'vue';
import type { ChartData } from 'chart.js';

interface UseChartOptions {
  maxPoints?: number;
}

export function useChartData<T extends 'bar' | 'line'>(
  initialData: ChartData<T>,
  options: UseChartOptions = {}
) {
  const { maxPoints = 50 } = options;

  const labels = ref<string[]>([...(initialData.labels as string[])]);
  const datasets = ref(
    initialData.datasets.map((ds) => ({
      ...ds,
      data: [...(ds.data as number[])]
    }))
  );

  const chartData = computed<ChartData<T>>(() => ({
    labels: labels.value,
    datasets: datasets.value
  })) as Ref<ChartData<T>>;

  function addPoint(label: string, values: number[]) {
    labels.value.push(label);
    datasets.value.forEach((ds, i) => {
      (ds.data as number[]).push(values[i] ?? 0);
    });

    if (labels.value.length > maxPoints) {
      labels.value.shift();
      datasets.value.forEach((ds) => {
        (ds.data as number[]).shift();
      });
    }
  }

  function updatePoint(datasetIndex: number, pointIndex: number, value: number) {
    (datasets.value[datasetIndex].data as number[])[pointIndex] = value;
  }

  function reset() {
    labels.value = [...(initialData.labels as string[])];
    datasets.value = initialData.datasets.map((ds) => ({
      ...ds,
      data: [...(ds.data as number[])]
    }));
  }

  return {
    chartData,
    addPoint,
    updatePoint,
    reset
  };
}
```

## Usage with Composable

```vue
<!-- views/Analytics.vue -->
<template>
  <div>
    <Line :data="chartData" :options="options" />
    <div class="mt-4">
      <span class="text-sm text-gray-500">Connected: {{ isConnected ? 'Yes' : 'No' }}</span>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, onUnmounted } from 'vue';
import { Line } from 'vue-chartjs';
import { useChartData } from '../composables/useChart';

const { chartData, addPoint } = useChartData({
  labels: [],
  datasets: [
    {
      label: 'Metric A',
      data: [],
      borderColor: 'rgb(75, 192, 192)'
    },
    {
      label: 'Metric B',
      data: [],
      borderColor: 'rgb(255, 99, 132)'
    }
  ]
}, { maxPoints: 30 });

const options = {
  responsive: true,
  animation: { duration: 0 }
};

const isConnected = ref(false);
let ws: WebSocket | null = null;

onMounted(() => {
  ws = new WebSocket('wss://api.example.com/metrics');

  ws.onopen = () => {
    isConnected.value = true;
  };

  ws.onmessage = (event) => {
    const { timestamp, metricA, metricB } = JSON.parse(event.data);
    addPoint(timestamp, [metricA, metricB]);
  };

  ws.onclose = () => {
    isConnected.value = false;
  };
});

onUnmounted(() => {
  ws?.close();
});
</script>
```

## Options API Example

```vue
<!-- components/OptionsApiChart.vue -->
<template>
  <div>
    <Bar :data="chartData" :options="chartOptions" />
    <button @click="updateData" class="mt-4 px-4 py-2 bg-blue-500 text-white rounded">
      Randomize Data
    </button>
  </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js';
import { Bar } from 'vue-chartjs';

ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
);

export default defineComponent({
  name: 'OptionsApiChart',
  components: { Bar },

  data() {
    return {
      chartData: {
        labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May'],
        datasets: [
          {
            label: 'Sales',
            data: [40, 20, 12, 39, 10],
            backgroundColor: 'rgba(54, 162, 235, 0.5)'
          }
        ]
      },
      chartOptions: {
        responsive: true,
        plugins: {
          legend: { position: 'top' }
        }
      }
    };
  },

  methods: {
    updateData() {
      this.chartData = {
        ...this.chartData,
        datasets: [
          {
            ...this.chartData.datasets[0],
            data: this.chartData.datasets[0].data.map(() =>
              Math.floor(Math.random() * 100)
            )
          }
        ]
      };
    }
  }
});
</script>
```

## Chart with Props Validation

```vue
<!-- components/ConfigurableChart.vue -->
<template>
  <component :is="chartComponent" :data="data" :options="mergedOptions" />
</template>

<script setup lang="ts">
import { computed } from 'vue';
import {
  Chart as ChartJS,
  registerables,
  type ChartType,
  type ChartData,
  type ChartOptions
} from 'chart.js';
import { Bar, Line, Pie, Doughnut, Radar, PolarArea } from 'vue-chartjs';

ChartJS.register(...registerables);

type SupportedType = 'bar' | 'line' | 'pie' | 'doughnut' | 'radar' | 'polarArea';

interface Props {
  type: SupportedType;
  data: ChartData;
  options?: ChartOptions;
  height?: number;
}

const props = withDefaults(defineProps<Props>(), {
  options: () => ({}),
  height: 400
});

const chartComponents = {
  bar: Bar,
  line: Line,
  pie: Pie,
  doughnut: Doughnut,
  radar: Radar,
  polarArea: PolarArea
} as const;

const chartComponent = computed(() => chartComponents[props.type]);

const mergedOptions = computed<ChartOptions>(() => ({
  responsive: true,
  maintainAspectRatio: true,
  ...props.options
}));
</script>
```

## Provide/Inject for Theme

```typescript
// plugins/chartTheme.ts
import { inject, provide, type InjectionKey } from 'vue';
import type { ChartOptions } from 'chart.js';

interface ChartTheme {
  colors: string[];
  fontFamily: string;
  gridColor: string;
}

const ChartThemeKey: InjectionKey<ChartTheme> = Symbol('ChartTheme');

const defaultTheme: ChartTheme = {
  colors: [
    'rgb(54, 162, 235)',
    'rgb(255, 99, 132)',
    'rgb(75, 192, 192)',
    'rgb(255, 205, 86)'
  ],
  fontFamily: 'Inter, system-ui, sans-serif',
  gridColor: 'rgba(0, 0, 0, 0.1)'
};

export function provideChartTheme(theme: Partial<ChartTheme> = {}) {
  provide(ChartThemeKey, { ...defaultTheme, ...theme });
}

export function useChartTheme(): ChartOptions {
  const theme = inject(ChartThemeKey, defaultTheme);

  return {
    plugins: {
      legend: {
        labels: {
          font: { family: theme.fontFamily }
        }
      }
    },
    scales: {
      x: {
        grid: { color: theme.gridColor },
        ticks: { font: { family: theme.fontFamily } }
      },
      y: {
        grid: { color: theme.gridColor },
        ticks: { font: { family: theme.fontFamily } }
      }
    }
  };
}
```

## Testing Vue Charts

```typescript
// __tests__/BarChart.spec.ts
import { describe, it, expect, vi, beforeAll } from 'vitest';
import { mount } from '@vue/test-utils';
import BarChart from '../components/BarChart.vue';

// Mock canvas for testing
beforeAll(() => {
  HTMLCanvasElement.prototype.getContext = vi.fn();
});

describe('BarChart', () => {
  const mockData = {
    labels: ['A', 'B', 'C'],
    datasets: [{ label: 'Test', data: [1, 2, 3] }]
  };

  it('renders without crashing', () => {
    const wrapper = mount(BarChart, {
      props: { chartData: mockData }
    });
    expect(wrapper.find('canvas').exists()).toBe(true);
  });

  it('accepts custom options', () => {
    const options = {
      plugins: { title: { display: true, text: 'Test' } }
    };

    const wrapper = mount(BarChart, {
      props: { chartData: mockData, chartOptions: options }
    });

    expect(wrapper.find('canvas').exists()).toBe(true);
  });

  it('updates when data changes', async () => {
    const wrapper = mount(BarChart, {
      props: { chartData: mockData }
    });

    await wrapper.setProps({
      chartData: {
        ...mockData,
        datasets: [{ label: 'Updated', data: [4, 5, 6] }]
      }
    });

    expect(wrapper.find('canvas').exists()).toBe(true);
  });
});
```
