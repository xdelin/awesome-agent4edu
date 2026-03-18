# Chart.js with React (react-chartjs-2) Examples

Complete examples for integrating Chart.js with React using react-chartjs-2.

## Installation

```bash
npm install react-chartjs-2 chart.js
```

## Basic Setup with Tree-Shaking

```tsx
// components/BarChart.tsx
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend,
  ChartOptions,
  ChartData
} from 'chart.js';
import { Bar } from 'react-chartjs-2';

// Register required components
ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
);

interface BarChartProps {
  data: ChartData<'bar'>;
  options?: ChartOptions<'bar'>;
}

export function BarChart({ data, options }: BarChartProps) {
  const defaultOptions: ChartOptions<'bar'> = {
    responsive: true,
    maintainAspectRatio: true,
    plugins: {
      legend: { position: 'top' },
      title: { display: false }
    },
    ...options
  };

  return <Bar data={data} options={defaultOptions} />;
}
```

## Usage Example

```tsx
// pages/Dashboard.tsx
import { BarChart } from '../components/BarChart';

export function Dashboard() {
  const salesData = {
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
  };

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Dashboard</h1>
      <div className="bg-white p-4 rounded-lg shadow">
        <BarChart data={salesData} />
      </div>
    </div>
  );
}
```

## Chart Reference and Click Events

```tsx
// components/InteractiveChart.tsx
import { useRef, MouseEvent } from 'react';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Tooltip,
  ChartData,
  ActiveElement,
  ChartEvent
} from 'chart.js';
import { Bar, getElementAtEvent } from 'react-chartjs-2';

ChartJS.register(CategoryScale, LinearScale, BarElement, Tooltip);

interface InteractiveChartProps {
  data: ChartData<'bar'>;
  onBarClick?: (datasetIndex: number, index: number, value: number) => void;
}

export function InteractiveChart({ data, onBarClick }: InteractiveChartProps) {
  const chartRef = useRef<ChartJS<'bar'>>(null);

  const handleClick = (event: MouseEvent<HTMLCanvasElement>) => {
    const chart = chartRef.current;
    if (!chart) return;

    const elements = getElementAtEvent(chart, event);
    if (elements.length > 0) {
      const { datasetIndex, index } = elements[0];
      const value = data.datasets[datasetIndex].data[index] as number;
      onBarClick?.(datasetIndex, index, value);
    }
  };

  return <Bar ref={chartRef} data={data} onClick={handleClick} />;
}
```

## Dynamic Data Updates with useState

```tsx
// components/LiveChart.tsx
import { useState, useEffect, useCallback } from 'react';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  ChartData
} from 'chart.js';
import { Line } from 'react-chartjs-2';

ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip
);

const MAX_POINTS = 20;

export function LiveChart() {
  const [chartData, setChartData] = useState<ChartData<'line'>>({
    labels: [],
    datasets: [
      {
        label: 'Live Data',
        data: [],
        borderColor: 'rgb(75, 192, 192)',
        backgroundColor: 'rgba(75, 192, 192, 0.5)',
        tension: 0.1
      }
    ]
  });

  const addDataPoint = useCallback(() => {
    const now = new Date().toLocaleTimeString();
    const value = Math.floor(Math.random() * 100);

    setChartData((prev) => ({
      ...prev,
      labels: [...prev.labels!.slice(-MAX_POINTS + 1), now],
      datasets: [
        {
          ...prev.datasets[0],
          data: [...(prev.datasets[0].data as number[]).slice(-MAX_POINTS + 1), value]
        }
      ]
    }));
  }, []);

  useEffect(() => {
    const interval = setInterval(addDataPoint, 1000);
    return () => clearInterval(interval);
  }, [addDataPoint]);

  return (
    <Line
      data={chartData}
      options={{
        responsive: true,
        animation: { duration: 0 },
        scales: {
          y: { beginAtZero: true, max: 100 }
        }
      }}
    />
  );
}
```

## Data Fetching with useEffect

```tsx
// components/FetchedChart.tsx
import { useState, useEffect } from 'react';
import { Line } from 'react-chartjs-2';
import { ChartData } from 'chart.js';

interface SalesDataPoint {
  month: string;
  revenue: number;
}

export function FetchedChart() {
  const [chartData, setChartData] = useState<ChartData<'line'> | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    async function fetchData() {
      try {
        const response = await fetch('/api/sales');
        const data: SalesDataPoint[] = await response.json();

        setChartData({
          labels: data.map((d) => d.month),
          datasets: [
            {
              label: 'Revenue',
              data: data.map((d) => d.revenue),
              borderColor: 'rgb(54, 162, 235)',
              backgroundColor: 'rgba(54, 162, 235, 0.5)'
            }
          ]
        });
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to fetch data');
      } finally {
        setLoading(false);
      }
    }

    fetchData();
  }, []);

  if (loading) return <div>Loading chart...</div>;
  if (error) return <div>Error: {error}</div>;
  if (!chartData) return null;

  return <Line data={chartData} />;
}
```

## Reusable Chart Component with TypeScript

```tsx
// components/Chart.tsx
import {
  Chart as ChartJS,
  ChartType,
  ChartData,
  ChartOptions,
  registerables
} from 'chart.js';
import {
  Chart as ChartComponent,
  Bar,
  Line,
  Pie,
  Doughnut,
  Radar,
  PolarArea
} from 'react-chartjs-2';

// Register all Chart.js components
ChartJS.register(...registerables);

type SupportedChartType = 'bar' | 'line' | 'pie' | 'doughnut' | 'radar' | 'polarArea';

interface ChartProps<T extends SupportedChartType> {
  type: T;
  data: ChartData<T>;
  options?: ChartOptions<T>;
  height?: number;
  width?: number;
}

const chartComponents = {
  bar: Bar,
  line: Line,
  pie: Pie,
  doughnut: Doughnut,
  radar: Radar,
  polarArea: PolarArea
} as const;

export function Chart<T extends SupportedChartType>({
  type,
  data,
  options,
  height,
  width
}: ChartProps<T>) {
  const ChartElement = chartComponents[type];

  return (
    <ChartElement
      data={data as ChartData<typeof type>}
      options={options as ChartOptions<typeof type>}
      height={height}
      width={width}
    />
  );
}
```

## Custom Hook for Chart Data

```tsx
// hooks/useChartData.ts
import { useState, useCallback } from 'react';
import { ChartData } from 'chart.js';

interface UseChartDataOptions<T extends 'bar' | 'line'> {
  initialData: ChartData<T>;
  maxPoints?: number;
}

export function useChartData<T extends 'bar' | 'line'>({
  initialData,
  maxPoints = 50
}: UseChartDataOptions<T>) {
  const [data, setData] = useState<ChartData<T>>(initialData);

  const addPoint = useCallback(
    (label: string, values: number[]) => {
      setData((prev) => {
        const newLabels = [...(prev.labels || []), label].slice(-maxPoints);
        const newDatasets = prev.datasets.map((dataset, i) => ({
          ...dataset,
          data: [...(dataset.data as number[]), values[i] ?? 0].slice(-maxPoints)
        }));

        return { ...prev, labels: newLabels, datasets: newDatasets };
      });
    },
    [maxPoints]
  );

  const updatePoint = useCallback((datasetIndex: number, pointIndex: number, value: number) => {
    setData((prev) => {
      const newDatasets = [...prev.datasets];
      const newData = [...(newDatasets[datasetIndex].data as number[])];
      newData[pointIndex] = value;
      newDatasets[datasetIndex] = { ...newDatasets[datasetIndex], data: newData };

      return { ...prev, datasets: newDatasets };
    });
  }, []);

  const reset = useCallback(() => {
    setData(initialData);
  }, [initialData]);

  return { data, addPoint, updatePoint, reset, setData };
}
```

## Usage with Custom Hook

```tsx
// pages/Analytics.tsx
import { useEffect } from 'react';
import { Line } from 'react-chartjs-2';
import { useChartData } from '../hooks/useChartData';

export function Analytics() {
  const { data, addPoint } = useChartData({
    initialData: {
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
    },
    maxPoints: 30
  });

  useEffect(() => {
    const ws = new WebSocket('wss://api.example.com/metrics');

    ws.onmessage = (event) => {
      const { timestamp, metricA, metricB } = JSON.parse(event.data);
      addPoint(timestamp, [metricA, metricB]);
    };

    return () => ws.close();
  }, [addPoint]);

  return <Line data={data} options={{ responsive: true }} />;
}
```

## Chart with Context Provider

```tsx
// context/ChartContext.tsx
import { createContext, useContext, ReactNode } from 'react';
import { ChartOptions } from 'chart.js';

interface ChartTheme {
  colors: string[];
  fontFamily: string;
  gridColor: string;
}

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

const ChartContext = createContext<ChartTheme>(defaultTheme);

export function ChartProvider({
  children,
  theme = defaultTheme
}: {
  children: ReactNode;
  theme?: ChartTheme;
}) {
  return <ChartContext.Provider value={theme}>{children}</ChartContext.Provider>;
}

export function useChartTheme(): ChartOptions {
  const theme = useContext(ChartContext);

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

## Testing React Charts

```tsx
// __tests__/BarChart.test.tsx
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { BarChart } from '../components/BarChart';

// Mock canvas for testing
beforeAll(() => {
  HTMLCanvasElement.prototype.getContext = jest.fn();
});

describe('BarChart', () => {
  const mockData = {
    labels: ['A', 'B', 'C'],
    datasets: [{ label: 'Test', data: [1, 2, 3] }]
  };

  it('renders without crashing', () => {
    render(<BarChart data={mockData} />);
    expect(screen.getByRole('img')).toBeInTheDocument();
  });

  it('applies custom options', () => {
    const options = {
      plugins: { title: { display: true, text: 'Test Chart' } }
    };

    render(<BarChart data={mockData} options={options} />);
    // Chart.js renders to canvas, verify container exists
    expect(document.querySelector('canvas')).toBeInTheDocument();
  });
});
```
