# Chart.js with Angular (ng2-charts) Examples

Complete examples for integrating Chart.js with Angular using ng2-charts.

## Installation

```bash
npm install ng2-charts chart.js
```

## Standalone Component Setup (Angular 17+)

```typescript
// components/bar-chart.component.ts
import { Component, Input, OnChanges, SimpleChanges } from '@angular/core';
import { BaseChartDirective } from 'ng2-charts';
import {
  Chart,
  ChartConfiguration,
  ChartData,
  ChartType,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js';

// Register required components
Chart.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
);

@Component({
  selector: 'app-bar-chart',
  standalone: true,
  imports: [BaseChartDirective],
  template: `
    <canvas baseChart
      [data]="chartData"
      [options]="chartOptions"
      [type]="chartType">
    </canvas>
  `
})
export class BarChartComponent {
  @Input() chartData!: ChartData<'bar'>;
  @Input() chartOptions: ChartConfiguration<'bar'>['options'] = {
    responsive: true,
    maintainAspectRatio: true,
    plugins: {
      legend: { position: 'top' }
    }
  };

  chartType: ChartType = 'bar';
}
```

## Usage Example

```typescript
// pages/dashboard.component.ts
import { Component } from '@angular/core';
import { BarChartComponent } from '../components/bar-chart.component';
import { LineChartComponent } from '../components/line-chart.component';
import { ChartData } from 'chart.js';

@Component({
  selector: 'app-dashboard',
  standalone: true,
  imports: [BarChartComponent, LineChartComponent],
  template: `
    <div class="p-6">
      <h1 class="text-2xl font-bold mb-4">Dashboard</h1>
      <div class="grid grid-cols-2 gap-6">
        <div class="bg-white p-4 rounded-lg shadow">
          <h2 class="text-lg font-semibold mb-2">Monthly Sales</h2>
          <app-bar-chart [chartData]="salesData"></app-bar-chart>
        </div>
        <div class="bg-white p-4 rounded-lg shadow">
          <h2 class="text-lg font-semibold mb-2">Revenue Trend</h2>
          <app-line-chart [chartData]="revenueData"></app-line-chart>
        </div>
      </div>
    </div>
  `
})
export class DashboardComponent {
  salesData: ChartData<'bar'> = {
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

  revenueData: ChartData<'line'> = {
    labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'],
    datasets: [
      {
        label: 'Revenue',
        data: [1200, 1900, 300, 500, 200, 300],
        borderColor: 'rgb(75, 192, 192)',
        tension: 0.1
      }
    ]
  };
}
```

## Module-Based Setup (Legacy)

```typescript
// app.module.ts
import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { NgChartsModule } from 'ng2-charts';
import { AppComponent } from './app.component';

@NgModule({
  declarations: [AppComponent],
  imports: [
    BrowserModule,
    NgChartsModule
  ],
  bootstrap: [AppComponent]
})
export class AppModule {}
```

```typescript
// chart.component.ts (NgModule-based)
import { Component } from '@angular/core';
import { ChartConfiguration, ChartType } from 'chart.js';

@Component({
  selector: 'app-chart',
  template: `
    <canvas baseChart
      [data]="barChartData"
      [options]="barChartOptions"
      [type]="barChartType">
    </canvas>
  `
})
export class ChartComponent {
  barChartType: ChartType = 'bar';

  barChartData: ChartConfiguration<'bar'>['data'] = {
    labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May'],
    datasets: [{
      label: 'Sales',
      data: [12, 19, 3, 5, 2],
      backgroundColor: 'rgba(54, 162, 235, 0.5)'
    }]
  };

  barChartOptions: ChartConfiguration<'bar'>['options'] = {
    responsive: true,
    plugins: {
      legend: { position: 'top' }
    }
  };
}
```

## Dynamic Data Updates

```typescript
// components/dynamic-chart.component.ts
import { Component, ViewChild, OnDestroy } from '@angular/core';
import { BaseChartDirective } from 'ng2-charts';
import {
  Chart,
  ChartConfiguration,
  ChartData,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement
} from 'chart.js';

Chart.register(CategoryScale, LinearScale, PointElement, LineElement);

@Component({
  selector: 'app-dynamic-chart',
  standalone: true,
  imports: [BaseChartDirective],
  template: `
    <canvas baseChart
      [data]="chartData"
      [options]="chartOptions"
      type="line">
    </canvas>
    <div class="mt-4 space-x-2">
      <button (click)="addData()" class="px-4 py-2 bg-blue-500 text-white rounded">
        Add Data
      </button>
      <button (click)="updateData()" class="px-4 py-2 bg-green-500 text-white rounded">
        Update Data
      </button>
      <button (click)="resetData()" class="px-4 py-2 bg-gray-500 text-white rounded">
        Reset
      </button>
    </div>
  `
})
export class DynamicChartComponent implements OnDestroy {
  @ViewChild(BaseChartDirective) chart?: BaseChartDirective;

  private readonly maxPoints = 10;

  chartData: ChartData<'line'> = {
    labels: ['Jan', 'Feb', 'Mar'],
    datasets: [{
      label: 'Data',
      data: [10, 20, 30],
      borderColor: 'rgb(75, 192, 192)',
      tension: 0.1
    }]
  };

  chartOptions: ChartConfiguration<'line'>['options'] = {
    responsive: true,
    animation: { duration: 300 }
  };

  addData(): void {
    const labels = this.chartData.labels as string[];
    const data = this.chartData.datasets[0].data as number[];

    labels.push(`Point ${labels.length + 1}`);
    data.push(Math.floor(Math.random() * 100));

    if (labels.length > this.maxPoints) {
      labels.shift();
      data.shift();
    }

    this.chart?.update();
  }

  updateData(): void {
    const data = this.chartData.datasets[0].data as number[];
    for (let i = 0; i < data.length; i++) {
      data[i] = Math.floor(Math.random() * 100);
    }
    this.chart?.update();
  }

  resetData(): void {
    this.chartData = {
      labels: ['Jan', 'Feb', 'Mar'],
      datasets: [{
        ...this.chartData.datasets[0],
        data: [10, 20, 30]
      }]
    };
    this.chart?.update();
  }

  ngOnDestroy(): void {
    // Chart cleanup is handled by ng2-charts
  }
}
```

## Live Data with RxJS

```typescript
// components/live-chart.component.ts
import { Component, OnInit, OnDestroy, ViewChild } from '@angular/core';
import { BaseChartDirective } from 'ng2-charts';
import { interval, Subject, takeUntil } from 'rxjs';
import { Chart, ChartData, CategoryScale, LinearScale, PointElement, LineElement } from 'chart.js';

Chart.register(CategoryScale, LinearScale, PointElement, LineElement);

@Component({
  selector: 'app-live-chart',
  standalone: true,
  imports: [BaseChartDirective],
  template: `
    <canvas baseChart
      [data]="chartData"
      [options]="chartOptions"
      type="line">
    </canvas>
    <div class="mt-4 space-x-2">
      <button (click)="startStream()" [disabled]="isStreaming"
        class="px-4 py-2 bg-green-500 text-white rounded disabled:opacity-50">
        Start
      </button>
      <button (click)="stopStream()" [disabled]="!isStreaming"
        class="px-4 py-2 bg-red-500 text-white rounded disabled:opacity-50">
        Stop
      </button>
    </div>
  `
})
export class LiveChartComponent implements OnInit, OnDestroy {
  @ViewChild(BaseChartDirective) chart?: BaseChartDirective;

  private destroy$ = new Subject<void>();
  private readonly maxPoints = 20;

  isStreaming = false;

  chartData: ChartData<'line'> = {
    labels: [],
    datasets: [{
      label: 'Live Data',
      data: [],
      borderColor: 'rgb(75, 192, 192)',
      fill: false
    }]
  };

  chartOptions = {
    responsive: true,
    animation: { duration: 0 },
    scales: {
      y: { beginAtZero: true, max: 100 }
    }
  };

  ngOnInit(): void {
    // Optional: auto-start
  }

  ngOnDestroy(): void {
    this.destroy$.next();
    this.destroy$.complete();
  }

  startStream(): void {
    this.isStreaming = true;

    interval(1000)
      .pipe(takeUntil(this.destroy$))
      .subscribe(() => {
        if (!this.isStreaming) return;
        this.addDataPoint();
      });
  }

  stopStream(): void {
    this.isStreaming = false;
    this.destroy$.next();
    this.destroy$ = new Subject<void>();
  }

  private addDataPoint(): void {
    const labels = this.chartData.labels as string[];
    const data = this.chartData.datasets[0].data as number[];

    const now = new Date().toLocaleTimeString();
    const value = Math.floor(Math.random() * 100);

    labels.push(now);
    data.push(value);

    if (labels.length > this.maxPoints) {
      labels.shift();
      data.shift();
    }

    this.chart?.update();
  }
}
```

## Data Fetching with HttpClient

### App Configuration (Angular 17+)

```typescript
// app.config.ts
import { ApplicationConfig } from '@angular/core';
import { provideHttpClient } from '@angular/common/http';

export const appConfig: ApplicationConfig = {
  providers: [
    provideHttpClient()
  ]
};
```

### Service

```typescript
// services/chart-data.service.ts
import { Injectable, inject } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable, map } from 'rxjs';
import { ChartData } from 'chart.js';

interface SalesDataPoint {
  month: string;
  revenue: number;
}

@Injectable({ providedIn: 'root' })
export class ChartDataService {
  private http = inject(HttpClient);

  getSalesData(): Observable<ChartData<'bar'>> {
    return this.http.get<SalesDataPoint[]>('/api/sales').pipe(
      map(data => ({
        labels: data.map(d => d.month),
        datasets: [{
          label: 'Revenue',
          data: data.map(d => d.revenue),
          backgroundColor: 'rgba(54, 162, 235, 0.5)',
          borderColor: 'rgb(54, 162, 235)',
          borderWidth: 1
        }]
      }))
    );
  }
}
```

### Component

```typescript
// components/fetched-chart.component.ts
import { Component, inject, OnInit } from '@angular/core';
import { CommonModule, AsyncPipe } from '@angular/common';
import { BaseChartDirective } from 'ng2-charts';
import { Observable } from 'rxjs';
import { ChartData } from 'chart.js';
import { ChartDataService } from '../services/chart-data.service';

@Component({
  selector: 'app-fetched-chart',
  standalone: true,
  imports: [CommonModule, AsyncPipe, BaseChartDirective],
  template: `
    @if (chartData$ | async; as chartData) {
      <canvas baseChart
        [data]="chartData"
        [options]="chartOptions"
        type="bar">
      </canvas>
    } @else {
      <div class="text-gray-500">Loading chart...</div>
    }
  `
})
export class FetchedChartComponent implements OnInit {
  private chartDataService = inject(ChartDataService);

  chartData$!: Observable<ChartData<'bar'>>;

  chartOptions = {
    responsive: true,
    plugins: {
      title: { display: true, text: 'Monthly Revenue' }
    }
  };

  ngOnInit(): void {
    this.chartData$ = this.chartDataService.getSalesData();
  }
}
```

## Reusable Chart Component

```typescript
// components/chart.component.ts
import { Component, Input, ViewChild, OnChanges, SimpleChanges } from '@angular/core';
import { BaseChartDirective } from 'ng2-charts';
import {
  Chart,
  ChartConfiguration,
  ChartData,
  ChartType,
  ChartOptions,
  registerables
} from 'chart.js';

Chart.register(...registerables);

type SupportedChartType = 'bar' | 'line' | 'pie' | 'doughnut' | 'radar' | 'polarArea';

@Component({
  selector: 'app-chart',
  standalone: true,
  imports: [BaseChartDirective],
  template: `
    <div [style.height.px]="height">
      <canvas baseChart
        [data]="data"
        [options]="mergedOptions"
        [type]="type">
      </canvas>
    </div>
  `
})
export class ChartComponent implements OnChanges {
  @ViewChild(BaseChartDirective) chart?: BaseChartDirective;

  @Input({ required: true }) type!: SupportedChartType;
  @Input({ required: true }) data!: ChartData;
  @Input() options: ChartOptions = {};
  @Input() height = 400;

  mergedOptions: ChartOptions = {};

  ngOnChanges(changes: SimpleChanges): void {
    if (changes['options'] || changes['type']) {
      this.mergedOptions = {
        responsive: true,
        maintainAspectRatio: false,
        ...this.options
      };
    }

    if (changes['data'] && !changes['data'].firstChange) {
      this.chart?.update();
    }
  }
}
```

## Chart with Signals (Angular 17+)

```typescript
// components/signal-chart.component.ts
import { Component, signal, computed, effect } from '@angular/core';
import { BaseChartDirective } from 'ng2-charts';
import { Chart, ChartData, CategoryScale, LinearScale, BarElement, Tooltip } from 'chart.js';

Chart.register(CategoryScale, LinearScale, BarElement, Tooltip);

@Component({
  selector: 'app-signal-chart',
  standalone: true,
  imports: [BaseChartDirective],
  template: `
    <canvas baseChart
      [data]="chartData()"
      [options]="chartOptions"
      type="bar">
    </canvas>
    <div class="mt-4 space-x-2">
      <button (click)="addValue()" class="px-4 py-2 bg-blue-500 text-white rounded">
        Add Value
      </button>
      <button (click)="randomize()" class="px-4 py-2 bg-green-500 text-white rounded">
        Randomize
      </button>
    </div>
    <p class="mt-2 text-sm text-gray-600">Total: {{ total() }}</p>
  `
})
export class SignalChartComponent {
  private labels = signal<string[]>(['A', 'B', 'C', 'D', 'E']);
  private values = signal<number[]>([10, 20, 30, 40, 50]);

  // Computed signal for chart data
  chartData = computed<ChartData<'bar'>>(() => ({
    labels: this.labels(),
    datasets: [{
      label: 'Values',
      data: this.values(),
      backgroundColor: 'rgba(54, 162, 235, 0.5)'
    }]
  }));

  // Computed signal for derived data
  total = computed(() =>
    this.values().reduce((sum, val) => sum + val, 0)
  );

  chartOptions = {
    responsive: true,
    plugins: {
      legend: { display: false }
    }
  };

  addValue(): void {
    this.labels.update(l => [...l, `Item ${l.length + 1}`]);
    this.values.update(v => [...v, Math.floor(Math.random() * 100)]);
  }

  randomize(): void {
    this.values.update(v => v.map(() => Math.floor(Math.random() * 100)));
  }
}
```

## Chart Event Handling

```typescript
// components/interactive-chart.component.ts
import { Component, Output, EventEmitter } from '@angular/core';
import { BaseChartDirective } from 'ng2-charts';
import { Chart, ChartData, ChartEvent, ActiveElement, CategoryScale, LinearScale, BarElement } from 'chart.js';

Chart.register(CategoryScale, LinearScale, BarElement);

interface ChartClickEvent {
  datasetIndex: number;
  index: number;
  label: string;
  value: number;
}

@Component({
  selector: 'app-interactive-chart',
  standalone: true,
  imports: [BaseChartDirective],
  template: `
    <canvas baseChart
      [data]="chartData"
      [options]="chartOptions"
      type="bar"
      (chartClick)="onChartClick($event)">
    </canvas>
    @if (selectedItem) {
      <div class="mt-4 p-4 bg-blue-100 rounded">
        Selected: {{ selectedItem.label }} - {{ selectedItem.value }}
      </div>
    }
  `
})
export class InteractiveChartComponent {
  @Output() barClick = new EventEmitter<ChartClickEvent>();

  selectedItem: ChartClickEvent | null = null;

  chartData: ChartData<'bar'> = {
    labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May'],
    datasets: [{
      label: 'Sales',
      data: [12, 19, 3, 5, 2],
      backgroundColor: 'rgba(54, 162, 235, 0.5)'
    }]
  };

  chartOptions = {
    responsive: true,
    onClick: (event: ChartEvent, elements: ActiveElement[]) => {
      // Alternative: handle click in options
    }
  };

  onChartClick(event: { event?: ChartEvent; active?: ActiveElement[] }): void {
    if (!event.active?.length) return;

    const element = event.active[0];
    const datasetIndex = element.datasetIndex;
    const index = element.index;
    const labels = this.chartData.labels as string[];
    const data = this.chartData.datasets[datasetIndex].data as number[];

    this.selectedItem = {
      datasetIndex,
      index,
      label: labels[index],
      value: data[index]
    };

    this.barClick.emit(this.selectedItem);
  }
}
```

## Chart Theme Service

```typescript
// services/chart-theme.service.ts
import { Injectable, signal, computed } from '@angular/core';
import { ChartOptions } from 'chart.js';

interface ChartTheme {
  colors: string[];
  fontFamily: string;
  gridColor: string;
  backgroundColor: string;
}

@Injectable({ providedIn: 'root' })
export class ChartThemeService {
  private theme = signal<ChartTheme>({
    colors: [
      'rgb(54, 162, 235)',
      'rgb(255, 99, 132)',
      'rgb(75, 192, 192)',
      'rgb(255, 205, 86)'
    ],
    fontFamily: 'Inter, system-ui, sans-serif',
    gridColor: 'rgba(0, 0, 0, 0.1)',
    backgroundColor: 'white'
  });

  readonly colors = computed(() => this.theme().colors);

  readonly defaultOptions = computed<ChartOptions>(() => ({
    responsive: true,
    maintainAspectRatio: true,
    plugins: {
      legend: {
        labels: {
          font: { family: this.theme().fontFamily }
        }
      }
    },
    scales: {
      x: {
        grid: { color: this.theme().gridColor },
        ticks: { font: { family: this.theme().fontFamily } }
      },
      y: {
        grid: { color: this.theme().gridColor },
        ticks: { font: { family: this.theme().fontFamily } }
      }
    }
  }));

  setTheme(theme: Partial<ChartTheme>): void {
    this.theme.update(current => ({ ...current, ...theme }));
  }

  getColor(index: number): string {
    const colors = this.colors();
    return colors[index % colors.length];
  }
}
```

## Testing Angular Charts

```typescript
// components/bar-chart.component.spec.ts
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { BarChartComponent } from './bar-chart.component';

describe('BarChartComponent', () => {
  let component: BarChartComponent;
  let fixture: ComponentFixture<BarChartComponent>;

  const mockData = {
    labels: ['A', 'B', 'C'],
    datasets: [{ label: 'Test', data: [1, 2, 3] }]
  };

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [BarChartComponent]
    }).compileComponents();

    fixture = TestBed.createComponent(BarChartComponent);
    component = fixture.componentInstance;
    component.chartData = mockData;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  it('should render canvas element', () => {
    const canvas = fixture.nativeElement.querySelector('canvas');
    expect(canvas).toBeTruthy();
  });

  it('should have baseChart directive', () => {
    const canvas = fixture.nativeElement.querySelector('canvas[baseChart]');
    expect(canvas).toBeTruthy();
  });

  it('should accept chart options', () => {
    component.chartOptions = {
      plugins: { title: { display: true, text: 'Test' } }
    };
    fixture.detectChanges();
    expect(component.chartOptions.plugins?.title?.text).toBe('Test');
  });
});
```

```typescript
// services/chart-data.service.spec.ts
import { TestBed } from '@angular/core/testing';
import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { ChartDataService } from './chart-data.service';

describe('ChartDataService', () => {
  let service: ChartDataService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [ChartDataService]
    });

    service = TestBed.inject(ChartDataService);
    httpMock = TestBed.inject(HttpTestingController);
  });

  afterEach(() => {
    httpMock.verify();
  });

  it('should fetch and transform sales data', () => {
    const mockResponse = [
      { month: 'Jan', revenue: 1000 },
      { month: 'Feb', revenue: 2000 }
    ];

    service.getSalesData().subscribe(chartData => {
      expect(chartData.labels).toEqual(['Jan', 'Feb']);
      expect(chartData.datasets[0].data).toEqual([1000, 2000]);
    });

    const req = httpMock.expectOne('/api/sales');
    expect(req.request.method).toBe('GET');
    req.flush(mockResponse);
  });
});
```
