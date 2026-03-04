# Chart.js Date Adapters

Date adapters are required for time-based axes in Chart.js. Install one based on your project's date library.

## Available Adapters

| Adapter | Date Library | Install Command |
|---------|--------------|-----------------|
| chartjs-adapter-date-fns | date-fns | `npm install chartjs-adapter-date-fns date-fns` |
| chartjs-adapter-moment | Moment.js | `npm install chartjs-adapter-moment moment` |
| chartjs-adapter-luxon | Luxon | `npm install chartjs-adapter-luxon luxon` |
| chartjs-adapter-dayjs-4 | Day.js | `npm install chartjs-adapter-dayjs-4 dayjs` |

## Recommended: date-fns Adapter

date-fns is recommended for its tree-shakable design and modern API.

### Installation

```bash
npm install chartjs-adapter-date-fns date-fns
```

### Usage

```javascript
import { Chart, TimeScale } from 'chart.js';
import 'chartjs-adapter-date-fns';

Chart.register(TimeScale);

new Chart(ctx, {
  type: 'line',
  data: {
    datasets: [{
      data: [
        { x: '2024-01-01', y: 10 },
        { x: '2024-02-01', y: 20 },
        { x: '2024-03-01', y: 15 }
      ]
    }]
  },
  options: {
    scales: {
      x: {
        type: 'time',
        time: {
          unit: 'month'
        }
      }
    }
  }
});
```

## Date Formats

### Input Formats

```javascript
// ISO 8601 strings (recommended)
data: ['2024-01-15', '2024-02-15', '2024-03-15']

// Date objects
data: [new Date(2024, 0, 15), new Date(2024, 1, 15)]

// Timestamps (milliseconds)
data: [1705276800000, 1707955200000]

// Object format
data: [
  { x: '2024-01-15', y: 10 },
  { x: new Date(2024, 1, 15), y: 20 }
]
```

### Display Formats

```javascript
options: {
  scales: {
    x: {
      type: 'time',
      time: {
        displayFormats: {
          millisecond: 'HH:mm:ss.SSS',
          second: 'HH:mm:ss',
          minute: 'HH:mm',
          hour: 'HH:mm',
          day: 'MMM d',
          week: 'PP',
          month: 'MMM yyyy',
          quarter: 'QQQ yyyy',
          year: 'yyyy'
        }
      }
    }
  }
}
```

## Time Units

Available units for the `time.unit` option:

| Unit | Description |
|------|-------------|
| `millisecond` | Millisecond precision |
| `second` | Second precision |
| `minute` | Minute precision |
| `hour` | Hourly intervals |
| `day` | Daily intervals |
| `week` | Weekly intervals |
| `month` | Monthly intervals |
| `quarter` | Quarterly intervals |
| `year` | Yearly intervals |

## Time Scale Configuration

```javascript
scales: {
  x: {
    type: 'time',
    time: {
      unit: 'day',              // Time unit to use
      stepSize: 1,              // Step between ticks
      round: 'day',             // Round timestamps to unit
      displayFormats: {         // Format per unit
        day: 'MMM d'
      },
      tooltipFormat: 'PPP',     // Tooltip date format
      parser: 'yyyy-MM-dd'      // Custom parse format
    },
    min: '2024-01-01',          // Min date (ISO string or Date)
    max: '2024-12-31',          // Max date
    ticks: {
      source: 'auto'            // auto, data, labels
    }
  }
}
```

## Time Series Scale

For evenly-spaced data points:

```javascript
scales: {
  x: {
    type: 'timeseries',         // Use timeseries instead of time
    time: {
      unit: 'day'
    }
  }
}
```

The `timeseries` scale is optimized for data where:
- Time intervals are consistent
- Data points are evenly distributed
- Performance is critical (large datasets)
