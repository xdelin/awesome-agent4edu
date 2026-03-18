# Chart.js with Rails 8 + Stimulus Examples

Complete examples for integrating Chart.js with Rails 8 using Hotwire.

## Basic Stimulus Controller

```javascript
// app/javascript/controllers/chart_controller.js
import { Controller } from "@hotwired/stimulus";
import Chart from "chart.js/auto";

export default class extends Controller {
  static targets = ["canvas"];
  static values = {
    type: { type: String, default: "bar" },
    data: Object,
    options: { type: Object, default: {} }
  };

  connect() {
    this.initChart();
  }

  disconnect() {
    this.destroyChart();
  }

  initChart() {
    if (!this.hasCanvasTarget) return;

    const ctx = this.canvasTarget.getContext("2d");
    this.chart = new Chart(ctx, {
      type: this.typeValue,
      data: this.dataValue,
      options: {
        responsive: true,
        maintainAspectRatio: true,
        ...this.optionsValue
      }
    });
  }

  destroyChart() {
    if (this.chart) {
      this.chart.destroy();
      this.chart = null;
    }
  }

  // Automatically called when data-chart-data-value changes
  dataValueChanged() {
    if (this.chart) {
      this.chart.data = this.dataValue;
      this.chart.update();
    }
  }

  // Refresh chart (useful for Turbo)
  refresh() {
    this.destroyChart();
    this.initChart();
  }
}
```

## View Helper

```ruby
# app/helpers/charts_helper.rb
module ChartsHelper
  def chart_tag(type:, data:, options: {}, height: nil, **html_options)
    wrapper_style = height ? "height: #{height}px" : nil

    content_tag :div,
      style: wrapper_style,
      data: {
        controller: "chart",
        chart_type_value: type,
        chart_data_value: data.to_json,
        chart_options_value: options.to_json
      },
      **html_options do
        content_tag :canvas, "", data: { chart_target: "canvas" }
      end
  end

  def bar_chart(data:, **options)
    chart_tag(type: "bar", data: data, **options)
  end

  def line_chart(data:, **options)
    chart_tag(type: "line", data: data, **options)
  end

  def pie_chart(data:, **options)
    chart_tag(type: "pie", data: data, **options)
  end

  def doughnut_chart(data:, **options)
    chart_tag(type: "doughnut", data: data, **options)
  end
end
```

## Controller Example

```ruby
# app/controllers/dashboard_controller.rb
class DashboardController < ApplicationController
  def index
    @sales_data = build_sales_chart_data
    @revenue_data = build_revenue_chart_data
  end

  private

  def build_sales_chart_data
    {
      labels: Date::MONTHNAMES.compact.first(6),
      datasets: [{
        label: "Sales 2024",
        data: Sale.monthly_totals(2024),
        backgroundColor: "rgba(54, 162, 235, 0.5)",
        borderColor: "rgb(54, 162, 235)",
        borderWidth: 1
      }, {
        label: "Sales 2023",
        data: Sale.monthly_totals(2023),
        backgroundColor: "rgba(255, 99, 132, 0.5)",
        borderColor: "rgb(255, 99, 132)",
        borderWidth: 1
      }]
    }
  end

  def build_revenue_chart_data
    {
      labels: %w[Products Services Subscriptions],
      datasets: [{
        data: Revenue.by_category,
        backgroundColor: [
          "rgb(255, 99, 132)",
          "rgb(54, 162, 235)",
          "rgb(255, 205, 86)"
        ]
      }]
    }
  end
end
```

## View Example

```erb
<!-- app/views/dashboard/index.html.erb -->
<div class="grid grid-cols-2 gap-6">
  <div class="bg-white p-6 rounded-lg shadow">
    <h2 class="text-lg font-semibold mb-4">Monthly Sales</h2>
    <%= bar_chart data: @sales_data,
                  options: { plugins: { legend: { position: "bottom" } } },
                  height: 300 %>
  </div>

  <div class="bg-white p-6 rounded-lg shadow">
    <h2 class="text-lg font-semibold mb-4">Revenue by Category</h2>
    <%= doughnut_chart data: @revenue_data, height: 300 %>
  </div>
</div>
```

## Turbo Stream Updates

### Controller

```ruby
# app/controllers/sales_controller.rb
class SalesController < ApplicationController
  def create
    @sale = Sale.create!(sale_params)

    respond_to do |format|
      format.turbo_stream
      format.html { redirect_to dashboard_path }
    end
  end

  private

  def sale_params
    params.require(:sale).permit(:amount, :category)
  end
end
```

### Turbo Stream View

```erb
<!-- app/views/sales/create.turbo_stream.erb -->
<%= turbo_stream.replace "sales-chart" do %>
  <%= render "dashboard/sales_chart", data: build_updated_chart_data %>
<% end %>

<%= turbo_stream.prepend "notifications" do %>
  <div class="alert alert-success">Sale recorded!</div>
<% end %>
```

### Partial

```erb
<!-- app/views/dashboard/_sales_chart.html.erb -->
<div id="sales-chart">
  <%= bar_chart data: data,
                options: { animation: { duration: 500 } } %>
</div>
```

## Live Updates with ActionCable

### Channel

```ruby
# app/channels/charts_channel.rb
class ChartsChannel < ApplicationCable::Channel
  def subscribed
    stream_from "charts_#{params[:chart_id]}"
  end
end
```

### Stimulus Controller with Cable

```javascript
// app/javascript/controllers/live_chart_controller.js
import { Controller } from "@hotwired/stimulus";
import { createConsumer } from "@rails/actioncable";
import Chart from "chart.js/auto";

export default class extends Controller {
  static targets = ["canvas"];
  static values = {
    chartId: String,
    config: Object
  };

  connect() {
    this.initChart();
    this.subscribeToChannel();
  }

  disconnect() {
    this.unsubscribe();
    this.destroyChart();
  }

  initChart() {
    const ctx = this.canvasTarget.getContext("2d");
    this.chart = new Chart(ctx, this.configValue);
  }

  destroyChart() {
    if (this.chart) {
      this.chart.destroy();
      this.chart = null;
    }
  }

  subscribeToChannel() {
    this.subscription = createConsumer().subscriptions.create(
      { channel: "ChartsChannel", chart_id: this.chartIdValue },
      {
        received: (data) => this.handleUpdate(data)
      }
    );
  }

  unsubscribe() {
    if (this.subscription) {
      this.subscription.unsubscribe();
    }
  }

  handleUpdate(data) {
    if (data.action === "add_point") {
      this.addDataPoint(data.label, data.value);
    } else if (data.action === "update") {
      this.chart.data = data.chartData;
      this.chart.update();
    }
  }

  addDataPoint(label, value) {
    this.chart.data.labels.push(label);
    this.chart.data.datasets[0].data.push(value);

    // Keep last 20 points
    if (this.chart.data.labels.length > 20) {
      this.chart.data.labels.shift();
      this.chart.data.datasets[0].data.shift();
    }

    this.chart.update();
  }
}
```

### Broadcasting Updates

```ruby
# app/models/sale.rb
class Sale < ApplicationRecord
  after_create_commit :broadcast_chart_update

  private

  def broadcast_chart_update
    ActionCable.server.broadcast(
      "charts_sales",
      {
        action: "add_point",
        label: created_at.strftime("%H:%M"),
        value: amount
      }
    )
  end
end
```

## Importmap Configuration

```ruby
# config/importmap.rb
pin "application", preload: true
pin "@hotwired/turbo-rails", to: "turbo.min.js", preload: true
pin "@hotwired/stimulus", to: "stimulus.min.js", preload: true
pin "@hotwired/stimulus-loading", to: "stimulus-loading.js", preload: true

# Chart.js
pin "chart.js", to: "https://cdn.jsdelivr.net/npm/chart.js@4.5.1/dist/chart.umd.min.js"

# Or use auto build
pin "chart.js/auto", to: "https://cdn.jsdelivr.net/npm/chart.js@4.5.1/auto/+esm"

pin_all_from "app/javascript/controllers", under: "controllers"
```

## Testing Charts

```ruby
# test/system/dashboard_test.rb
require "application_system_test_case"

class DashboardTest < ApplicationSystemTestCase
  test "displays sales chart" do
    visit dashboard_path

    assert_selector "[data-controller='chart']"
    assert_selector "canvas[data-chart-target='canvas']"
  end

  test "updates chart when new sale is created" do
    visit dashboard_path

    # Create sale via form
    fill_in "Amount", with: 100
    click_button "Add Sale"

    # Chart should update
    assert_selector "[data-controller='chart']"
  end
end
```
