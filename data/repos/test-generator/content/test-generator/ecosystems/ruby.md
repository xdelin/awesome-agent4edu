# Ruby Testing Ecosystem

Comprehensive reference for testing Ruby projects.

## Detection

**Manifest files:** `Gemfile`, `*.gemspec`
**Test frameworks:** RSpec, Minitest

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `rspec` in Gemfile | RSpec |
| `spec/` directory | RSpec |
| `.rspec` file | RSpec |
| `minitest` in Gemfile | Minitest |
| `test/` directory with `*_test.rb` | Minitest |
| `rails` in Gemfile | Rails (RSpec or Minitest) |

---

## File Structure

### RSpec Naming Conventions
```
*_spec.rb           # Standard RSpec
```

### Minitest Naming Conventions
```
*_test.rb           # Standard Minitest
test_*.rb           # Alternative
```

### Directory Patterns

**RSpec structure:**
```
spec/
├── spec_helper.rb
├── rails_helper.rb          # Rails projects
├── models/
│   └── user_spec.rb
├── controllers/
│   └── users_controller_spec.rb
├── services/
│   └── payment_service_spec.rb
├── requests/                 # API/integration specs
│   └── users_spec.rb
├── features/                 # Feature/acceptance specs
│   └── user_registration_spec.rb
└── support/
    ├── factory_bot.rb
    └── helpers.rb
```

**Minitest structure:**
```
test/
├── test_helper.rb
├── models/
│   └── user_test.rb
├── controllers/
│   └── users_controller_test.rb
├── integration/
│   └── user_flows_test.rb
└── fixtures/
    └── users.yml
```

---

## RSpec Patterns

### Basic Structure
```ruby
require 'spec_helper'

RSpec.describe UserService do
  describe '#find_by_id' do
    context 'when user exists' do
      let(:user) { create(:user) }

      it 'returns the user' do
        result = described_class.find_by_id(user.id)
        expect(result).to eq(user)
      end
    end

    context 'when user does not exist' do
      it 'raises an error' do
        expect { described_class.find_by_id(-1) }
          .to raise_error(UserNotFoundError)
      end
    end
  end
end
```

### Matchers
```ruby
# Equality
expect(actual).to eq(expected)
expect(actual).to eql(expected)      # Type-strict
expect(actual).to equal(expected)    # Object identity
expect(actual).to be(expected)       # Object identity

# Boolean
expect(value).to be true
expect(value).to be false
expect(value).to be_truthy
expect(value).to be_falsy
expect(value).to be_nil

# Comparisons
expect(value).to be > 5
expect(value).to be >= 5
expect(value).to be < 10
expect(value).to be_between(1, 10).inclusive
expect(value).to be_within(0.1).of(expected)

# Strings
expect(string).to include('substring')
expect(string).to start_with('prefix')
expect(string).to end_with('suffix')
expect(string).to match(/regex/)
expect(string).to be_empty

# Collections
expect(array).to include(item)
expect(array).to include(1, 2, 3)
expect(array).to contain_exactly(1, 2, 3)  # Order-independent
expect(array).to match_array([1, 2, 3])
expect(array).to have_attributes(size: 3)
expect(array).to be_empty
expect(array).to all(be > 0)

# Hashes
expect(hash).to include(key: 'value')
expect(hash).to have_key(:key)
expect(hash).to have_value('value')

# Types
expect(object).to be_a(ClassName)
expect(object).to be_an(Array)
expect(object).to be_an_instance_of(ExactClass)
expect(object).to respond_to(:method_name)

# Predicate matchers (be_*)
expect(user).to be_active        # calls user.active?
expect(list).to be_empty         # calls list.empty?
expect(value).to be_present      # Rails: calls value.present?

# Changes
expect { action }.to change { object.value }.from(1).to(2)
expect { action }.to change { object.value }.by(1)
expect { action }.to change { Model.count }.by(1)
expect { action }.not_to change { object.value }

# Errors
expect { action }.to raise_error
expect { action }.to raise_error(ErrorClass)
expect { action }.to raise_error(ErrorClass, 'message')
expect { action }.to raise_error(/regex/)
expect { action }.not_to raise_error

# Output
expect { puts 'hello' }.to output('hello').to_stdout
expect { warn 'error' }.to output(/error/).to_stderr

# Compound matchers
expect(value).to be > 5 and be < 10
expect(string).to start_with('a').or start_with('b')
```

### Let and Let!
```ruby
RSpec.describe User do
  # Lazy evaluation - created when first accessed
  let(:user) { create(:user, name: 'John') }

  # Immediate evaluation - created before each example
  let!(:admin) { create(:user, :admin) }

  # Subject
  subject(:service) { described_class.new(user) }

  it 'has a name' do
    expect(user.name).to eq('John')
  end
end
```

### Hooks
```ruby
RSpec.describe User do
  before(:all) do
    # Run once before all examples in this group
  end

  after(:all) do
    # Run once after all examples
  end

  before(:each) do
    # Run before each example (default)
  end

  after(:each) do
    # Run after each example
  end

  around(:each) do |example|
    # Wrap each example
    DatabaseCleaner.cleaning do
      example.run
    end
  end
end
```

### Shared Examples
```ruby
# Define shared examples
RSpec.shared_examples 'a valid model' do
  it { is_expected.to be_valid }
  it { is_expected.to respond_to(:save) }
end

RSpec.shared_examples 'requires authentication' do |method, path|
  it 'returns 401 when not authenticated' do
    send(method, path)
    expect(response).to have_http_status(:unauthorized)
  end
end

# Use shared examples
RSpec.describe User do
  subject { build(:user) }
  it_behaves_like 'a valid model'
end

RSpec.describe 'Users API' do
  it_behaves_like 'requires authentication', :get, '/api/users'
end
```

### Mocking with RSpec Mocks
```ruby
RSpec.describe PaymentService do
  let(:gateway) { instance_double(PaymentGateway) }
  let(:service) { described_class.new(gateway) }

  describe '#process' do
    it 'charges the payment' do
      # Stub
      allow(gateway).to receive(:charge).and_return(true)

      result = service.process(100)

      expect(result).to be true
      # Verify
      expect(gateway).to have_received(:charge).with(100)
    end

    it 'handles errors' do
      allow(gateway).to receive(:charge)
        .and_raise(PaymentError, 'declined')

      expect { service.process(100) }
        .to raise_error(PaymentError, 'declined')
    end
  end
end

# Spies
RSpec.describe NotificationService do
  let(:mailer) { spy('mailer') }

  it 'sends notification' do
    service = described_class.new(mailer)
    service.notify(user)

    expect(mailer).to have_received(:send).with(user.email)
  end
end

# Message expectations
RSpec.describe Logger do
  it 'logs messages' do
    logger = instance_double(Logger)
    expect(logger).to receive(:info).with('message').once

    service = Service.new(logger)
    service.perform
  end
end
```

---

## Minitest Patterns

### Basic Structure
```ruby
require 'test_helper'

class UserServiceTest < Minitest::Test
  def setup
    @service = UserService.new
  end

  def teardown
    # Cleanup
  end

  def test_find_by_id_returns_user
    user = create(:user)
    
    result = @service.find_by_id(user.id)
    
    assert_equal user, result
  end

  def test_find_by_id_raises_when_not_found
    assert_raises(UserNotFoundError) do
      @service.find_by_id(-1)
    end
  end
end
```

### Assertions
```ruby
# Equality
assert_equal expected, actual
refute_equal unexpected, actual
assert_same expected, actual          # Object identity

# Boolean
assert value
refute value
assert_nil value
refute_nil value

# Comparisons
assert_operator value, :>, 5
assert_in_delta expected, actual, 0.001
assert_in_epsilon expected, actual, 0.001

# Strings/Collections
assert_includes collection, item
refute_includes collection, item
assert_match /regex/, string
assert_empty collection
refute_empty collection

# Types
assert_instance_of ExactClass, object
assert_kind_of ParentClass, object
assert_respond_to object, :method_name

# Exceptions
assert_raises(ErrorClass) { action }
error = assert_raises(ErrorClass) { action }
assert_equal 'message', error.message

# Output
assert_output('expected') { puts 'expected' }
assert_silent { quiet_method }
```

### Minitest::Spec (BDD Style)
```ruby
require 'minitest/autorun'

describe UserService do
  before do
    @service = UserService.new
  end

  describe '#find_by_id' do
    it 'returns the user when found' do
      user = create(:user)
      result = @service.find_by_id(user.id)
      _(result).must_equal user
    end

    it 'raises when not found' do
      _ { @service.find_by_id(-1) }.must_raise UserNotFoundError
    end
  end
end
```

### Spec Expectations
```ruby
_(value).must_equal expected
_(value).wont_equal unexpected
_(value).must_be :>, 5
_(value).must_be_nil
_(value).wont_be_nil
_(value).must_include item
_(value).must_match /regex/
_(value).must_be_instance_of Class
_(value).must_respond_to :method
_ { action }.must_raise ErrorClass
_ { action }.must_output 'text'
```

---

## Rails Testing

### Model Tests (RSpec)
```ruby
require 'rails_helper'

RSpec.describe User, type: :model do
  describe 'validations' do
    it { is_expected.to validate_presence_of(:email) }
    it { is_expected.to validate_uniqueness_of(:email) }
    it { is_expected.to validate_length_of(:name).is_at_most(100) }
  end

  describe 'associations' do
    it { is_expected.to have_many(:posts).dependent(:destroy) }
    it { is_expected.to belong_to(:organization).optional }
  end

  describe 'scopes' do
    describe '.active' do
      it 'returns only active users' do
        active = create(:user, active: true)
        inactive = create(:user, active: false)

        expect(User.active).to include(active)
        expect(User.active).not_to include(inactive)
      end
    end
  end

  describe '#full_name' do
    it 'combines first and last name' do
      user = build(:user, first_name: 'John', last_name: 'Doe')
      expect(user.full_name).to eq('John Doe')
    end
  end
end
```

### Request/API Tests (RSpec)
```ruby
require 'rails_helper'

RSpec.describe 'Users API', type: :request do
  describe 'GET /api/users' do
    it 'returns all users' do
      create_list(:user, 3)

      get '/api/users'

      expect(response).to have_http_status(:ok)
      expect(json_response.size).to eq(3)
    end

    it 'requires authentication' do
      get '/api/users'
      expect(response).to have_http_status(:unauthorized)
    end

    context 'when authenticated' do
      let(:user) { create(:user) }
      let(:headers) { auth_headers(user) }

      it 'returns users' do
        get '/api/users', headers: headers
        expect(response).to have_http_status(:ok)
      end
    end
  end

  describe 'POST /api/users' do
    let(:valid_params) do
      { user: { name: 'John', email: 'john@example.com' } }
    end

    it 'creates a new user' do
      expect {
        post '/api/users', params: valid_params, headers: headers
      }.to change(User, :count).by(1)

      expect(response).to have_http_status(:created)
    end

    it 'returns errors for invalid data' do
      post '/api/users', params: { user: { name: '' } }, headers: headers

      expect(response).to have_http_status(:unprocessable_entity)
      expect(json_response['errors']).to include('name')
    end
  end

  private

  def json_response
    JSON.parse(response.body)
  end
end
```

### Controller Tests (RSpec)
```ruby
require 'rails_helper'

RSpec.describe UsersController, type: :controller do
  describe 'GET #index' do
    it 'returns a success response' do
      get :index
      expect(response).to be_successful
    end

    it 'assigns @users' do
      user = create(:user)
      get :index
      expect(assigns(:users)).to include(user)
    end
  end

  describe 'POST #create' do
    context 'with valid params' do
      let(:valid_attributes) { { name: 'John', email: 'john@example.com' } }

      it 'creates a new User' do
        expect {
          post :create, params: { user: valid_attributes }
        }.to change(User, :count).by(1)
      end

      it 'redirects to the created user' do
        post :create, params: { user: valid_attributes }
        expect(response).to redirect_to(User.last)
      end
    end

    context 'with invalid params' do
      let(:invalid_attributes) { { name: '' } }

      it 'does not create a User' do
        expect {
          post :create, params: { user: invalid_attributes }
        }.not_to change(User, :count)
      end

      it 'renders the new template' do
        post :create, params: { user: invalid_attributes }
        expect(response).to render_template(:new)
      end
    end
  end
end
```

### Feature/System Tests (Capybara)
```ruby
require 'rails_helper'

RSpec.describe 'User Registration', type: :system do
  before do
    driven_by(:selenium_chrome_headless)
  end

  it 'allows users to register' do
    visit new_user_registration_path

    fill_in 'Email', with: 'test@example.com'
    fill_in 'Password', with: 'password123'
    fill_in 'Password confirmation', with: 'password123'
    click_button 'Sign up'

    expect(page).to have_content('Welcome!')
    expect(page).to have_current_path(dashboard_path)
  end

  it 'shows validation errors' do
    visit new_user_registration_path

    fill_in 'Email', with: 'invalid'
    click_button 'Sign up'

    expect(page).to have_content('Email is invalid')
  end
end
```

---

## FactoryBot

```ruby
# spec/factories/users.rb
FactoryBot.define do
  factory :user do
    sequence(:email) { |n| "user#{n}@example.com" }
    name { 'John Doe' }
    password { 'password123' }
    active { true }

    trait :admin do
      role { 'admin' }
    end

    trait :inactive do
      active { false }
    end

    trait :with_posts do
      transient do
        posts_count { 3 }
      end

      after(:create) do |user, evaluator|
        create_list(:post, evaluator.posts_count, user: user)
      end
    end

    factory :admin_user, traits: [:admin]
  end
end

# Usage
user = create(:user)
admin = create(:user, :admin)
user_with_posts = create(:user, :with_posts, posts_count: 5)
users = create_list(:user, 3)
built_user = build(:user)  # Not persisted
attributes = attributes_for(:user)
```

---

## Configuration

### RSpec (.rspec)
```
--require spec_helper
--format documentation
--color
```

### spec_helper.rb
```ruby
RSpec.configure do |config|
  config.expect_with :rspec do |expectations|
    expectations.include_chain_clauses_in_custom_matcher_descriptions = true
  end

  config.mock_with :rspec do |mocks|
    mocks.verify_partial_doubles = true
  end

  config.shared_context_metadata_behavior = :apply_to_host_groups
  config.filter_run_when_matching :focus
  config.example_status_persistence_file_path = 'spec/examples.txt'
  config.disable_monkey_patching!
  config.order = :random
  Kernel.srand config.seed
end
```

### rails_helper.rb
```ruby
require 'spec_helper'
ENV['RAILS_ENV'] ||= 'test'
require_relative '../config/environment'

abort('Running in production!') if Rails.env.production?

require 'rspec/rails'
require 'capybara/rspec'

Dir[Rails.root.join('spec/support/**/*.rb')].each { |f| require f }

RSpec.configure do |config|
  config.fixture_path = Rails.root.join('spec/fixtures')
  config.use_transactional_fixtures = true
  config.infer_spec_type_from_file_location!
  config.filter_rails_from_backtrace!
  
  config.include FactoryBot::Syntax::Methods
  config.include Devise::Test::IntegrationHelpers, type: :request
end
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
