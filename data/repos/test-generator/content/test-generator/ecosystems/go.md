# Go Testing Ecosystem

Comprehensive reference for testing Go projects.

## Detection

**Manifest files:** `go.mod`, `go.sum`
**Test frameworks:** Built-in `testing` package, testify, ginkgo

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `go.mod` exists | Go (built-in testing) |
| `github.com/stretchr/testify` in go.mod | testify assertions |
| `github.com/onsi/ginkgo` in go.mod | Ginkgo BDD |
| `github.com/onsi/gomega` in go.mod | Gomega matchers |

---

## File Structure

### Naming Conventions
```
*_test.go           # REQUIRED by Go toolchain
```

### Directory Patterns

**Colocated (standard Go convention):**
```
mypackage/
├── user.go
├── user_test.go
├── service.go
├── service_test.go
└── helpers_test.go    # Test-only helpers
```

**Integration tests in separate directory:**
```
mypackage/
├── user.go
├── user_test.go       # Unit tests
integration/
├── api_test.go        # Integration tests
└── db_test.go
```

### Function Naming
```go
func TestFunctionName(t *testing.T)          // Standard test
func TestFunctionName_Scenario(t *testing.T) // Scenario-specific
func BenchmarkFunctionName(b *testing.B)     // Benchmark
func ExampleFunctionName()                    // Example (documentation)
```

---

## Standard Library Testing

### Basic Structure
```go
package mypackage

import "testing"

func TestAdd(t *testing.T) {
    // Arrange
    a, b := 2, 3
    expected := 5

    // Act
    result := Add(a, b)

    // Assert
    if result != expected {
        t.Errorf("Add(%d, %d) = %d; want %d", a, b, result, expected)
    }
}

func TestAdd_NegativeNumbers(t *testing.T) {
    result := Add(-1, -2)
    if result != -3 {
        t.Errorf("Add(-1, -2) = %d; want -3", result)
    }
}
```

### Test Methods
```go
// Report failure but continue
t.Error("message")
t.Errorf("formatted %s", "message")

// Report failure and stop test
t.Fatal("message")
t.Fatalf("formatted %s", "message")

// Log information
t.Log("info message")
t.Logf("formatted %s", "info")

// Skip test
t.Skip("reason to skip")
t.Skipf("skip because %s", "reason")

// Run subtest
t.Run("subtest name", func(t *testing.T) {
    // subtest code
})

// Parallel execution
t.Parallel()

// Cleanup (runs after test completes)
t.Cleanup(func() {
    // cleanup code
})

// Temporary directory (auto-cleaned)
dir := t.TempDir()
```

### Table-Driven Tests
```go
func TestAdd(t *testing.T) {
    tests := []struct {
        name     string
        a, b     int
        expected int
    }{
        {"positive numbers", 2, 3, 5},
        {"negative numbers", -1, -2, -3},
        {"mixed numbers", -1, 3, 2},
        {"zeros", 0, 0, 0},
    }

    for _, tt := range tests {
        t.Run(tt.name, func(t *testing.T) {
            result := Add(tt.a, tt.b)
            if result != tt.expected {
                t.Errorf("Add(%d, %d) = %d; want %d",
                    tt.a, tt.b, result, tt.expected)
            }
        })
    }
}
```

### Table-Driven with Error Cases
```go
func TestDivide(t *testing.T) {
    tests := []struct {
        name      string
        a, b      int
        expected  int
        wantErr   bool
        errString string
    }{
        {"valid division", 10, 2, 5, false, ""},
        {"division by zero", 10, 0, 0, true, "division by zero"},
    }

    for _, tt := range tests {
        t.Run(tt.name, func(t *testing.T) {
            result, err := Divide(tt.a, tt.b)

            if tt.wantErr {
                if err == nil {
                    t.Errorf("expected error, got nil")
                }
                if err.Error() != tt.errString {
                    t.Errorf("error = %v; want %v", err.Error(), tt.errString)
                }
                return
            }

            if err != nil {
                t.Errorf("unexpected error: %v", err)
            }
            if result != tt.expected {
                t.Errorf("result = %d; want %d", result, tt.expected)
            }
        })
    }
}
```

### Setup and Teardown
```go
func TestMain(m *testing.M) {
    // Setup before all tests
    setup()

    // Run tests
    code := m.Run()

    // Teardown after all tests
    teardown()

    os.Exit(code)
}

// Per-test setup using t.Cleanup
func TestWithCleanup(t *testing.T) {
    resource := createResource()
    t.Cleanup(func() {
        resource.Close()
    })

    // Test using resource
}
```

---

## Testify Patterns

### Assertions
```go
import (
    "testing"
    "github.com/stretchr/testify/assert"
)

func TestWithAssert(t *testing.T) {
    // Equality
    assert.Equal(t, expected, actual)
    assert.NotEqual(t, unexpected, actual)
    assert.EqualValues(t, expected, actual)  // Type conversion

    // Boolean
    assert.True(t, condition)
    assert.False(t, condition)
    assert.Nil(t, value)
    assert.NotNil(t, value)

    // Comparisons
    assert.Greater(t, 5, value)
    assert.GreaterOrEqual(t, 5, value)
    assert.Less(t, 10, value)

    // Strings
    assert.Contains(t, "hello world", "world")
    assert.NotContains(t, "hello", "world")
    assert.Regexp(t, `\d+`, "123")

    // Collections
    assert.Len(t, slice, 3)
    assert.Empty(t, slice)
    assert.NotEmpty(t, slice)
    assert.ElementsMatch(t, expected, actual)  // Order-independent

    // Errors
    assert.Error(t, err)
    assert.NoError(t, err)
    assert.EqualError(t, err, "expected message")
    assert.ErrorContains(t, err, "partial message")
    assert.ErrorIs(t, err, ErrNotFound)

    // Panics
    assert.Panics(t, func() { panicFunc() })
    assert.NotPanics(t, func() { safeFunc() })
    assert.PanicsWithValue(t, "panic message", func() { panicFunc() })
}
```

### Require (Fatal on Failure)
```go
import "github.com/stretchr/testify/require"

func TestWithRequire(t *testing.T) {
    result, err := SomeFunction()

    // Stops test if fails (useful for setup)
    require.NoError(t, err)
    require.NotNil(t, result)

    // Continue with assertions
    assert.Equal(t, "expected", result.Field)
}
```

### Testify Suites
```go
import (
    "testing"
    "github.com/stretchr/testify/suite"
)

type UserServiceTestSuite struct {
    suite.Suite
    service *UserService
    db      *sql.DB
}

func (s *UserServiceTestSuite) SetupSuite() {
    // Run once before all tests
    s.db = setupTestDB()
}

func (s *UserServiceTestSuite) TearDownSuite() {
    // Run once after all tests
    s.db.Close()
}

func (s *UserServiceTestSuite) SetupTest() {
    // Run before each test
    s.service = NewUserService(s.db)
    s.db.Exec("DELETE FROM users")
}

func (s *UserServiceTestSuite) TestCreateUser() {
    user, err := s.service.Create("John", "john@example.com")

    s.NoError(err)
    s.NotNil(user)
    s.Equal("John", user.Name)
}

func (s *UserServiceTestSuite) TestFindUser() {
    // Uses setup from SetupTest
    s.service.Create("Jane", "jane@example.com")

    user, err := s.service.FindByEmail("jane@example.com")

    s.NoError(err)
    s.Equal("Jane", user.Name)
}

func TestUserServiceSuite(t *testing.T) {
    suite.Run(t, new(UserServiceTestSuite))
}
```

---

## Mocking

### Testify Mock
```go
import (
    "github.com/stretchr/testify/mock"
)

// Define mock
type MockRepository struct {
    mock.Mock
}

func (m *MockRepository) FindByID(id int) (*User, error) {
    args := m.Called(id)
    if args.Get(0) == nil {
        return nil, args.Error(1)
    }
    return args.Get(0).(*User), args.Error(1)
}

func (m *MockRepository) Save(user *User) error {
    args := m.Called(user)
    return args.Error(0)
}

// Use in tests
func TestUserService_GetUser(t *testing.T) {
    mockRepo := new(MockRepository)
    service := NewUserService(mockRepo)

    expectedUser := &User{ID: 1, Name: "John"}
    mockRepo.On("FindByID", 1).Return(expectedUser, nil)

    user, err := service.GetUser(1)

    assert.NoError(t, err)
    assert.Equal(t, "John", user.Name)
    mockRepo.AssertExpectations(t)
}

func TestUserService_GetUser_NotFound(t *testing.T) {
    mockRepo := new(MockRepository)
    service := NewUserService(mockRepo)

    mockRepo.On("FindByID", 999).Return(nil, ErrNotFound)

    user, err := service.GetUser(999)

    assert.Nil(t, user)
    assert.ErrorIs(t, err, ErrNotFound)
    mockRepo.AssertCalled(t, "FindByID", 999)
}
```

### Interface-Based Mocking
```go
// Define interface
type UserRepository interface {
    FindByID(id int) (*User, error)
    Save(user *User) error
}

// Production implementation
type PostgresUserRepository struct {
    db *sql.DB
}

// Test implementation
type InMemoryUserRepository struct {
    users map[int]*User
}

func (r *InMemoryUserRepository) FindByID(id int) (*User, error) {
    user, ok := r.users[id]
    if !ok {
        return nil, ErrNotFound
    }
    return user, nil
}

// Test uses in-memory implementation
func TestUserService(t *testing.T) {
    repo := &InMemoryUserRepository{
        users: map[int]*User{
            1: {ID: 1, Name: "John"},
        },
    }
    service := NewUserService(repo)

    user, err := service.GetUser(1)

    assert.NoError(t, err)
    assert.Equal(t, "John", user.Name)
}
```

---

## HTTP Testing

### Handler Testing
```go
import (
    "net/http"
    "net/http/httptest"
    "testing"
)

func TestGetUserHandler(t *testing.T) {
    // Create request
    req := httptest.NewRequest(http.MethodGet, "/users/1", nil)
    rec := httptest.NewRecorder()

    // Call handler
    handler := NewUserHandler(mockService)
    handler.GetUser(rec, req)

    // Assert response
    assert.Equal(t, http.StatusOK, rec.Code)
    assert.Contains(t, rec.Body.String(), "John")
}

func TestPostUserHandler(t *testing.T) {
    body := strings.NewReader(`{"name":"John","email":"john@example.com"}`)
    req := httptest.NewRequest(http.MethodPost, "/users", body)
    req.Header.Set("Content-Type", "application/json")
    rec := httptest.NewRecorder()

    handler.CreateUser(rec, req)

    assert.Equal(t, http.StatusCreated, rec.Code)
}
```

### Server Testing
```go
func TestAPIEndpoint(t *testing.T) {
    // Create test server
    server := httptest.NewServer(NewRouter())
    defer server.Close()

    // Make request
    resp, err := http.Get(server.URL + "/api/users")
    require.NoError(t, err)
    defer resp.Body.Close()

    // Assert
    assert.Equal(t, http.StatusOK, resp.StatusCode)

    var users []User
    json.NewDecoder(resp.Body).Decode(&users)
    assert.Len(t, users, 3)
}
```

---

## Benchmarks

```go
func BenchmarkFunction(b *testing.B) {
    // Setup outside the loop
    data := prepareTestData()

    b.ResetTimer() // Reset timer after setup

    for i := 0; i < b.N; i++ {
        FunctionToTest(data)
    }
}

func BenchmarkParallel(b *testing.B) {
    b.RunParallel(func(pb *testing.PB) {
        for pb.Next() {
            FunctionToTest()
        }
    })
}

// Run: go test -bench=. -benchmem
```

---

## Common Patterns

### Testing Errors
```go
func TestFunction_ReturnsError(t *testing.T) {
    _, err := Function(invalidInput)

    // Check error exists
    assert.Error(t, err)

    // Check specific error
    assert.ErrorIs(t, err, ErrInvalidInput)

    // Check error type
    var validationErr *ValidationError
    assert.ErrorAs(t, err, &validationErr)
    assert.Equal(t, "field", validationErr.Field)
}
```

### Testing Context Cancellation
```go
func TestFunction_RespectsContext(t *testing.T) {
    ctx, cancel := context.WithCancel(context.Background())
    cancel() // Cancel immediately

    _, err := FunctionWithContext(ctx)

    assert.ErrorIs(t, err, context.Canceled)
}
```

### Testing Concurrent Code
```go
func TestConcurrentAccess(t *testing.T) {
    counter := NewSafeCounter()

    var wg sync.WaitGroup
    for i := 0; i < 100; i++ {
        wg.Add(1)
        go func() {
            defer wg.Done()
            counter.Increment()
        }()
    }
    wg.Wait()

    assert.Equal(t, 100, counter.Value())
}
```

---

## Configuration

### Run Tests
```bash
# Run all tests
go test ./...

# Verbose output
go test -v ./...

# Run specific test
go test -run TestFunctionName ./...

# With coverage
go test -cover ./...
go test -coverprofile=coverage.out ./...
go tool cover -html=coverage.out

# Race detection
go test -race ./...
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
