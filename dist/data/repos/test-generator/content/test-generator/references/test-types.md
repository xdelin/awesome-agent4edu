# Test Types Reference

Comprehensive guide for selecting the right test type based on code characteristics.

## The Testing Trophy

Prioritize tests that give maximum confidence with minimum cost:

```
        ╱╲
       ╱  ╲         E2E Tests (few, critical paths only)
      ╱────╲
     ╱      ╲
    ╱        ╲      Integration Tests (most value)
   ╱──────────╲
  ╱            ╲
 ╱              ╲   Unit Tests (fast, focused)
╱────────────────╲
══════════════════  Static Analysis (linting, types)
```

**Principle:** The more tests resemble real usage, the more confidence they provide.

---

## Unit Tests

### When to Use
- Pure functions with clear inputs/outputs
- Utility functions and helpers
- Complex algorithms with many edge cases
- State transformations
- Parsers and formatters
- Validation logic

### Indicators in Code
```
✓ No external dependencies (database, API, filesystem)
✓ Deterministic output for given input
✓ High cyclomatic complexity (many branches)
✓ Mathematical or string operations
✓ Data transformation logic
```

### Structure: AAA Pattern
```
Arrange  → Set up test data and conditions
Act      → Execute the function under test
Assert   → Verify the expected outcome
```

### What to Test
- Happy path (expected inputs)
- Edge cases (empty, null, undefined, boundary values)
- Error cases (invalid inputs, exceptions)
- Type coercion edge cases

### What NOT to Test
- Implementation details (private methods)
- Framework code
- Third-party library behavior

---

## Integration Tests

### When to Use
- Service layers with multiple dependencies
- Database operations (repositories, DAOs)
- API client wrappers
- Code crossing module boundaries
- Message queue handlers
- Cache interactions

### Indicators in Code
```
✓ Class/function depends on other services
✓ Uses ORM or database client
✓ Makes HTTP calls to external services
✓ Reads/writes files
✓ Interacts with message brokers
```

### Mocking Strategy
- Mock at boundaries (HTTP, database, filesystem)
- Use test containers for database tests when possible
- Mock external APIs with tools like MSW, nock, or responses
- Don't mock what you own (test it integrated)

### What to Test
- Correct interaction between components
- Data flows through the system
- Error propagation across boundaries
- Transaction behavior
- Connection handling

---

## API/Endpoint Tests

### When to Use
- REST API endpoints
- GraphQL resolvers
- RPC handlers
- Webhooks
- Any public-facing HTTP interface

### Indicators in Code
```
✓ Route handlers or controllers
✓ Request/response processing
✓ Authentication/authorization checks
✓ Input validation at API boundary
✓ Response formatting
```

### What to Test

**Status Codes:**
- 200/201 for success
- 400 for validation errors
- 401 for unauthorized
- 403 for forbidden
- 404 for not found
- 500 for server errors

**Response Shape:**
- Correct JSON structure
- Required fields present
- Proper data types
- Pagination metadata

**Authentication:**
- Rejects missing tokens
- Rejects invalid tokens
- Accepts valid tokens
- Proper scope/permission checks

**Validation:**
- Rejects missing required fields
- Rejects invalid formats
- Accepts valid input

### Contract Testing
For APIs consumed by other services:
- Define expected request/response shapes
- Version contracts appropriately
- Test backward compatibility

---

## Validation Tests

### When to Use
- Form validation
- Input sanitization
- Schema validation
- Business rule validation
- Configuration validation

### Indicators in Code
```
✓ Validation functions or decorators
✓ Schema definitions (Zod, Yup, Joi, etc.)
✓ Form handling logic
✓ Input constraints (min/max, patterns)
```

### Boundary Value Analysis

For any constrained input, test at boundaries:

```
Example: Age field (18-65)

Test cases:
├── Below minimum: 17 (invalid)
├── At minimum: 18 (valid)
├── Normal value: 30 (valid)
├── At maximum: 65 (valid)
├── Above maximum: 66 (invalid)
├── Edge: 0 (invalid)
├── Edge: negative (-1) (invalid)
├── Edge: very large (999999) (invalid)
└── Edge: non-integer (18.5) (depends on rules)
```

### Equivalence Partitioning

Group inputs into classes that should behave the same:

```
Email validation:
├── Valid emails (any valid format)
├── Invalid format (missing @)
├── Invalid domain (no TLD)
├── Empty string
├── Null/undefined
└── Non-string types
```

### What to Test
- Required field enforcement
- Format validation (email, phone, URL)
- Length constraints (min/max)
- Numeric constraints (range, precision)
- Pattern matching (regex)
- Cross-field validation
- Custom business rules

---

## Component Tests (Frontend)

### When to Use
- React/Vue/Angular/Svelte components
- UI widgets
- Interactive elements
- Forms

### Indicators in Code
```
✓ Component functions/classes
✓ JSX/template code
✓ Event handlers
✓ State management
✓ Props/inputs handling
```

### Testing Library Philosophy
Test behavior, not implementation:
- Find elements by accessible roles/labels
- Simulate user interactions
- Assert on visible outcomes
- Avoid testing internal state

### What to Test
- Renders without crashing
- Displays correct content
- Responds to user interactions
- Shows loading/error states
- Handles props correctly
- Accessibility (a11y)

### What NOT to Test
- CSS styling (use visual regression)
- Internal component state
- Implementation details
- Third-party component behavior

---

## E2E Tests

### When to Use (SPARINGLY)
- Critical user journeys only
- Checkout/payment flows
- Authentication flows
- Cross-page workflows
- Smoke tests for deployment

### Indicators
```
✓ Multi-step user flow
✓ Involves multiple pages/screens
✓ Business-critical path
✓ Requires real browser behavior
```

### Page Object Model
Encapsulate page interactions:

```
LoginPage
├── navigate()
├── fillUsername(value)
├── fillPassword(value)
├── submit()
└── getErrorMessage()

Test uses:
loginPage.navigate()
loginPage.fillUsername('user')
loginPage.fillPassword('pass')
loginPage.submit()
expect(dashboardPage.isVisible()).toBe(true)
```

### Best Practices
- Test happy paths primarily
- Keep E2E suite small (<50 tests)
- Run on CI, not every commit
- Use stable selectors (data-testid)
- Handle flakiness with retries
- Isolate test data

### What NOT to E2E Test
- Every edge case (use unit tests)
- Visual styling (use visual regression)
- Error messages (use integration tests)
- Input validation (use validation tests)

---

## Snapshot Tests

### When to Use (Carefully)
- UI component output
- Serializable data structures
- API response shapes
- Configuration objects

### Caution
Snapshots can become:
- Noisy (frequent updates)
- Meaningless (auto-accepted)
- Large (hard to review)

### Best Practices
- Keep snapshots small and focused
- Review snapshot changes carefully
- Don't snapshot dynamic content (timestamps, IDs)
- Consider inline snapshots for small outputs

---

## Test Selection Decision Tree

```
Is it a pure function with no dependencies?
├── YES → Unit Test
└── NO
    │
    Does it cross system boundaries (DB, API, filesystem)?
    ├── YES → Integration Test
    └── NO
        │
        Is it an HTTP endpoint?
        ├── YES → API Test
        └── NO
            │
            Is it input validation?
            ├── YES → Validation Test (boundary values)
            └── NO
                │
                Is it a UI component?
                ├── YES → Component Test
                └── NO
                    │
                    Is it a critical multi-page user flow?
                    ├── YES → E2E Test (keep minimal)
                    └── NO → Probably Unit or Integration Test
```

---

## Anti-Patterns to Avoid

### 1. Testing Implementation
❌ Testing private methods
❌ Asserting on internal state
❌ Testing that specific functions were called (without behavior check)

### 2. Over-Mocking
❌ Mocking everything including code you own
❌ Mocks that return exactly what tests expect
❌ Tests that pass when code is broken

### 3. Meaningless Assertions
❌ `expect(true).toBe(true)`
❌ `expect(result).toBeTruthy()` (when you know the expected value)
❌ Only checking that no error was thrown

### 4. Test Duplication
❌ Same logic tested in unit AND integration AND E2E
❌ Multiple tests checking the same behavior

### 5. Flaky Tests
❌ Tests depending on timing
❌ Tests depending on external services
❌ Tests with race conditions
