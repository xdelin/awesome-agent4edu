# TypeScript Testing Ecosystem

Comprehensive reference for testing TypeScript projects.

## Detection

**Manifest files:** `package.json` with TypeScript dependency, `tsconfig.json`
**Test frameworks:** Jest, Vitest, Mocha, AVA

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `jest.config.ts`, `jest` in devDeps | Jest |
| `vitest.config.ts`, `vitest` in devDeps | Vitest |
| `mocha` in devDeps, `.mocharc.*` | Mocha |
| `ava` in devDeps | AVA |

---

## File Structure

### Naming Conventions
```
*.test.ts     # Most common (Jest, Vitest)
*.spec.ts     # Alternative (common in Angular)
*.test.tsx    # For React component tests
```

### Directory Patterns

**Colocated (recommended for components):**
```
src/
├── components/
│   ├── Button.tsx
│   └── Button.test.tsx
├── services/
│   ├── AuthService.ts
│   └── AuthService.test.ts
```

**Separate directory:**
```
src/
├── components/
│   └── Button.tsx
tests/
├── unit/
│   └── Button.test.tsx
├── integration/
│   └── AuthService.test.ts
```

**__tests__ folders (Jest default):**
```
src/
├── components/
│   ├── __tests__/
│   │   └── Button.test.tsx
│   └── Button.tsx
```

---

## Jest Patterns

### Basic Structure
```typescript
describe('ServiceName', () => {
  // Setup
  beforeEach(() => {
    jest.clearAllMocks();
  });

  // Cleanup
  afterEach(() => {
    jest.restoreAllMocks();
  });

  describe('methodName', () => {
    it('should do expected behavior when condition', () => {
      // Arrange
      const input = 'test';
      
      // Act
      const result = service.method(input);
      
      // Assert
      expect(result).toBe('expected');
    });

    it('should throw when invalid input', () => {
      expect(() => service.method(null)).toThrow('Invalid input');
    });
  });
});
```

### Assertions
```typescript
// Equality
expect(value).toBe(primitive);           // Strict equality (===)
expect(value).toEqual(object);           // Deep equality
expect(value).toStrictEqual(object);     // Deep + type equality

// Truthiness
expect(value).toBeTruthy();
expect(value).toBeFalsy();
expect(value).toBeNull();
expect(value).toBeUndefined();
expect(value).toBeDefined();

// Numbers
expect(value).toBeGreaterThan(3);
expect(value).toBeGreaterThanOrEqual(3);
expect(value).toBeLessThan(5);
expect(value).toBeCloseTo(0.3, 5);       // Floating point

// Strings
expect(string).toMatch(/regex/);
expect(string).toContain('substring');

// Arrays
expect(array).toContain(item);
expect(array).toHaveLength(3);
expect(array).toContainEqual({id: 1});   // Deep equality in array

// Objects
expect(object).toHaveProperty('key');
expect(object).toHaveProperty('key', 'value');
expect(object).toMatchObject({partial: 'match'});

// Exceptions
expect(() => fn()).toThrow();
expect(() => fn()).toThrow('message');
expect(() => fn()).toThrow(ErrorType);

// Async
await expect(promise).resolves.toBe('value');
await expect(promise).rejects.toThrow('error');
```

### Mocking

**Module mocking:**
```typescript
// At top of file
jest.mock('./module');

// With implementation
jest.mock('./module', () => ({
  fn: jest.fn(() => 'mocked'),
}));

// Access mocked module
import { fn } from './module';
const mockedFn = fn as jest.MockedFunction<typeof fn>;
```

**Function mocking:**
```typescript
const mockFn = jest.fn();
mockFn.mockReturnValue('value');
mockFn.mockReturnValueOnce('first call');
mockFn.mockResolvedValue('async value');
mockFn.mockRejectedValue(new Error('failed'));
mockFn.mockImplementation((arg) => arg * 2);

// Assertions
expect(mockFn).toHaveBeenCalled();
expect(mockFn).toHaveBeenCalledTimes(2);
expect(mockFn).toHaveBeenCalledWith('arg1', 'arg2');
expect(mockFn).toHaveBeenLastCalledWith('last');
```

**Spy on existing methods:**
```typescript
const spy = jest.spyOn(object, 'method');
spy.mockReturnValue('mocked');

// Restore original
spy.mockRestore();
```

**Timer mocking:**
```typescript
jest.useFakeTimers();

setTimeout(callback, 1000);
jest.advanceTimersByTime(1000);
expect(callback).toHaveBeenCalled();

jest.useRealTimers();
```

---

## Vitest Patterns

### Basic Structure
```typescript
import { describe, it, expect, beforeEach, vi } from 'vitest';

describe('ServiceName', () => {
  beforeEach(() => {
    vi.clearAllMocks();
  });

  it('should work', () => {
    expect(true).toBe(true);
  });
});
```

### Key Differences from Jest
```typescript
// Mocking
vi.mock('./module');                    // Instead of jest.mock
const mockFn = vi.fn();                 // Instead of jest.fn
vi.spyOn(object, 'method');             // Instead of jest.spyOn
vi.useFakeTimers();                     // Instead of jest.useFakeTimers

// Type-safe mocks
import { MockedFunction } from 'vitest';
```

---

## React Testing (Testing Library)

### Setup
```typescript
import { render, screen, fireEvent, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
```

### Component Testing
```typescript
describe('Button', () => {
  it('renders with text', () => {
    render(<Button>Click me</Button>);
    expect(screen.getByRole('button', { name: /click me/i })).toBeInTheDocument();
  });

  it('calls onClick when clicked', async () => {
    const handleClick = jest.fn();
    const user = userEvent.setup();
    
    render(<Button onClick={handleClick}>Click</Button>);
    await user.click(screen.getByRole('button'));
    
    expect(handleClick).toHaveBeenCalledTimes(1);
  });

  it('shows loading state', () => {
    render(<Button loading>Submit</Button>);
    expect(screen.getByRole('button')).toBeDisabled();
    expect(screen.getByTestId('spinner')).toBeInTheDocument();
  });
});
```

### Queries Priority (most to least preferred)
```typescript
// Accessible to everyone
screen.getByRole('button', { name: /submit/i });
screen.getByLabelText('Email');
screen.getByPlaceholderText('Enter email');
screen.getByText('Hello World');
screen.getByDisplayValue('current value');

// Semantic
screen.getByAltText('Profile picture');
screen.getByTitle('Close');

// Test IDs (last resort)
screen.getByTestId('custom-element');
```

### Async Testing
```typescript
it('loads data', async () => {
  render(<DataComponent />);
  
  // Wait for element to appear
  await waitFor(() => {
    expect(screen.getByText('Loaded')).toBeInTheDocument();
  });
  
  // Or use findBy (combines getBy + waitFor)
  const element = await screen.findByText('Loaded');
  expect(element).toBeInTheDocument();
});
```

### Testing with Providers
```typescript
// test-utils.tsx
const AllProviders = ({ children }: { children: React.ReactNode }) => (
  <ThemeProvider theme={theme}>
    <QueryClientProvider client={queryClient}>
      <Router>
        {children}
      </Router>
    </QueryClientProvider>
  </ThemeProvider>
);

export const renderWithProviders = (ui: React.ReactElement) =>
  render(ui, { wrapper: AllProviders });

// In tests
import { renderWithProviders } from '../test-utils';

it('works with context', () => {
  renderWithProviders(<MyComponent />);
});
```

---

## API Mocking

### MSW (Mock Service Worker)
```typescript
import { rest } from 'msw';
import { setupServer } from 'msw/node';

const server = setupServer(
  rest.get('/api/users', (req, res, ctx) => {
    return res(ctx.json([{ id: 1, name: 'John' }]));
  }),
  
  rest.post('/api/users', async (req, res, ctx) => {
    const body = await req.json();
    return res(ctx.status(201), ctx.json({ id: 2, ...body }));
  })
);

beforeAll(() => server.listen());
afterEach(() => server.resetHandlers());
afterAll(() => server.close());

it('fetches users', async () => {
  render(<UserList />);
  await screen.findByText('John');
});

it('handles error', async () => {
  server.use(
    rest.get('/api/users', (req, res, ctx) => {
      return res(ctx.status(500));
    })
  );
  
  render(<UserList />);
  await screen.findByText('Error loading users');
});
```

---

## Common Patterns

### Testing Hooks
```typescript
import { renderHook, act } from '@testing-library/react';

it('increments counter', () => {
  const { result } = renderHook(() => useCounter());
  
  expect(result.current.count).toBe(0);
  
  act(() => {
    result.current.increment();
  });
  
  expect(result.current.count).toBe(1);
});
```

### Testing Forms
```typescript
it('submits form with values', async () => {
  const onSubmit = jest.fn();
  const user = userEvent.setup();
  
  render(<LoginForm onSubmit={onSubmit} />);
  
  await user.type(screen.getByLabelText('Email'), 'test@example.com');
  await user.type(screen.getByLabelText('Password'), 'password123');
  await user.click(screen.getByRole('button', { name: /submit/i }));
  
  expect(onSubmit).toHaveBeenCalledWith({
    email: 'test@example.com',
    password: 'password123',
  });
});
```

### Testing Error Boundaries
```typescript
it('catches errors', () => {
  const ThrowError = () => {
    throw new Error('Test error');
  };
  
  render(
    <ErrorBoundary fallback={<div>Error occurred</div>}>
      <ThrowError />
    </ErrorBoundary>
  );
  
  expect(screen.getByText('Error occurred')).toBeInTheDocument();
});
```

---

## Configuration

### Jest Config (jest.config.ts)
```typescript
export default {
  preset: 'ts-jest',
  testEnvironment: 'jsdom',
  setupFilesAfterEnv: ['<rootDir>/src/setupTests.ts'],
  moduleNameMapper: {
    '^@/(.*)$': '<rootDir>/src/$1',
  },
  collectCoverageFrom: [
    'src/**/*.{ts,tsx}',
    '!src/**/*.d.ts',
  ],
};
```

### Vitest Config (vitest.config.ts)
```typescript
import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    environment: 'jsdom',
    setupFiles: ['./src/setupTests.ts'],
    globals: true,
  },
});
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
