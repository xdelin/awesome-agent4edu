# JavaScript Testing Ecosystem

Comprehensive reference for testing JavaScript projects.

## Detection

**Manifest files:** `package.json` (no TypeScript)
**Test frameworks:** Jest, Vitest, Mocha, AVA, Jasmine

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `jest.config.js`, `jest` in devDeps | Jest |
| `vitest.config.js`, `vitest` in devDeps | Vitest |
| `mocha` in devDeps, `.mocharc.*` | Mocha + Chai |
| `jasmine` in devDeps | Jasmine |
| `ava` in devDeps | AVA |

---

## File Structure

### Naming Conventions
```
*.test.js     # Most common
*.spec.js     # Alternative
*.test.jsx    # React components
*.test.mjs    # ES Modules
```

### Directory Patterns

**Colocated:**
```
src/
├── utils/
│   ├── format.js
│   └── format.test.js
├── components/
│   ├── Button.jsx
│   └── Button.test.jsx
```

**Separate directory:**
```
src/
├── utils/
│   └── format.js
tests/
├── unit/
│   └── format.test.js
├── integration/
│   └── api.test.js
```

---

## Jest Patterns

### Basic Structure
```javascript
describe('moduleName', () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  describe('functionName', () => {
    it('should return expected value when given valid input', () => {
      // Arrange
      const input = 'test';
      
      // Act
      const result = functionName(input);
      
      // Assert
      expect(result).toBe('expected');
    });

    it('should throw when input is invalid', () => {
      expect(() => functionName(null)).toThrow('Invalid');
    });
  });
});
```

### Assertions
```javascript
// Equality
expect(value).toBe(5);                   // Strict equality
expect(obj).toEqual({ a: 1 });           // Deep equality

// Truthiness
expect(value).toBeTruthy();
expect(value).toBeFalsy();
expect(value).toBeNull();
expect(value).toBeUndefined();
expect(value).toBeDefined();

// Numbers
expect(value).toBeGreaterThan(3);
expect(value).toBeLessThan(10);
expect(0.1 + 0.2).toBeCloseTo(0.3);

// Strings
expect(str).toMatch(/pattern/);
expect(str).toContain('substring');

// Arrays
expect(arr).toContain('item');
expect(arr).toHaveLength(3);
expect(arr).toEqual(expect.arrayContaining([1, 2]));

// Objects
expect(obj).toHaveProperty('key');
expect(obj).toHaveProperty('nested.key', 'value');
expect(obj).toMatchObject({ partial: true });

// Exceptions
expect(() => fn()).toThrow();
expect(() => fn()).toThrow('message');
expect(() => fn()).toThrow(TypeError);

// Promises
await expect(asyncFn()).resolves.toBe('value');
await expect(asyncFn()).rejects.toThrow('error');
```

### Mocking

**Module mocking:**
```javascript
// Hoist to top
jest.mock('./module');

// With factory
jest.mock('./api', () => ({
  fetchData: jest.fn(() => Promise.resolve({ data: 'mocked' })),
}));

// Partial mock
jest.mock('./utils', () => ({
  ...jest.requireActual('./utils'),
  specificFn: jest.fn(),
}));
```

**Function mocking:**
```javascript
const mockFn = jest.fn();

// Return values
mockFn.mockReturnValue('value');
mockFn.mockReturnValueOnce('first').mockReturnValueOnce('second');

// Async
mockFn.mockResolvedValue('async result');
mockFn.mockRejectedValue(new Error('failed'));

// Implementation
mockFn.mockImplementation((x) => x * 2);

// Assertions
expect(mockFn).toHaveBeenCalled();
expect(mockFn).toHaveBeenCalledTimes(2);
expect(mockFn).toHaveBeenCalledWith('arg');
expect(mockFn).toHaveBeenNthCalledWith(1, 'first arg');
```

**Spying:**
```javascript
const spy = jest.spyOn(object, 'method');
spy.mockImplementation(() => 'mocked');

object.method();
expect(spy).toHaveBeenCalled();

spy.mockRestore(); // Restore original
```

**Timer mocking:**
```javascript
jest.useFakeTimers();

const callback = jest.fn();
setTimeout(callback, 1000);

jest.advanceTimersByTime(1000);
expect(callback).toHaveBeenCalled();

// Or run all timers
jest.runAllTimers();

jest.useRealTimers();
```

---

## Mocha + Chai Patterns

### Basic Structure
```javascript
const { expect } = require('chai');

describe('moduleName', function() {
  beforeEach(function() {
    // Setup
  });

  describe('functionName', function() {
    it('should return expected value', function() {
      const result = functionName('input');
      expect(result).to.equal('expected');
    });
  });
});
```

### Chai Assertions
```javascript
// Equality
expect(value).to.equal(5);
expect(obj).to.deep.equal({ a: 1 });
expect(obj).to.eql({ a: 1 });            // Shorthand for deep.equal

// Truthiness
expect(value).to.be.true;
expect(value).to.be.false;
expect(value).to.be.null;
expect(value).to.be.undefined;
expect(value).to.exist;

// Numbers
expect(value).to.be.above(5);
expect(value).to.be.below(10);
expect(value).to.be.within(5, 10);
expect(value).to.be.closeTo(0.3, 0.01);

// Strings
expect(str).to.include('substring');
expect(str).to.match(/pattern/);

// Arrays
expect(arr).to.include('item');
expect(arr).to.have.length(3);
expect(arr).to.include.members([1, 2]);

// Objects
expect(obj).to.have.property('key');
expect(obj).to.have.property('key', 'value');
expect(obj).to.include({ partial: true });

// Exceptions
expect(() => fn()).to.throw();
expect(() => fn()).to.throw('message');
expect(() => fn()).to.throw(TypeError);

// Promises (chai-as-promised)
await expect(promise).to.eventually.equal('value');
await expect(promise).to.be.rejected;
```

### Sinon for Mocking
```javascript
const sinon = require('sinon');

// Stubs
const stub = sinon.stub(object, 'method');
stub.returns('value');
stub.resolves('async value');
stub.rejects(new Error('failed'));
stub.callsFake((arg) => arg * 2);

// Spies
const spy = sinon.spy(object, 'method');
expect(spy.calledOnce).to.be.true;
expect(spy.calledWith('arg')).to.be.true;

// Restore
sinon.restore();
```

---

## React Testing (Testing Library)

### Basic Pattern
```javascript
import { render, screen, fireEvent, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';

describe('Button', () => {
  it('renders correctly', () => {
    render(<Button>Click me</Button>);
    expect(screen.getByRole('button', { name: /click me/i })).toBeInTheDocument();
  });

  it('handles click', async () => {
    const handleClick = jest.fn();
    const user = userEvent.setup();
    
    render(<Button onClick={handleClick}>Click</Button>);
    await user.click(screen.getByRole('button'));
    
    expect(handleClick).toHaveBeenCalled();
  });
});
```

### Query Priority
```javascript
// Most accessible (prefer these)
screen.getByRole('button', { name: /submit/i });
screen.getByLabelText('Email');
screen.getByText('Hello');

// Less accessible
screen.getByTestId('custom-id');
```

### Async Testing
```javascript
it('loads data', async () => {
  render(<DataLoader />);
  
  // Wait for element
  await waitFor(() => {
    expect(screen.getByText('Loaded')).toBeInTheDocument();
  });
  
  // Or use findBy
  const element = await screen.findByText('Loaded');
});
```

---

## Node.js Backend Testing

### Express Route Testing (Supertest)
```javascript
const request = require('supertest');
const app = require('./app');

describe('GET /api/users', () => {
  it('returns users list', async () => {
    const response = await request(app)
      .get('/api/users')
      .expect('Content-Type', /json/)
      .expect(200);
    
    expect(response.body).toEqual(
      expect.arrayContaining([
        expect.objectContaining({ id: expect.any(Number) })
      ])
    );
  });

  it('requires authentication', async () => {
    await request(app)
      .get('/api/users')
      .expect(401);
  });
});

describe('POST /api/users', () => {
  it('creates new user', async () => {
    const response = await request(app)
      .post('/api/users')
      .send({ name: 'John', email: 'john@example.com' })
      .expect(201);
    
    expect(response.body.id).toBeDefined();
  });

  it('validates required fields', async () => {
    const response = await request(app)
      .post('/api/users')
      .send({ name: 'John' }) // Missing email
      .expect(400);
    
    expect(response.body.error).toMatch(/email/i);
  });
});
```

### Database Testing
```javascript
const { db } = require('./database');

describe('UserRepository', () => {
  beforeEach(async () => {
    await db.migrate.latest();
    await db.seed.run();
  });

  afterEach(async () => {
    await db.migrate.rollback();
  });

  afterAll(async () => {
    await db.destroy();
  });

  it('finds user by id', async () => {
    const user = await userRepository.findById(1);
    expect(user.name).toBe('Test User');
  });
});
```

---

## ES Modules Testing

### Jest with ESM
```javascript
// jest.config.js
export default {
  transform: {},
  testEnvironment: 'node',
  moduleFileExtensions: ['js', 'mjs'],
};

// Test file (*.test.mjs)
import { functionToTest } from './module.mjs';

describe('ESM module', () => {
  it('works', () => {
    expect(functionToTest()).toBe('result');
  });
});
```

---

## Configuration

### Jest Config (jest.config.js)
```javascript
module.exports = {
  testEnvironment: 'node', // or 'jsdom' for browser
  setupFilesAfterEnv: ['./setupTests.js'],
  collectCoverageFrom: [
    'src/**/*.js',
    '!src/**/*.test.js',
  ],
  coverageThreshold: {
    global: {
      branches: 80,
      functions: 80,
      lines: 80,
    },
  },
};
```

### Mocha Config (.mocharc.json)
```json
{
  "spec": "test/**/*.test.js",
  "timeout": 5000,
  "recursive": true,
  "require": ["chai/register-expect"]
}
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
