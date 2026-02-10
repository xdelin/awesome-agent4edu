# Rust Testing Ecosystem

Comprehensive reference for testing Rust projects.

## Detection

**Manifest files:** `Cargo.toml`, `Cargo.lock`
**Test frameworks:** Built-in `#[test]`, proptest, criterion (benchmarks)

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `Cargo.toml` exists | Rust (built-in testing) |
| `proptest` in dependencies | Property-based testing |
| `criterion` in dev-dependencies | Benchmarking |
| `mockall` in dev-dependencies | Mocking |
| `rstest` in dev-dependencies | Fixtures/Parameterized tests |

---

## File Structure

### Naming Conventions
```
*_test.rs            # Integration tests in tests/
tests.rs             # Inline test module
mod tests { }        # Test module within source file
```

### Directory Patterns

**Standard structure:**
```
my_project/
├── Cargo.toml
├── src/
│   ├── lib.rs              # Library code + unit tests
│   ├── main.rs             # Binary entry point
│   ├── user.rs             # Module with inline tests
│   └── service/
│       ├── mod.rs
│       └── payment.rs      # Inline tests at bottom
└── tests/                  # Integration tests
    ├── common/
    │   └── mod.rs          # Shared test utilities
    ├── user_test.rs
    └── api_test.rs
```

### Inline Unit Tests
```rust
// src/calculator.rs
pub fn add(a: i32, b: i32) -> i32 {
    a + b
}

pub fn divide(a: i32, b: i32) -> Result<i32, String> {
    if b == 0 {
        return Err("Division by zero".to_string());
    }
    Ok(a / b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        assert_eq!(add(2, 3), 5);
    }

    #[test]
    fn test_add_negative() {
        assert_eq!(add(-1, -2), -3);
    }

    #[test]
    fn test_divide() {
        assert_eq!(divide(10, 2).unwrap(), 5);
    }

    #[test]
    fn test_divide_by_zero() {
        assert!(divide(10, 0).is_err());
    }
}
```

---

## Built-in Testing Patterns

### Basic Tests
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        // Arrange
        let input = "hello";
        
        // Act
        let result = process(input);
        
        // Assert
        assert_eq!(result, "HELLO");
    }

    #[test]
    fn test_with_message() {
        let result = calculate(5);
        assert_eq!(result, 25, "Expected 5² = 25, got {}", result);
    }

    #[test]
    #[should_panic]
    fn test_panic() {
        panic_function();
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn test_panic_message() {
        let v = vec![1, 2, 3];
        v[99]; // Panics
    }

    #[test]
    fn test_result() -> Result<(), String> {
        let result = fallible_function()?;
        assert_eq!(result, expected);
        Ok(())
    }

    #[test]
    #[ignore]
    fn expensive_test() {
        // Run with: cargo test -- --ignored
    }
}
```

### Assertions
```rust
// Equality
assert_eq!(actual, expected);
assert_eq!(actual, expected, "Custom message: {}", info);
assert_ne!(actual, unexpected);

// Boolean
assert!(condition);
assert!(condition, "Failed because: {}", reason);
assert!(!condition);

// Pattern matching
assert!(matches!(value, Some(_)));
assert!(matches!(result, Ok(x) if x > 0));
assert!(matches!(error, Err(MyError::NotFound)));

// Debug assertions (only in debug builds)
debug_assert!(condition);
debug_assert_eq!(a, b);
```

### Testing Results and Options
```rust
#[test]
fn test_result_ok() {
    let result: Result<i32, &str> = Ok(42);
    
    assert!(result.is_ok());
    assert_eq!(result.unwrap(), 42);
    assert_eq!(result.ok(), Some(42));
}

#[test]
fn test_result_err() {
    let result: Result<i32, &str> = Err("error");
    
    assert!(result.is_err());
    assert_eq!(result.unwrap_err(), "error");
    assert_eq!(result.err(), Some("error"));
}

#[test]
fn test_option() {
    let some_value: Option<i32> = Some(42);
    let none_value: Option<i32> = None;
    
    assert!(some_value.is_some());
    assert!(none_value.is_none());
    assert_eq!(some_value.unwrap(), 42);
}

// Using ? operator in tests
#[test]
fn test_with_question_mark() -> Result<(), Box<dyn std::error::Error>> {
    let data = read_file("test.txt")?;
    let parsed = parse_data(&data)?;
    assert_eq!(parsed.len(), 10);
    Ok(())
}
```

### Testing Panics
```rust
use std::panic;

#[test]
fn test_catch_panic() {
    let result = panic::catch_unwind(|| {
        panic!("test panic");
    });
    assert!(result.is_err());
}

#[test]
fn test_panic_payload() {
    let result = panic::catch_unwind(|| {
        panic!("specific message");
    });
    
    if let Err(err) = result {
        if let Some(s) = err.downcast_ref::<&str>() {
            assert_eq!(*s, "specific message");
        }
    }
}
```

---

## Integration Tests

```rust
// tests/integration_test.rs
use my_crate::{Config, Server};

mod common;

#[test]
fn test_full_workflow() {
    // Setup
    let config = common::setup_test_config();
    let server = Server::new(config);
    
    // Test
    let result = server.process_request("test");
    
    // Assert
    assert!(result.is_ok());
    assert_eq!(result.unwrap().status, 200);
}

// tests/common/mod.rs
pub fn setup_test_config() -> Config {
    Config {
        debug: true,
        port: 0, // Random port
        ..Default::default()
    }
}

pub fn create_test_user() -> User {
    User {
        id: 1,
        name: "Test User".to_string(),
    }
}
```

---

## rstest (Fixtures & Parameterized)

### Fixtures
```rust
use rstest::*;

struct Database {
    // ...
}

#[fixture]
fn database() -> Database {
    Database::new_test_instance()
}

#[fixture]
fn user(database: Database) -> User {
    database.create_user("test@example.com")
}

#[rstest]
fn test_user_creation(user: User) {
    assert_eq!(user.email, "test@example.com");
}

#[rstest]
fn test_with_database(database: Database) {
    assert!(database.is_connected());
}
```

### Parameterized Tests
```rust
use rstest::*;

#[rstest]
#[case(0, 0, 0)]
#[case(1, 2, 3)]
#[case(-1, 1, 0)]
#[case(100, 200, 300)]
fn test_add(#[case] a: i32, #[case] b: i32, #[case] expected: i32) {
    assert_eq!(add(a, b), expected);
}

#[rstest]
#[case("", false)]
#[case("test@example.com", true)]
#[case("invalid", false)]
fn test_email_validation(#[case] input: &str, #[case] expected: bool) {
    assert_eq!(is_valid_email(input), expected);
}

// Matrix testing
#[rstest]
fn test_matrix(
    #[values(1, 2, 3)] x: i32,
    #[values(10, 20)] y: i32,
) {
    assert!(multiply(x, y) > 0);
}
```

---

## proptest (Property-Based Testing)

```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn test_add_commutative(a: i32, b: i32) {
        prop_assert_eq!(add(a, b), add(b, a));
    }

    #[test]
    fn test_add_associative(a: i32, b: i32, c: i32) {
        prop_assert_eq!(add(add(a, b), c), add(a, add(b, c)));
    }

    #[test]
    fn test_string_reverse(s: String) {
        let reversed: String = s.chars().rev().collect();
        let double_reversed: String = reversed.chars().rev().collect();
        prop_assert_eq!(s, double_reversed);
    }

    #[test]
    fn test_parse_roundtrip(value in 0i32..1000) {
        let s = value.to_string();
        let parsed: i32 = s.parse().unwrap();
        prop_assert_eq!(value, parsed);
    }
}

// Custom strategies
fn valid_email_strategy() -> impl Strategy<Value = String> {
    "[a-z]{1,10}@[a-z]{1,5}\\.(com|org|net)"
}

proptest! {
    #[test]
    fn test_email_validation(email in valid_email_strategy()) {
        prop_assert!(is_valid_email(&email));
    }
}

// Constrained values
proptest! {
    #[test]
    fn test_with_constraints(
        age in 18..65i32,
        name in "[A-Z][a-z]{2,10}",
    ) {
        let user = User::new(&name, age);
        prop_assert!(user.is_valid());
    }
}
```

---

## Mocking with mockall

```rust
use mockall::*;
use mockall::predicate::*;

// Define trait
#[automock]
trait Database {
    fn get_user(&self, id: i32) -> Option<User>;
    fn save_user(&self, user: &User) -> Result<(), DbError>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_user_service() {
        let mut mock_db = MockDatabase::new();
        
        // Setup expectations
        mock_db.expect_get_user()
            .with(eq(1))
            .times(1)
            .returning(|_| Some(User { id: 1, name: "John".into() }));
        
        let service = UserService::new(Box::new(mock_db));
        let user = service.find_user(1);
        
        assert!(user.is_some());
        assert_eq!(user.unwrap().name, "John");
    }

    #[test]
    fn test_save_user() {
        let mut mock_db = MockDatabase::new();
        
        mock_db.expect_save_user()
            .withf(|user| user.name == "Jane")
            .times(1)
            .returning(|_| Ok(()));
        
        let service = UserService::new(Box::new(mock_db));
        let result = service.create_user("Jane");
        
        assert!(result.is_ok());
    }

    #[test]
    fn test_sequence() {
        let mut mock_db = MockDatabase::new();
        let mut seq = Sequence::new();
        
        mock_db.expect_get_user()
            .times(1)
            .in_sequence(&mut seq)
            .returning(|_| None);
        
        mock_db.expect_save_user()
            .times(1)
            .in_sequence(&mut seq)
            .returning(|_| Ok(()));
        
        // ... test code that calls get then save
    }
}
```

---

## Async Testing

```rust
// With tokio
#[tokio::test]
async fn test_async_function() {
    let result = async_fetch_data().await;
    assert!(result.is_ok());
}

#[tokio::test]
async fn test_async_timeout() {
    let result = tokio::time::timeout(
        std::time::Duration::from_secs(5),
        slow_operation()
    ).await;
    
    assert!(result.is_ok());
}

// With async-std
#[async_std::test]
async fn test_with_async_std() {
    let result = fetch().await;
    assert_eq!(result, expected);
}

// Testing async errors
#[tokio::test]
async fn test_async_error() {
    let result = fallible_async().await;
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), MyError::NotFound));
}
```

---

## Test Organization

### Test Modules
```rust
// src/lib.rs
pub mod user;
pub mod service;

#[cfg(test)]
mod tests {
    use super::*;

    mod user_tests {
        use super::*;

        #[test]
        fn test_user_creation() {
            // ...
        }
    }

    mod service_tests {
        use super::*;

        #[test]
        fn test_service_init() {
            // ...
        }
    }
}
```

### Conditional Compilation
```rust
#[cfg(test)]
mod tests {
    use super::*;

    // Test-only helpers
    fn setup_test_env() -> TestEnv {
        // ...
    }

    // Test-only trait implementations
    impl TestHelpers for MyStruct {
        fn test_method(&self) -> bool {
            true
        }
    }
}

// Test-only dependencies
#[cfg(test)]
use test_utilities::*;
```

---

## Benchmarking (criterion)

```rust
// benches/my_benchmark.rs
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use my_crate::expensive_function;

fn benchmark_expensive_function(c: &mut Criterion) {
    c.bench_function("expensive_function", |b| {
        b.iter(|| expensive_function(black_box(1000)))
    });
}

fn benchmark_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("Algorithms");
    
    group.bench_function("algorithm_a", |b| {
        b.iter(|| algorithm_a(black_box(&data)))
    });
    
    group.bench_function("algorithm_b", |b| {
        b.iter(|| algorithm_b(black_box(&data)))
    });
    
    group.finish();
}

criterion_group!(benches, benchmark_expensive_function, benchmark_comparison);
criterion_main!(benches);
```

---

## Configuration

### Cargo.toml
```toml
[package]
name = "my_project"
version = "0.1.0"
edition = "2021"

[dependencies]
# ...

[dev-dependencies]
rstest = "0.18"
proptest = "1.4"
mockall = "0.12"
criterion = "0.5"
tokio = { version = "1", features = ["test-util", "macros", "rt-multi-thread"] }

[[bench]]
name = "my_benchmark"
harness = false
```

### Running Tests
```bash
# Run all tests
cargo test

# Run specific test
cargo test test_name

# Run tests in specific module
cargo test module_name::

# Run ignored tests
cargo test -- --ignored

# Run tests with output
cargo test -- --nocapture

# Run benchmarks
cargo bench

# Run with features
cargo test --features "feature_name"

# Run only doc tests
cargo test --doc

# Run only unit tests
cargo test --lib

# Run only integration tests
cargo test --test integration_test
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
